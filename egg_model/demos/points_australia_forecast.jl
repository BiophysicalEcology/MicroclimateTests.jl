# Stage 8, scaled: history -> forecast splice over a land-masked regular grid
# of points (points_australia.jl's Pass A setup) and a user-chosen number of
# ACCESS-S2 ensemble members. Microclimate results are cached to disk per
# (leg, member, batch) -- raw multi-depth, Float32-trimmed to
# CACHED_DEPTH_RANGE, batched like points_australia.jl since a full grid
# won't fit one MicroVectorProblem call. The egg-model runs share a
# points_australia-style cache pool (one integrator per thread, reinit! per
# point/member work item).

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_properties_model, example_soil_hydraulic_model
using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Rasters, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using Dates, Unitful, Statistics
using Serialization
using NaturalEarth, GeoInterface

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))
include(joinpath(@__DIR__, "..", "src", "access_s2.jl"))

include(joinpath(@__DIR__, "..", "params", "chortoicetes.jl"))

output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

n_ensembles = 25    # how many of the up-to-99 ACCESS-S2 members to run
# no trajectories needed -- outputs are maps + a stats CSV, not per-point
# development/mass/temperature curves.
save_trajectory = false
diapause = true
oviposition_date = Date(2026, 4, 25)
use_cache = true
batch_size = 100

if diapause
    nest_depth = 5.0u"cm"
    cold_hour_threshold = 30u"d"
    diapause_hour_threshold = 240.0u"d"
else
    nest_depth = 10.0u"cm"
    cold_hour_threshold = 0.0u"d"
    diapause_hour_threshold = 0.0u"d"
end

# raw caches only keep depths down to this value -- set past whatever
# nest_depth values you want to compare (deeper needs a fresh solve).
cache_max_depth = 15.0u"cm"
CACHED_DEPTH_RANGE = 1:nearest_node(cache_max_depth, depths)
nest_node = nearest_node(nest_depth, depths)
nest_node in CACHED_DEPTH_RANGE || error(
    "nest_depth=$nest_depth is deeper than cache_max_depth=$cache_max_depth -- raise cache_max_depth and re-solve")

# Just what egg_nest_forcing (forcing.jl) actually reads, plus soil_moisture
# (kept for a possible future moisture threshold).
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
)

issue_date = Date(2026, 7, 1)   # ACCESS-S2 issue date -- also used below for the ACCESS-S2 land-mask probe

# ── regular grid over south-eastern Australia, land-masked (points_australia.jl) ──

lon_range = range(135.0, 153.0; length=50)
lat_range = range(-39.0, -24.0; length=40)
all_grid_points = vec([(lon, lat) for lon in lon_range, lat in lat_range])   # column-major (lon fastest)

# TODO work out what DEM SILO uses and use that instead, probably only need
# that and not the CRUCL2_ELV filter.
const CRUCL2_ELV = read(Raster(RasterDataSources.getraster(CRUCL2); name=:elv, lazy=true))
has_crucl2_land(lon, lat) = !ismissing(CRUCL2_ELV[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_crucl2_land(p...), all_grid_points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the CRUCL2 land-mask pre-check.")

const SILO_MAXTEMP_PROBE = read(Raster(RasterDataSources.getraster(SILO, :max_temp; date=Date(2020, 1, 1)); name=:max_temp, lazy=true)[Ti(1)])
has_silo_land(lon, lat) = !ismissing(SILO_MAXTEMP_PROBE[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_silo_land(p...), points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the SILO land-mask pre-check.")

# ACCESS-S2's own land/sea representation doesn't always agree with
# CRUCL2/SILO's -- a point that passes both those masks can still be
# `missing` in ACCESS-S2 at its (coarser) resolution, which propagates into
# the soil ODE as degenerate forcing rather than a clean error.
const ACCESS_S2_TMAX_PROBE = read(Raster(joinpath(ENV["RASTERDATASOURCES_PATH"], "ACCESS-S2", "$(Dates.format(issue_date, "yyyymmdd"))_tmax.nc");
    name=:tmax, lazy=true)[Dim{:ensemble}(Rasters.At(1)), Ti(1)])
has_access_s2_land(lon, lat) = !ismissing(ACCESS_S2_TMAX_PROBE[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_access_s2_land(p...), points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the ACCESS-S2 land-mask pre-check.")

# sanity check before committing to the full (expensive) solve below
using Plots
let
    domain = Extent(X=(minimum(lon_range) - 1, maximum(lon_range) + 1),
                     Y=(minimum(lat_range) - 1, maximum(lat_range) + 1))
    elv_domain = crop(CRUCL2_ELV; to=domain, touches=true)
    p = plot(elv_domain; title="CRUCL2 elevation + grid points (green=kept, red=rejected)")
    rejected = setdiff(all_grid_points, points)
    scatter!(p, first.(points), last.(points); markersize=2, markerstrokewidth=0, color=:green, label="kept")
    scatter!(p, first.(rejected), last.(rejected); markersize=2, markerstrokewidth=0, color=:red, label="rejected")
    savefig(p, joinpath(output_dir, "history_forecast_splice_domain_check.png"))
    display(p)
end

# ── uniform soil texture (Campbell & Norman) -- Pass A, no per-point texture ──

const CAMPBELL_NORMAN_TEXTURES = (
    sand             = (air_entry=0.7u"J/kg", b=1.7, Ksat=5.8e-3u"kg*s/m^3", field_capacity=0.09, wilting_point=0.03),
    loamy_sand       = (air_entry=0.9u"J/kg", b=2.1, Ksat=1.7e-3u"kg*s/m^3", field_capacity=0.13, wilting_point=0.06),
    sandy_loam       = (air_entry=1.5u"J/kg", b=3.1, Ksat=7.2e-4u"kg*s/m^3", field_capacity=0.21, wilting_point=0.10),
    loam             = (air_entry=1.1u"J/kg", b=4.5, Ksat=3.7e-4u"kg*s/m^3", field_capacity=0.27, wilting_point=0.12),
    silt_loam        = (air_entry=2.1u"J/kg", b=4.7, Ksat=1.9e-4u"kg*s/m^3", field_capacity=0.33, wilting_point=0.13),
    sandy_clay_loam  = (air_entry=2.8u"J/kg", b=4.0, Ksat=1.2e-3u"kg*s/m^3", field_capacity=0.26, wilting_point=0.15),
    clay_loam        = (air_entry=2.6u"J/kg", b=5.2, Ksat=6.4e-5u"kg*s/m^3", field_capacity=0.32, wilting_point=0.20),
    silty_clay_loam  = (air_entry=3.3u"J/kg", b=6.6, Ksat=4.2e-5u"kg*s/m^3", field_capacity=0.37, wilting_point=0.32),
    sandy_clay       = (air_entry=2.9u"J/kg", b=6.0, Ksat=3.3e-5u"kg*s/m^3", field_capacity=0.34, wilting_point=0.24),
    silty_clay       = (air_entry=3.4u"J/kg", b=7.9, Ksat=2.5e-5u"kg*s/m^3", field_capacity=0.39, wilting_point=0.25),
    clay             = (air_entry=3.7u"J/kg", b=7.6, Ksat=1.7e-5u"kg*s/m^3", field_capacity=0.40, wilting_point=0.27),
)
soil_source = :sandy_loam

function soil_profile_from_texture(texture::NamedTuple, depths;
    bulk_density=1.3u"Mg/m^3", mineral_density=2.560u"Mg/m^3",
    mineral_conductivity=1.25u"W/m/K", mineral_heat_capacity=870.0u"J/kg/K",
    root_density=Microclimate.example_campbell_hydraulic_profile(depths).root_density,
)
    n = length(depths)
    Microclimate.SoilProfile(;
        bulk_density=fill(bulk_density, n), mineral_density=fill(mineral_density, n),
        mineral_conductivity=fill(mineral_conductivity, n), mineral_heat_capacity=fill(mineral_heat_capacity, n),
        hydraulics=Microclimate.CampbellHydraulicProfile(;
            air_entry_water_potential=fill(-texture.air_entry, n),
            saturated_hydraulic_conductivity=fill(texture.Ksat, n),
            campbell_b_parameter=fill(texture.b, n), root_density,
        ),
    )
end
soil_profile = soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[soil_source], depths)
soil_hydraulics = (;
    air_entry_potential    = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b             = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)
environment_pars = example_environment_pars()

# ── egg model, identical config to point_silo_deterministic.jl ──

geometry = Ellipsoid(initial_egg_mass, egg_density, axis_ratio, axis_ratio)
arrest = ProportionWindowArrest(;
    cold_temperature, diapause_window, quiescence_windows,
    cold_hour_threshold, diapause_hour_threshold,
    desiccation_tolerance,
)
dm = arrhenius_development_model(;
    T_A, T_AL, T_AH, T_L, T_H, T_ref, rate_at_reference, rate_unit,
)
stage = SteppedHydricStage(;
    conductance_threshold, wetness_threshold, dormant_conductance, active_conductance,
    dormant_wetness, active_wetness,
)
pars = EggParameters(;
    hydraulic_conductance, specific_hydration, conduction_fraction, skin_wetness,
    initial_egg_mass, minimum_egg_mass,
)
survival_model = CombinedSurvival(
    HardTemperatureLimit(; lower_lethal_temperature=u"K"(-5.0u"°C"), upper_lethal_temperature=u"K"(52.0u"°C")),
    DesiccationLimit(; dry_mass=0.1 * initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

historical_dates = oviposition_date:Day(1):issue_date
members = 1:n_ensembles
n = length(points)
n_batches = cld(n, batch_size)

_raw_layers(output, i) = (;
    soil_temperature          = Float32.(collect(output.soil_temperature[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_moisture             = Float32.(collect(output.soil_moisture[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_water_potential      = Float32.(collect(output.soil_water_potential[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_thermal_conductivity = Float32.(collect(output.soil_thermal_conductivity[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_humidity             = Float32.(collect(output.soil_humidity[point=i, depth=CACHED_DEPTH_RANGE])),
)

# solves `model` over `points` in batches, caching each batch to disk under
# `label` -- shared by the historical leg and every forecast member.
function solve_batched(model, label, points, dates, init)
    n = length(points)
    n_batches = cld(n, batch_size)
    per_point = Vector{Any}(undef, n)
    final_soil_temperature = Vector{Any}(undef, n)
    final_soil_moisture = Vector{Any}(undef, n)
    for b in 1:n_batches
        i_start, i_end = (b - 1) * batch_size + 1, min(b * batch_size, n)
        batch_points = points[i_start:i_end]
        batch_cache_path = joinpath(output_dir, "$(label)_batch$(b)of$(n_batches)_n$(n).jls")
        batch = if isfile(batch_cache_path) && use_cache
            println("Loading cached $label batch $b/$n_batches ($(length(batch_points)) points)...")
            deserialize(batch_cache_path)
        else
            println("Solving $label batch $b/$n_batches ($(length(batch_points)) points)...")
            batch_problem = MicroVectorProblem(; model, points=batch_points, dates, soil_profile, init)
            @time batch_output = solve(batch_problem)
            last_hour = size(batch_output.soil_temperature, 1)
            result = (;
                per_point = [_raw_layers(batch_output, i) for i in 1:length(batch_points)],
                # full (untrimmed) depth profile at the last hour, for splice continuity --
                # kept Float64 (unlike per_point's Float32 trim): this seeds
                # MicroVectorProblem's init=, which expects Float64 internally.
                final_soil_temperature = [collect(batch_output.soil_temperature[point=i, Ti=last_hour]) for i in 1:length(batch_points)],
                final_soil_moisture    = [collect(batch_output.soil_moisture[point=i, Ti=last_hour]) for i in 1:length(batch_points)],
            )
            serialize(batch_cache_path, result)
            batch_output = nothing
            batch_problem = nothing
            GC.gc()
            result
        end
        per_point[i_start:i_end] .= batch.per_point
        final_soil_temperature[i_start:i_end] .= batch.final_soil_temperature
        final_soil_moisture[i_start:i_end] .= batch.final_soil_moisture
    end
    (; per_point, final_soil_temperature, final_soil_moisture)
end

# ── historical microclimate: SILO, oviposition_date -> issue_date ──

historical_model = MicroMapModel(;
    micro_model=MicroModel(;
        depths, heights,
        soil_properties_model=example_soil_properties_model(),
        soil_hydraulic_model=example_soil_hydraulic_model(),
        snow_model=NoSnow(),
        config=MicroConfig(soil_moisture_strategy=DynamicSoilMoisture()),
    ),
    dem_source=CRUCL2, weather_source=SILO,
    surface_albedo_source=0.15, roughness_height_source=0.004u"m",
    compute_terrain=false, output_layers,
)
println("Solving historical SILO microclimate for $n points: $oviposition_date to $issue_date...")
historical_raw = solve_batched(historical_model, "splice_historical_lay$(oviposition_date)_issue$(issue_date)",
    points, historical_dates, (; soil_moisture=fill(0.2, length(depths))))

historical_day_range = 1:size(historical_raw.per_point[1].soil_temperature, 1)
historical_tspan = (0.0u"hr", length(historical_day_range) * 1.0u"hr")
historical_forcings = [egg_nest_forcing(historical_raw.per_point[i], historical_day_range, nest_node, environment_pars)
                        for i in 1:n]

# MicroVectorProblem takes one shared init.soil_moisture/soil_temperature
# vector for all points, not a point x depth matrix -- average the per-point
# historical endpoints across points (a known approximation, not per-point
# continuity) to seed the forecast leg.
now_soil_temperature = reduce(+, historical_raw.final_soil_temperature) ./ n
now_soil_moisture = reduce(+, historical_raw.final_soil_moisture) ./ n

# ── egg-model cache pool: one integrator per thread, reinit! per work item ──

nworkers = min(Threads.nthreads(), n)
build_cache() = init_egg_cache(egg_model, pars, initial_state, soil_hydraulics,
    historical_forcings[1], historical_tspan; save_trajectory)
cache_pool = Channel{typeof(build_cache())}(nworkers)
for _ in 1:nworkers
    put!(cache_pool, build_cache())
end

function run_threaded!(results, indices, f)
    work = Channel{Int}(length(indices))
    for i in indices
        put!(work, i)
    end
    close(work)
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            cache = take!(cache_pool)
            for i in work
                results[i] = f(cache, i)
            end
            put!(cache_pool, cache)
        end
    end
end

println("Running historical egg model at $n points, lay date $oviposition_date...")
historical_egg_results = Vector{Any}(undef, n)
@time run_threaded!(historical_egg_results, 1:n, (cache, i) ->
    simulate_egg!(cache, initial_state, soil_hydraulics, historical_forcings[i], historical_tspan))

n_hatched_historically = count(r -> r.hatched, historical_egg_results)
n_died_historically = count(r -> r.died, historical_egg_results)
println("$n_hatched_historically/$n hatched and $n_died_historically/$n died before $issue_date; " *
        "$(n - n_hatched_historically - n_died_historically)/$n continue into the forecast.")
forecast_point_indices = [i for i in 1:n
                           if !historical_egg_results[i].hatched && !historical_egg_results[i].died]

# ── forecast leg: one ACCESS-S2 member at a time (all points, batched), egg
# model threaded across points within each member ──

forecast_horizon_days = 214
forecast_dates = issue_date:Day(1):(issue_date + Day(forecast_horizon_days - 1))
forecast_tspan = (0.0u"hr", length(forecast_dates) * 24.0u"hr")
forecast_outcomes = Matrix{Any}(nothing, n, length(members))

if !isempty(forecast_point_indices)
    for (member_index, member) in pairs(members)
        println("\nACCESS-S2 member $member/$(length(members)): solving microclimate for $n points...")
        forecast_model = MicroMapModel(;
            micro_model=MicroModel(;
                depths, heights,
                soil_properties_model=example_soil_properties_model(),
                soil_hydraulic_model=example_soil_hydraulic_model(),
                snow_model=NoSnow(),
                config=MicroConfig(
                    soil_moisture_strategy=DynamicSoilMoisture(),
                    convergence=FixedSoilTemperatureIterations(1),
                ),
            ),
            dem_source=CRUCL2, weather_source=AccessS2(issue_date, member),
            surface_albedo_source=0.15, roughness_height_source=0.004u"m",
            compute_terrain=false, output_layers,
        )
        forecast_raw = solve_batched(forecast_model, "splice_forecast_member$(member)_issue$(issue_date)",
            points, forecast_dates, (; soil_moisture=now_soil_moisture, soil_temperature=now_soil_temperature))

        forecast_day_range = 1:size(forecast_raw.per_point[1].soil_temperature, 1)
        forecast_forcings = Dict(i => egg_nest_forcing(forecast_raw.per_point[i], forecast_day_range, nest_node, environment_pars)
                                  for i in forecast_point_indices)

        member_outcomes = Vector{Any}(undef, n)
        run_threaded!(member_outcomes, forecast_point_indices, (cache, i) ->
            simulate_egg!(cache, historical_egg_results[i].final_state, soil_hydraulics,
                forecast_forcings[i], forecast_tspan))
        forecast_outcomes[forecast_point_indices, member_index] .= member_outcomes[forecast_point_indices]
    end
end

# ── per-point summary: median/earliest/latest hatching member, death-cause
# counts, saved to one combined CSV ──

function summarize_point(historical_egg_result, point_outcomes, members)
    if historical_egg_result.hatched
        return (; status=:hatched_historically, counts=Dict(:hatched => 1),
            earliest=nothing, median=nothing, latest=nothing)
    elseif historical_egg_result.died
        return (; status=Symbol("died_historically_$(historical_egg_result.death_cause)"),
            counts=Dict(historical_egg_result.death_cause => 1),
            earliest=nothing, median=nothing, latest=nothing)
    end
    counts = Dict{Symbol,Int}()
    for r in point_outcomes
        code = r.hatched ? :hatched : r.died ? r.death_cause : :timeout
        counts[code] = get(counts, code, 0) + 1
    end
    hatched = [(member, r) for (member, r) in zip(members, point_outcomes) if r.hatched]
    if isempty(hatched)
        return (; status=:no_hatch, counts, earliest=nothing, median=nothing, latest=nothing)
    end
    sorted = sort(hatched; by=x -> x[2].hatch_time)
    (; status=:forecast, counts,
        earliest=sorted[1], median=sorted[cld(length(sorted), 2)], latest=sorted[end])
end

hatch_date_of(member_result, issue_date) =
    issue_date + Day(round(Int, ustrip(u"d", member_result[2].hatch_time)))
mass_mg_of(member_result) = ustrip(u"mg", member_result[2].final_state.egg_mass)

point_summaries = [summarize_point(historical_egg_results[i],
    [forecast_outcomes[i, m] for m in eachindex(members) if forecast_outcomes[i, m] !== nothing], members)
    for i in 1:n]

stats_path = joinpath(output_dir, "history_forecast_splice_stats.csv")
open(stats_path, "w") do io
    println(io, "lon,lat,status,n_members,n_hatched,n_died_cold,n_died_heat,n_died_desiccation,n_timeout," *
                 "earliest_hatch_date,earliest_hatch_mass_mg,median_hatch_date,median_hatch_mass_mg," *
                 "latest_hatch_date,latest_hatch_mass_mg")
    for i in 1:n
        s = point_summaries[i]
        n_run = s.earliest === nothing && s.status != :no_hatch ? 0 : length(members)
        c(sym) = get(s.counts, sym, 0)
        fields = if s.earliest === nothing
            ("", "", "", "", "", "")
        else
            (string(hatch_date_of(s.earliest, issue_date)), string(round(mass_mg_of(s.earliest); digits=2)),
             string(hatch_date_of(s.median, issue_date)), string(round(mass_mg_of(s.median); digits=2)),
             string(hatch_date_of(s.latest, issue_date)), string(round(mass_mg_of(s.latest); digits=2)))
        end
        println(io, join((points[i][1], points[i][2], s.status, n_run, c(:hatched), c(:cold), c(:heat),
            c(:desiccation), c(:timeout), fields...), ","))
    end
end
println("Saved $stats_path")

# ── maps: heatmap over the regular grid + basemap (points_australia.jl) ──

point_to_index = Dict(p => i for (i, p) in enumerate(all_grid_points))
to_heatmap_z(vals) = permutedims(reshape(vals, length(lon_range), length(lat_range)))

const GI = GeoInterface
function _collect_boundary!(xs, ys, geom)
    trait = GI.geomtrait(geom)
    if trait isa GI.AbstractPointTrait
        push!(xs, GI.x(geom)); push!(ys, GI.y(geom))
    elseif trait isa GI.LinearRingTrait || trait isa GI.LineStringTrait
        for pt in GI.getpoint(geom)
            push!(xs, GI.x(pt)); push!(ys, GI.y(pt))
        end
        push!(xs, NaN); push!(ys, NaN)
    else
        for sub in GI.getgeom(geom)
            _collect_boundary!(xs, ys, sub)
        end
    end
end
au_states = naturalearth("admin_1_states_provinces", 10)
state_xs, state_ys = Float64[], Float64[]
for i in 1:length(au_states)
    au_states[i].iso_a2 == "AU" || continue
    _collect_boundary!(state_xs, state_ys, au_states[i].geometry)
end

map_towns = ["Birdsville, Queensland", "Roma, Queensland", "Charleville, Queensland",
             "Dubbo, New South Wales", "Broken Hill, New South Wales", "Bourke, New South Wales",
             "Mildura, Victoria", "Canberra, Australia"]
towns_cache_path = joinpath(output_dir, "points_australia_towns.jls")
towns = if isfile(towns_cache_path)
    deserialize(towns_cache_path)
else
    t = map(name -> geocode(name), map_towns)
    serialize(towns_cache_path, t)
    t
end

function add_basemap!(p)
    plot!(p, state_xs, state_ys; color=:black, linewidth=0.75, label=nothing)
    scatter!(p, [t.lon for t in towns], [t.lat for t in towns]; color=:black, markersize=2, label=nothing)
    for t in towns
        annotate!(p, t.lon, t.lat, text(first(split(t.display_name, ",")), 6, :left, :bottom))
    end
    xlims!(p, extrema(lon_range)...)
    ylims!(p, extrema(lat_range)...)
    p
end

historical_hatch_date(r) = oviposition_date + Day(round(Int, ustrip(u"d", r.hatch_time)))

median_dates = fill(NaN, length(all_grid_points))
earliest_dates = fill(NaN, length(all_grid_points))
latest_dates = fill(NaN, length(all_grid_points))
median_mass = fill(NaN, length(all_grid_points))
earliest_mass = fill(NaN, length(all_grid_points))
latest_mass = fill(NaN, length(all_grid_points))
frac_hatched = fill(NaN, length(all_grid_points))
frac_cold = fill(NaN, length(all_grid_points))
frac_heat = fill(NaN, length(all_grid_points))
frac_desiccation = fill(NaN, length(all_grid_points))
for i in 1:n
    gi = point_to_index[points[i]]
    s = point_summaries[i]
    if s.status == :hatched_historically
        d = Dates.value(historical_hatch_date(historical_egg_results[i]))
        median_dates[gi] = earliest_dates[gi] = latest_dates[gi] = d
        m = ustrip(u"mg", historical_egg_results[i].final_state.egg_mass)
        median_mass[gi] = earliest_mass[gi] = latest_mass[gi] = m
    elseif s.earliest !== nothing
        median_dates[gi] = Dates.value(hatch_date_of(s.median, issue_date))
        earliest_dates[gi] = Dates.value(hatch_date_of(s.earliest, issue_date))
        latest_dates[gi] = Dates.value(hatch_date_of(s.latest, issue_date))
        median_mass[gi] = mass_mg_of(s.median)
        earliest_mass[gi] = mass_mg_of(s.earliest)
        latest_mass[gi] = mass_mg_of(s.latest)
    end
    # historically-resolved points are a single outcome (n=1), not an ensemble
    total = s.status == :hatched_historically || startswith(string(s.status), "died_historically") ? 1 : length(members)
    frac_hatched[gi] = get(s.counts, :hatched, 0) / total
    frac_cold[gi] = get(s.counts, :cold, 0) / total
    frac_heat[gi] = get(s.counts, :heat, 0) / total
    frac_desiccation[gi] = get(s.counts, :desiccation, 0) / total
end

# Colourbars can't show date strings (Plots.jl GR limitation) -- build a
# series legend from dummy points instead, same trick as points_australia.jl.
function date_heatmap(title, ordinals)
    valid = filter(!isnan, ordinals)
    isempty(valid) && return add_basemap!(heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude"))
    lo, hi = extrema(valid)
    gradient = cgrad(:plasma)
    p = heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude", color=gradient, colorbar=false)
    tick_vals = round.(Int, range(lo, hi; length=min(6, length(unique(round.(Int, valid))))))
    for v in unique(tick_vals)
        scatter!(p, [NaN], [NaN]; color=gradient[(v - lo) / max(hi - lo, 1)], markersize=6,
            markerstrokewidth=0, label=string(Date(Dates.UTD(v))))
    end
    add_basemap!(p)
end

function mass_heatmap(title, vals)
    p = heatmap(lon_range, lat_range, to_heatmap_z(vals);
        title, xlabel="Longitude", ylabel="Latitude", color=cgrad(:viridis))
    add_basemap!(p)
end
function frac_heatmap(title, vals)
    p = heatmap(lon_range, lat_range, to_heatmap_z(vals);
        title, xlabel="Longitude", ylabel="Latitude", color=cgrad(:viridis), clims=(0, 1))
    add_basemap!(p)
end

plot_title = "Lay date $oviposition_date, issue date $issue_date, diapause=$diapause"

hatch_date_panel = plot(
    date_heatmap("Median hatch date", median_dates),
    date_heatmap("Earliest hatch date", earliest_dates),
    date_heatmap("Latest hatch date", latest_dates);
    layout=(2, 2), size=(1400, 1050), legendfontsize=6, plot_title,
)
hatch_date_path = joinpath(output_dir, "history_forecast_splice_hatch_dates.png")
savefig(hatch_date_panel, hatch_date_path)
println("Saved $hatch_date_path")
display(hatch_date_panel)

mass_panel = plot(
    mass_heatmap("Median egg mass at hatch (mg)", median_mass),
    mass_heatmap("Earliest egg mass at hatch (mg)", earliest_mass),
    mass_heatmap("Latest egg mass at hatch (mg)", latest_mass);
    layout=(2, 2), size=(1400, 1050), plot_title,
)
mass_path = joinpath(output_dir, "history_forecast_splice_mass.png")
savefig(mass_panel, mass_path)
println("Saved $mass_path")
display(mass_panel)

mortality_panel = plot(
    frac_heatmap("Fraction hatched", frac_hatched),
    frac_heatmap("Fraction died of heat", frac_heat),
    frac_heatmap("Fraction died of cold", frac_cold),
    frac_heatmap("Fraction died of desiccation", frac_desiccation);
    layout=(2, 2), size=(1400, 1050), plot_title,
)
mortality_path = joinpath(output_dir, "history_forecast_splice_mortality.png")
savefig(mortality_panel, mortality_path)
println("Saved $mortality_path")
display(mortality_panel)
