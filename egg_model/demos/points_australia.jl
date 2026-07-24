# Stage 6 Pass A: broad first-look egg-development sweep over a regular grid
# of points across south-eastern Australia, via MicroVectorProblem (SILO),
# assuming one uniform/optimal Campbell & Norman soil texture rather than
# fetching real per-point texture (that's Pass B, raster-mode, once this
# broad pass identifies a region worth zooming into).
#
# Reuses point_silo_deterministic.jl's exact model/egg-model setup, just with
# many points instead of one, and reuses egg_model/src/phases.jl's
# init_egg_cache/simulate_egg! cache-reuse API (one integrator cache per
# thread, reinit! per point) instead of building a fresh integrator per point
# -- the same Channel-based worker-pool pattern MicroclimateMapper.jl's own
# grid solve uses (src/common.jl).

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_properties_model, example_soil_hydraulic_model
using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Rasters, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using DataInterpolations
using Dates, Unitful
using Serialization
using NaturalEarth, GeoInterface

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))

# all cached/serialized run output goes here
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"#"z:/"

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
nest_depth = 10.0u"cm"

soil_source = :sandy_loam

dates = Date(2010, 1, 1):Day(1):Date(2011, 12, 31) # for microclimate
max_duration = 720.0u"d"
forcing_end_hr = length(day_range) * 1.0u"hr"
oviposition_date = Date(2010, 4, 1)

# layers wanted from the microclimate model runs
# egg_nest_forcing (forcing.jl) needs only these
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
        config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
    ),
    dem_source              = CRUCL2,
    weather_source          = SILO,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
    output_layers,
)

# regular grid over south-eastern Australia
lon_range = range(135.0, 153.0; length=50)
lat_range = range(-39.0, -24.0; length=40)
all_grid_points = vec([(lon, lat) for lon in lon_range, lat in lat_range])   # 2000 points, column-major (lon fastest)

# Reject points that fall in the ocean based on the DEM and the silo grid
# TODO work out what DEM SILO uses anduse that instead, probably only need than and not the CRUCL2_ELV filter
const CRUCL2_ELV = read(Raster(RasterDataSources.getraster(CRUCL2); name=:elv, lazy=true))
has_crucl2_land(lon, lat) = !ismissing(CRUCL2_ELV[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_crucl2_land(p...), all_grid_points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the CRUCL2 land-mask pre-check.")

# using reference year (2020)
const SILO_MAXTEMP_PROBE = read(Raster(RasterDataSources.getraster(SILO, :max_temp; date=Date(2020,1,1)); name=:max_temp, lazy=true)[Ti(1)])
has_silo_land(lon, lat) = !ismissing(SILO_MAXTEMP_PROBE[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_silo_land(p...), points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the SILO land-mask pre-check.")

# plot planned points to simulate
using Plots
let
    domain = Extent(X = (minimum(lon_range) - 1, maximum(lon_range) + 1),
                     Y = (minimum(lat_range) - 1, maximum(lat_range) + 1))
    elv_domain = crop(CRUCL2_ELV; to = domain, touches = true)
    p = plot(elv_domain; title = "CRUCL2 elevation + grid points (green=kept, red=rejected)")
    rejected = setdiff(all_grid_points, points)
    scatter!(p, first.(points), last.(points); markersize = 2, markerstrokewidth = 0, color = :green, label = "kept")
    scatter!(p, first.(rejected), last.(rejected); markersize = 2, markerstrokewidth = 0, color = :red, label = "rejected")
    savefig(p, joinpath(@__DIR__, "points_australia_domain_check.png"))
    display(p)
end

n_points_to_run = length(points)

# whole calendar year(s) only, per current SILO/point-query constraints.
use_cache = true

# ── uniform soil texture options ──
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
            campbell_b_parameter=fill(texture.b, n),
            root_density,
        ),
    )
end

soil_profile = soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[soil_source], depths)

nest_node = nearest_node(nest_depth, depths)
environment_pars = example_environment_pars()

# Batched solve MicroVectorProblem over all n_points_to_run. Each batch is
# also cached to disk.
batch_size = 100

n = length(points)
n_batches = cld(n, batch_size)
forcings = Vector{Any}(undef, n)
day_range = 1:(length(dates) * 24)

forcings_cache_path = joinpath(output_dir, "points_australia_forcings_n$(n_points_to_run).jls")
if isfile(forcings_cache_path) && use_cache
    println("Loading cached forcings from $forcings_cache_path...")
    forcings = deserialize(forcings_cache_path)
else
    for b in 1:n_batches
        i_start, i_end = (b - 1) * batch_size + 1, min(b * batch_size, n)
        batch_points = points[i_start:i_end]
        batch_cache_path = joinpath(output_dir, "points_australia_batch$(b)of$(n_batches)_n$(n_points_to_run).jls")

        batch_forcings = if isfile(batch_cache_path) && use_cache
            println("Loading cached batch $b/$n_batches ($(length(batch_points)) points)...")
            deserialize(batch_cache_path)
        else
            println("Solving SILO microclimate for batch $b/$n_batches ($(length(batch_points)) points)...")
            batch_problem = MicroVectorProblem(;
                model, points=batch_points, dates, soil_profile,
                init = (; soil_moisture = fill(0.2, length(depths))),
            )
            @time batch_output = solve(batch_problem)
            result = map(1:length(batch_points)) do i
                result_i = (;
                    soil_temperature          = collect(batch_output.soil_temperature[point=i]),
                    soil_moisture             = collect(batch_output.soil_moisture[point=i]),
                    soil_water_potential      = collect(batch_output.soil_water_potential[point=i]),
                    soil_thermal_conductivity = collect(batch_output.soil_thermal_conductivity[point=i]),
                    soil_humidity             = collect(batch_output.soil_humidity[point=i]),
                )
                egg_nest_forcing(result_i, day_range, nest_node, environment_pars)
            end
            serialize(batch_cache_path, result)
            batch_output = nothing
            batch_problem = nothing
            GC.gc()
            result
        end
        forcings[i_start:i_end] .= batch_forcings
    end
    serialize(forcings_cache_path, forcings)
end

soil_hydraulics = (;
    air_entry_potential    = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b             = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)

# ── egg model, identical config to point_silo_deterministic.jl ──

shape_b = 0.69 / 1.82
geometry = Ellipsoid(0.0036u"g", 1000.0u"kg/m^3", 1 / shape_b, 1 / shape_b)
arrest = ProportionWindowArrest(;
    cold_temperature=u"K"(0.0u"°C"), diapause_window=(0.45, 0.50),
    quiescence_windows=((0.25, 0.30), (0.45, 0.50)),
    cold_hour_threshold=0.0u"hr", diapause_hour_threshold=0.0u"hr",
    desiccation_tolerance=0.6,
)
dm = arrhenius_development_model(;
    T_A=6641.6175, T_AL=33600.0, T_AH=48000.0, T_L=289.15, T_H=314.65, T_ref=301.65,
    rate_at_reference=1 / 17.4, rate_unit=1.0u"d^-1",
)
base_K_e = 2.347802e-9 * u"kg/m^2/s/(J/kg)"
stage = SteppedHydricStage(;
    conductance_threshold=0.25, wetness_threshold=0.45,
    dormant_conductance=0.0u"kg/m^2/s/(J/kg)", active_conductance=base_K_e * 3,
    dormant_wetness=0.35 / 100, active_wetness=0.35,
)
pars = EggParameters(;
    hydraulic_conductance=base_K_e, specific_hydration=0.000304u"m^3/m^3/(J/kg)",
    conduction_fraction=0.5, skin_wetness=0.35 / 100,
    initial_egg_mass=0.0036u"g", minimum_egg_mass=0.0026u"g",
)
survival_model = CombinedSurvival(
    HardTemperatureLimit(; lower_lethal_temperature=u"K"(-5.0u"°C"), upper_lethal_temperature=u"K"(52.0u"°C")),
    DesiccationLimit(; dry_mass=0.1 * pars.initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

start_hr = oviposition_offset(oviposition_date, dates)
tspan = (start_hr, min(start_hr + max_duration, forcing_end_hr))
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

# ── grid loop: one egg-model integrator cache per thread, reused across all
# points via reinit! (init_egg_cache/simulate_egg!, phases.jl) -- the same
# per-thread-not-per-cell pattern MicroclimateMapper.jl's own grid solve uses
# (src/common.jl's Channel-based worker pool).

println("Running egg model at $(length(points)) points, lay date $oviposition_date...")
n = length(points)
# in test_mode, cap workers well below n so cache reuse across *different*
# locations is actually exercised (with nworkers>=n every point would get its
# own dedicated cache, never reinit!'d for a second location) -- the real
# correctness property Stage 6 depends on.
nworkers = test_mode ? min(3, n) : min(Threads.nthreads(), n)

build_cache() = init_egg_cache(egg_model, pars, initial_state, soil_hydraulics, forcings[1], tspan)
cache_pool = Channel{typeof(build_cache())}(nworkers)
put!(cache_pool, build_cache())
for _ in 2:nworkers
    put!(cache_pool, build_cache())
end

work = Channel{Int}(n)
for i in 1:n
    put!(work, i)
end
close(work)

results = Vector{Any}(undef, n)
# report ~20 times over the whole run regardless of n -- Threads.Atomic since
# multiple worker tasks increment this concurrently.
progress = Threads.Atomic{Int}(0)
report_every = max(1, n ÷ 20)
@time @sync for _ in 1:nworkers
    Threads.@spawn begin
        cache = take!(cache_pool)
        for i in work
            results[i] = simulate_egg!(cache, initial_state, soil_hydraulics, forcings[i], tspan)
            done = Threads.atomic_add!(progress, 1) + 1
            (done % report_every == 0 || done == n) && println("  $done/$n egg simulations done")
        end
        put!(cache_pool, cache)
    end
end

# for (i, r) in enumerate(results)
#     lon, lat = points[i]
#     outcome = if r.hatched
#         hatch_date = first(dates) + Day(round(Int, ustrip(u"d", r.hatch_time)))
#         "hatched on $hatch_date"
#     elseif r.died
#         "died of $(r.death_cause)"
#     else
#         "did not hatch in $(max_duration)"
#     end
#     println("  ($lon, $lat) -> $outcome")
# end

# ── rasterize for plotting: `all_grid_points` is the full regular 50x40 grid
# (lon fastest), but `points`/`points` is a filtered subset (ocean points
# removed) -- not every grid cell has a result, so this can't be a plain
# reshape (that failed: "new dimensions (50, 40) must be consistent with
# array length 1640"). Place each result back at its original grid index,
# leaving filtered-out cells as NaN, then reshape the *full-length* array.
if n_points_to_run == length(points)
    using Plots
    point_to_index = Dict(p => i for (i, p) in enumerate(all_grid_points))

    # heatmap(x, y, z) wants z sized (length(y), length(x)) -- transpose of
    # the (lon, lat)-shaped grid built everywhere else in this file.
    to_heatmap_z(vals) = permutedims(reshape(vals, length(lon_range), length(lat_range)))

    # ── background: state/territory boundaries + a few reference towns,
    # matching the general look of the APLC forecaster's own hatching-
    # prediction maps (state lines, coastline, labelled towns). Natural
    # Earth doesn't carry the APLC's own named sub-regions (Riverina,
    # Mallee, etc. -- those are BOM/APLC-internal management zones, not a
    # standard public boundary dataset), so this is state-level only.
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

    # overlays state boundaries + labelled town markers onto any heatmap `p`.
    # Adding the (whole-country-spanning) boundary lines expands Plots.jl's
    # auto axis limits to fit all of Australia -- reset to the actual data
    # domain afterward so the heatmap itself stays the visual focus.
    function add_basemap!(p)
        plot!(p, state_xs, state_ys; color=:black, linewidth=0.75, label=nothing)
        scatter!(p, [t.lon for t in towns], [t.lat for t in towns];
            color=:black, markersize=2, label=nothing)
        for t in towns
            town_name = first(split(t.display_name, ","))
            annotate!(p, t.lon, t.lat, text(town_name, 6, :left, :bottom))
        end
        xlims!(p, extrema(lon_range)...)
        ylims!(p, extrema(lat_range)...)
        p
    end

    hatch_days = fill(NaN, length(all_grid_points))
    for (i, r) in enumerate(results)
        r.hatched && (hatch_days[point_to_index[points[i]]] = ustrip(u"d", r.hatch_time))
    end
    hatch_plot = heatmap(lon_range, lat_range, to_heatmap_z(hatch_days);
        title="Hatch time (days since $(first(dates)))",
        xlabel="Longitude", ylabel="Latitude",
    )
    add_basemap!(hatch_plot)

    # Same data, plotted as actual calendar hatch dates instead of an elapsed
    # duration. Plots.jl's GR backend colourbar doesn't support custom string
    # tick labels (confirmed: `colorbar_ticks=(vals, labels)` silently falls
    # back to plain numbers) -- so instead the colourbar is hidden and a
    # normal *series* legend is built from a handful of invisible dummy
    # points, one per representative date, each coloured by sampling the
    # same gradient the heatmap uses. Series legends render arbitrary
    # strings fine; colourbars don't.
    hatch_date_ordinals = fill(NaN, length(all_grid_points))
    for (i, r) in enumerate(results)
        if r.hatched
            hatch_date = first(dates) + Day(round(Int, ustrip(u"d", r.hatch_time)))
            hatch_date_ordinals[point_to_index[points[i]]] = Dates.value(hatch_date)
        end
    end
    valid_ordinals = filter(!isnan, hatch_date_ordinals)
    lo, hi = extrema(valid_ordinals)
    tick_vals = round.(Int, range(lo, hi; length=6))
    tick_labels = string.(Date.(Dates.UTD.(tick_vals)))
    date_gradient = cgrad(:plasma)
    hatch_date_plot = heatmap(lon_range, lat_range, to_heatmap_z(hatch_date_ordinals);
        title="Hatch date",
        xlabel="Longitude", ylabel="Latitude",
        color=date_gradient, colorbar=false,
    )
    for (v, label) in zip(tick_vals, tick_labels)
        scatter!(hatch_date_plot, [NaN], [NaN];
            color=date_gradient[(v - lo) / (hi - lo)], markersize=6, markerstrokewidth=0, label=label)
    end
    add_basemap!(hatch_date_plot)

    # Outcome map: hatched / died (by cause) / survived-but-didn't-hatch in
    # max_duration, as a small integer code per cell (categorical, not a
    # continuous scale) -- same NaN-for-filtered-out convention as above.
    # death_cause values come from egg_model/src/types.jl's cause_of_death;
    # extend both dicts below if egg_model's survival_model ever gains
    # another cause. Colours are the Okabe-Ito colourblind-safe palette.
    # Same dummy-series-legend trick as the date plot above, for the same
    # reason (colourbar can't show the category name strings).
    outcome_codes = Dict(:hatched => 0, :cold => 1, :heat => 2, :desiccation => 3, :timeout => 4)
    outcome_labels = ["hatched", "died: cold", "died: heat", "died: desiccation", "no hatch (timeout)"]
    outcome_colors = ["#009E73", "#0072B2", "#D55E00", "#E69F00", "#999999"]
    outcomes = fill(NaN, length(all_grid_points))
    for (i, r) in enumerate(results)
        code = r.hatched ? :hatched : r.died ? r.death_cause : :timeout
        outcomes[point_to_index[points[i]]] = outcome_codes[code]
    end
    present_codes = Int.(sort(unique(filter(!isnan, outcomes))))
    outcome_plot = heatmap(lon_range, lat_range, to_heatmap_z(outcomes);
        title="Outcome",
        xlabel="Longitude", ylabel="Latitude",
        color=cgrad(outcome_colors[present_codes.+1]; categorical=true),
        clims=(minimum(present_codes) - 0.5, maximum(present_codes) + 0.5),
        colorbar=false,
    )
    for code in present_codes
        scatter!(outcome_plot, [NaN], [NaN];
            color=outcome_colors[code+1], markersize=6, markerstrokewidth=0, label=outcome_labels[code+1])
    end
    add_basemap!(outcome_plot)

    # Egg mass at hatch (mg) -- final_state.egg_mass at the moment of
    # hatching; NaN for cells that didn't hatch, same convention as above.
    hatch_mass_mg = fill(NaN, length(all_grid_points))
    for (i, r) in enumerate(results)
        r.hatched && (hatch_mass_mg[point_to_index[points[i]]] = ustrip(u"mg", r.final_state.egg_mass))
    end
    mass_plot = heatmap(lon_range, lat_range, to_heatmap_z(hatch_mass_mg);
        title="Egg mass at hatch (mg)",
        xlabel="Longitude", ylabel="Latitude",
        color=cgrad(:viridis),
    )
    add_basemap!(mass_plot)

    panel = plot(hatch_plot, hatch_date_plot, outcome_plot, mass_plot;
        layout=(2, 2), size=(1400, 1050),
        plot_title="Lay date $oviposition_date", legendfontsize=6,
    )
    savefig(panel, joinpath(@__DIR__, "points_australia_panel.png"))
    display(panel)
else
    println("\nSkipping the hatch-date raster plot -- only $(n_points_to_run)/$(length(points)) points were run " *
            "(set test_mode = false for the full Pass A sweep).")
end
