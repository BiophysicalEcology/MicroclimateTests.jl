# "ct_egg_point.jl"
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
using Dates, Unitful, Plots, Statistics
using NCDatasets

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

# ── microclimate: SILO point run, extended with the soil layers the egg model needs ──

site = geocode("Mildura, Vic, Australia", buffer = 0.04)
site_name = split(site.display_name, ",")[1]
points = [site]

diapause = true

soil_type = :sandy_loam
bulk_density = 1.3u"Mg/m^3"
depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

dates = Date(2025, 1, 1):Day(1):Date(2027, 01, 31)
oviposition_dates = [Date(2026, 4, 25)]

# ACCESS-S2 data
issue_date = Date(2026, 7, 1)   # ACCESS-S2 issue date -- "now"
forecast_horizon_days = 214
forecast_dates = issue_date:Day(1):(issue_date + Day(forecast_horizon_days - 1))
historical_dates = dates[1]:Day(1):issue_date
members = collect(1:5)   # ACCESS-S2 ensemble members to simulate and plot

save_trajectory = true

if diapause
    nest_depth = 5.0u"cm"
    cold_hour_threshold = 30u"d"
    diapause_hour_threshold = 240.0u"d"
else
    nest_depth = 10.0u"cm"
    cold_hour_threshold = 0.0u"d"
    diapause_hour_threshold = 0.0u"d"
end
nest_node = nearest_node(nest_depth, depths)

output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
    LayerSpec(:reference_temperature, :scalar),
)

# Generic thread-pool runner: distributes `indices` across `nworkers` tasks,
# each borrowing a cache from `cache_pool`, applying `f(cache, i)` to every
# index it's given, and always returning its cache to the pool (even on
# error, so other workers don't deadlock waiting on take!).
function run_threaded!(f, results, indices, cache_pool, nworkers)
    work = Channel{eltype(indices)}(length(indices))
    foreach(i -> put!(work, i), indices)
    close(work)

    @sync for _ in 1:nworkers
        Threads.@spawn begin
            cache = take!(cache_pool)
            try
                for i in work
                    results[i] = f(cache, i)
                end
            finally
                put!(cache_pool, cache)
            end
        end
    end
    return results
end

# ── microclimate caching (NetCDF) ────────────────────────────────────────────
# Caches each leg's raw microclimate output (soil-profile time series, air
# temperature, and final soil state) so re-running with different egg-model
# parameters doesn't require re-solving the (expensive) microclimate model.
# Stored under a dedicated subdirectory with a `ctpoint_` prefix so these
# never collide with other demo scripts' outputs or the HPC batch pipeline's
# own cache files.

use_microclimate_cache = true

microclimate_cache_dir = joinpath(output_dir, "microclimate_cache")
mkpath(microclimate_cache_dir)

# Included in every cache label so changing the soil profile/grid invalidates
# old caches automatically, even though it's not otherwise part of the
# (site, dates, member) identity below. Egg-model parameters (pars, arrest,
# dm, survival_model, ...) deliberately do NOT appear here -- that's the
# whole point of this cache.
microclimate_config_tag() =
    string(hash((depths, heights, bulk_density, soil_type, nest_depth)); base=16)

safe_name(s) = replace(string(s), r"[^A-Za-z0-9_-]" => "_")

historical_cache_label() = "ctpoint_historical_$(safe_name(site_name))_" *
    "$(first(historical_dates))_$(last(historical_dates))_$(microclimate_config_tag())"

forecast_cache_label(m) = "ctpoint_forecast_$(safe_name(site_name))_" *
    "issue$(issue_date)_member$(m)_$(forecast_horizon_days)d_$(microclimate_config_tag())"

microclimate_cache_path(label) = joinpath(microclimate_cache_dir, label * ".nc")

# Fixed, known units per field -- more robust than round-tripping arbitrary
# printed unit strings through NetCDF attributes.
const MICROCLIMATE_UNITS = (
    soil_temperature          = u"K",
    soil_moisture             = NoUnits,
    soil_water_potential      = u"J/kg",
    soil_thermal_conductivity = u"W/m/K",
    soil_humidity             = NoUnits,
    reference_temperature     = u"K",
)

_to_plain(x, unit)   = unit === NoUnits ? collect(x) : ustrip.(unit, x)
_from_plain(x, unit) = unit === NoUnits ? x : x .* unit

datetime_to_hours(t) = Dates.value(t - DateTime(1970, 1, 1)) / 3_600_000
hours_to_datetime(h) = DateTime(1970, 1, 1) + Millisecond(round(Int, h * 3_600_000))

function save_microclimate_cache(path, result, air_times, air_values,
                                  final_soil_temperature, final_soil_moisture)
    isfile(path) && rm(path)
    ntime, ndepth = size(result.soil_temperature)

    NCDataset(path, "c") do ds
        defDim(ds, "time", ntime)
        defDim(ds, "depth", ndepth)
        defDim(ds, "air_time", length(air_times))

        for name in (:soil_temperature, :soil_moisture, :soil_water_potential,
                     :soil_thermal_conductivity, :soil_humidity)
            v = defVar(ds, string(name), Float64, ("time", "depth"))
            v[:, :] = _to_plain(getfield(result, name), getfield(MICROCLIMATE_UNITS, name))
        end

        va = defVar(ds, "reference_temperature", Float64, ("air_time",))
        va[:] = _to_plain(air_values, MICROCLIMATE_UNITS.reference_temperature)

        vt = defVar(ds, "reference_temperature_time", Float64, ("air_time",))
        vt[:] = datetime_to_hours.(DateTime.(air_times))

        vf1 = defVar(ds, "final_soil_temperature", Float64, ("depth",))
        vf1[:] = _to_plain(final_soil_temperature, MICROCLIMATE_UNITS.soil_temperature)

        vf2 = defVar(ds, "final_soil_moisture", Float64, ("depth",))
        vf2[:] = _to_plain(final_soil_moisture, MICROCLIMATE_UNITS.soil_moisture)
    end
    return path
end

function load_microclimate_cache(path)
    NCDataset(path, "r") do ds
        result = (;
            soil_temperature          = _from_plain(Array(ds["soil_temperature"]), MICROCLIMATE_UNITS.soil_temperature),
            soil_moisture             = _from_plain(Array(ds["soil_moisture"]), MICROCLIMATE_UNITS.soil_moisture),
            soil_water_potential      = _from_plain(Array(ds["soil_water_potential"]), MICROCLIMATE_UNITS.soil_water_potential),
            soil_thermal_conductivity = _from_plain(Array(ds["soil_thermal_conductivity"]), MICROCLIMATE_UNITS.soil_thermal_conductivity),
            soil_humidity             = _from_plain(Array(ds["soil_humidity"]), MICROCLIMATE_UNITS.soil_humidity),
        )

        air_values = _from_plain(Array(ds["reference_temperature"]), MICROCLIMATE_UNITS.reference_temperature)
        air_times  = hours_to_datetime.(Array(ds["reference_temperature_time"]))

        final_soil_temperature = _from_plain(Array(ds["final_soil_temperature"]), MICROCLIMATE_UNITS.soil_temperature)
        final_soil_moisture    = _from_plain(Array(ds["final_soil_moisture"]), MICROCLIMATE_UNITS.soil_moisture)

        (; result, air_times, air_values, final_soil_temperature, final_soil_moisture)
    end
end

environment_pars = example_environment_pars()

# Campbell & Norman (1998) Table 9.1 texture-class hydraulic parameters
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
    bulk_density, mineral_density=2.560u"Mg/m^3",
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
soil_profile = soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[soil_type], depths; bulk_density)
soil_hydraulics = (;
    air_entry_potential    = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b             = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)

# ── egg model ──

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
    HardTemperatureLimit(; lower_lethal_temperature, upper_lethal_temperature),
    DesiccationLimit(; dry_mass, critical_water_ratio),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

# ── historical leg: SILO, oviposition_date through issue_date ("now") ──
historical_model = MicroMapModel(;
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

historical_cache_file = microclimate_cache_path(historical_cache_label())

if use_microclimate_cache && isfile(historical_cache_file)
    println("Loading cached historical microclimate from $historical_cache_file")
    cached = load_microclimate_cache(historical_cache_file)
    historical_result       = cached.result
    historical_air_times    = cached.air_times
    historical_air_values   = cached.air_values
    now_soil_temperature    = cached.final_soil_temperature
    now_soil_moisture       = cached.final_soil_moisture
else
    println("Solving historical (SILO) microclimate: $oviposition_dates to $issue_date...")
    historical_problem = MicroVectorProblem(;
        model = historical_model, points, dates=historical_dates, soil_profile,
        init = (; soil_moisture = fill(0.2, length(depths))),
    )
    @time historical_output = solve(historical_problem)

    historical_day_range_tmp = 1:size(historical_output.soil_temperature[point=1], 1)
    historical_result = (;
        soil_temperature          = collect(historical_output.soil_temperature[point=1]),
        soil_moisture             = collect(historical_output.soil_moisture[point=1]),
        soil_water_potential      = collect(historical_output.soil_water_potential[point=1]),
        soil_thermal_conductivity = collect(historical_output.soil_thermal_conductivity[point=1]),
        soil_humidity             = collect(historical_output.soil_humidity[point=1]),
    )
    historical_air_times  = collect(lookup(historical_output.reference_temperature, Ti))
    historical_air_values = collect(historical_output.reference_temperature[point=1])
    now_soil_temperature  = collect(historical_output.soil_temperature[point=1, Ti=lastindex(historical_day_range_tmp)])
    now_soil_moisture     = collect(historical_output.soil_moisture[point=1, Ti=lastindex(historical_day_range_tmp)])

    save_microclimate_cache(historical_cache_file, historical_result, historical_air_times, historical_air_values,
        now_soil_temperature, now_soil_moisture)
    println("Cached historical microclimate to $historical_cache_file")
end

historical_day_range = 1:size(historical_result.soil_temperature, 1)
historical_forcing = egg_nest_forcing(historical_result, historical_day_range, nest_node, environment_pars)

# ── historical egg-model leg: oviposition_dates -> issue_date ──
max_duration = 720.0u"d"
day_range = 1:size(historical_result.soil_temperature, 1)
forcing_end_hr = length(day_range) * 1.0u"hr"

function run_from_oviposition(oviposition_date, forcing; save_trajectory=false)
    start_hr = oviposition_offset(oviposition_date, dates)

    if start_hr >= forcing_end_hr
        return (; hatched=false, hatch_time=missing, died=false, death_time=missing,
                  death_cause=:none, final_state=missing, trajectory=missing)
    end

    tspan = (start_hr, min(start_hr + max_duration, forcing_end_hr))
    initial_state = EggState(;
        egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
        maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
    )
    r = simulate_egg(egg_model, pars, initial_state, soil_hydraulics, forcing, tspan; save_trajectory)

    r isa NamedTuple ||
        error("simulate_egg returned unexpected result ($(typeof(r))) for oviposition_date=$oviposition_date")

    return r
end

println("Trying oviposition dates: ", oviposition_dates)
historical_egg_result = map(oviposition_dates) do oviposition_date
    r = run_from_oviposition(oviposition_date, historical_forcing; save_trajectory=true)
    start_hr = oviposition_offset(oviposition_date, dates)
    outcome = if r.hatched
        hatch_date = first(dates) + Day(round(Int, ustrip(u"d", r.hatch_time)))
        "hatched after $(round(typeof(1.0u"d"), r.hatch_time - start_hr, digits=1)) (on $hatch_date)"
    elseif r.died
        "died of $(r.death_cause) after $(round(typeof(1.0u"d"), r.death_time - start_hr, digits=1))"
    else
        "did not hatch in $(max_duration)"
    end
    println("  lay date $oviposition_date -> $outcome")
    r
end

# ── forecast leg: ACCESS-S2, issue_date + chosen ensemble members, for a given lay date ──
function run_forecast_ensemble(oviposition_date, historical_state, members)
    isempty(members) &&
        return (; forecast_outcomes=NamedTuple[], forecast_air=NamedTuple[], label="no members")

    forecast_tspan = (0.0u"hr", length(forecast_dates) * 24.0u"hr")

    # ── microclimate solves: sequential, one member at a time, cached to
    # NetCDF so re-running with different egg parameters skips these. ──
    forcings = Dict{eltype(members), Any}()
    airs     = Dict{eltype(members), NamedTuple}()

    for m in members
        cache_file = microclimate_cache_path(forecast_cache_label(m))

        if use_microclimate_cache && isfile(cache_file)
            println("[member $m] loading cached forecast microclimate")
            cached = load_microclimate_cache(cache_file)
            forecast_result = cached.result
            airs[m] = (; member=m, times=cached.air_times, temperature=cached.air_values)
        else
            forecast_model = MicroMapModel(;
                micro_model = MicroModel(;
                    depths, heights,
                    soil_properties_model = example_soil_properties_model(),
                    soil_hydraulic_model  = example_soil_hydraulic_model(),
                    snow_model            = NoSnow(),
                    config = MicroConfig(
                        soil_moisture_strategy = DynamicSoilMoisture(),
                        convergence = FixedSoilTemperatureIterations(1),
                    ),
                ),
                dem_source              = CRUCL2,
                weather_source          = AccessS2(issue_date, m),
                surface_albedo_source   = 0.15,
                roughness_height_source = 0.004u"m",
                compute_terrain         = false,
                output_layers,
            )

            forecast_problem = MicroVectorProblem(;
                model = forecast_model, points, dates=forecast_dates, soil_profile,
                init = (; soil_moisture = now_soil_moisture, soil_temperature = now_soil_temperature),
            )

            println("simulate microclimate (member $m)")
            @time forecast_output = solve(forecast_problem)

            air_times  = collect(lookup(forecast_output.reference_temperature, Ti))
            air_values = collect(forecast_output.reference_temperature[point=1])
            airs[m] = (; member=m, times=air_times, temperature=air_values)

            forecast_day_range_tmp = 1:size(forecast_output.soil_temperature[point=1], 1)
            forecast_result = (;
                soil_temperature          = collect(forecast_output.soil_temperature[point=1]),
                soil_moisture             = collect(forecast_output.soil_moisture[point=1]),
                soil_water_potential      = collect(forecast_output.soil_water_potential[point=1]),
                soil_thermal_conductivity = collect(forecast_output.soil_thermal_conductivity[point=1]),
                soil_humidity             = collect(forecast_output.soil_humidity[point=1]),
            )

            final_soil_temperature = collect(forecast_output.soil_temperature[point=1, Ti=lastindex(forecast_day_range_tmp)])
            final_soil_moisture    = collect(forecast_output.soil_moisture[point=1, Ti=lastindex(forecast_day_range_tmp)])

            save_microclimate_cache(cache_file, forecast_result, air_times, air_values,
                final_soil_temperature, final_soil_moisture)
            println("[member $m] cached forecast microclimate")
        end

        forecast_day_range = 1:size(forecast_result.soil_temperature, 1)
        forcings[m] = egg_nest_forcing(forecast_result, forecast_day_range, nest_node, environment_pars)
    end

    # ── egg model: cheap ODE integrations, threaded via a small pool of
    # reusable caches (one per worker), reinit'd per member inside simulate_egg!. ──
    nworkers = min(Threads.nthreads(), length(members))
    println("Running egg model for $(length(members)) members across $nworkers threads...")

    build_cache() = init_egg_cache(egg_model, pars, historical_state, soil_hydraulics,
        forcings[first(members)], forecast_tspan; save_trajectory=true)

    cache_pool = Channel{typeof(build_cache())}(nworkers)
    for _ in 1:nworkers
        put!(cache_pool, build_cache())
    end

    progress = Threads.Atomic{Int}(0)
    report_every = max(1, length(members) ÷ 20)
    outcomes = Dict{eltype(members), Any}()

    run_threaded!(outcomes, members, cache_pool, nworkers) do cache, m
        outcome = simulate_egg!(cache, historical_state, soil_hydraulics,
            forcings[m], forecast_tspan)

        if !hasproperty(outcome, :hatched) || !hasproperty(outcome, :died) || !hasproperty(outcome, :final_state)
            error("Invalid simulate_egg! result for lay date $oviposition_date, member $m: " *
                  "type=$(typeof(outcome)), value=$(repr(outcome))")
        end

        n = Threads.atomic_add!(progress, 1) + 1
        (n % report_every == 0 || n == length(members)) &&
            println("  $n/$(length(members)) ensemble members done")

        outcome
    end

    forecast_outcomes = [outcomes[m] for m in members]
    forecast_air = [airs[m] for m in members]

    println("\nPer-member outcomes (lay date $oviposition_date, issue date $issue_date):")
    hatch_dates = Date[]
    for (m, fr) in zip(members, forecast_outcomes)
        if fr.hatched
            hatch_date = issue_date + Day(round(Int, ustrip(u"d", fr.hatch_time)))
            push!(hatch_dates, hatch_date)
            println("  member $m -> hatched on $hatch_date")
        elseif fr.died
            println("  member $m -> died of $(fr.death_cause)")
        else
            println("  member $m -> did not hatch within $forecast_horizon_days days of the forecast")
        end
    end

    label = if isempty(hatch_dates)
        "0/$(length(members)) members hatched."
    else
        hatch_days = Dates.value.(hatch_dates .- issue_date)
        text = "$(length(hatch_dates))/$(length(members)) members hatched. " *
            "Median hatch: $(issue_date + Day(round(Int, median(hatch_days)))), " *
            "range: $(minimum(hatch_dates)) to $(maximum(hatch_dates))"
        println(text)
        text
    end

    return (; forecast_outcomes, forecast_air, label)
end

# ── drive the forecast leg for every oviposition date ──
results_by_lay_date = Dict{Date, NamedTuple}()

for (oviposition_date, r) in zip(oviposition_dates, historical_egg_result)
    if r.hatched
        println("Lay date $oviposition_date: hatched historically.")
    elseif r.died
        println("Lay date $oviposition_date: died historically ($(r.death_cause)).")
    elseif ismissing(r.final_state)
        println("Lay date $oviposition_date: no historical final state.")
    else
        println("Lay date $oviposition_date: still developing at issue date $issue_date " *
                "(development_fraction=$(r.final_state.development_fraction)); running forecast ensemble.")
        results_by_lay_date[oviposition_date] =
            run_forecast_ensemble(oviposition_date, r.final_state, members)
    end
end

# ── Plot helpers ─────────────────────────────────────────────────────────────

dtimes(base, t) = DateTime(base) .+ Millisecond.(round.(Int, ustrip.(u"ms", t)))
degC(x) = ustrip.(u"°C", collect(x))
mg(x)   = ustrip.(u"mg", collect(x))

function ensemble_series(x, nt)
    A = collect(x)
    ndims(A) == 1 && return [vec(A)]
    size(A, 1) == nt && return [vec(A[:, i]) for i in axes(A, 2)]
    size(A, 2) == nt && return [vec(A[i, :]) for i in axes(A, 1)]
    error("No dimension of $(size(A)) matches $nt forecast dates")
end

ensemble_median(x) = [median(getindex.(x, i)) for i in eachindex(first(x))]

# Historical and forecast series use different colours (not just linestyle)
# so the two legs remain distinguishable even where they overlap in time.
function add_comparison!(p, silo_t, silo_v, access_t, access_v;
                          historical_color, forecast_color,
                          historical_label="SILO", forecast_label="ACCESS-S2 median")
    plot!(p, silo_t, silo_v; color=historical_color, linewidth=2, label=historical_label)
    for v in access_v
        plot!(p, access_t, v; color=forecast_color, alpha=0.15, linewidth=0.8, label="")
    end
    plot!(p, access_t, ensemble_median(access_v);
          color=forecast_color, linewidth=2.5, label=forecast_label)
end

# ── Weather data ─────────────────────────────────────────────────────────────
# SILO is loaded through the forecast period so it overlaps ACCESS-S2.

dates_historical_weather =
    Date(first(historical_dates)):Day(1):(Date(issue_date) + Day(forecast_horizon_days - 1))

weather_historical = getpoint(
    SILO,
    (:daily_rain, :max_temp, :min_temp, :radiation, :vp);
    lon=site.lon, lat=site.lat,
    date=dates_historical_weather,
    username="m.kearney@unimelb.edu.au",
)

nh = length(weather_historical.daily_rain.values)
silo_dates = collect(dates_historical_weather)[1:nh]
silo_times = DateTime.(silo_dates) .+ Hour(12)

silo = (
    rain = weather_historical.daily_rain.values,
    tmin = weather_historical.min_temp.values,
    tmax = weather_historical.max_temp.values,
    radn = weather_historical.radiation.values,
    vapr = weather_historical.vp.values,
)

years = MicroclimateMapper._years_from_dates(forecast_dates)

start_ti = Dates.value(Date(first(forecast_dates)) - Date(first(years), 1, 1)) + 1
ti = start_ti:(start_ti + length(forecast_dates) - 1)
access_dates = Date.(collect(forecast_dates))
access_times = DateTime.(access_dates) .+ Hour(12)

# AccessS2 takes a single member, so load each requested member separately
# and collect the results into one vector-of-series per variable.
access_by_member = map(members) do m
    weather_forecast = MicroclimateMapper._load_layers(
        AccessS2Loader(),
        AccessS2(issue_date, m),
        (:tmin, :tmax, :rain, :radn, :vapr),
        site.extent,
        years,
    )
    (;
        rain = vec(collect(weather_forecast[:rain][X=1, Y=1, Ti(ti)])),
        tmin = vec(collect(weather_forecast[:tmin][X=1, Y=1, Ti(ti)])),
        tmax = vec(collect(weather_forecast[:tmax][X=1, Y=1, Ti(ti)])),
        radn = vec(collect(weather_forecast[:radn][X=1, Y=1, Ti(ti)])),
        vapr = vec(collect(weather_forecast[:vapr][X=1, Y=1, Ti(ti)])),
    )
end

access = (
    rain = [s.rain for s in access_by_member],
    tmin = [s.tmin for s in access_by_member],
    tmax = [s.tmax for s in access_by_member],
    radn = [s.radn for s in access_by_member],
    vapr = [s.vapr for s in access_by_member],
)

# ── Figure 1: SILO versus ACCESS-S2 ─────────────────────────────────────────

w1 = plot(ylabel="temperature, °C",
          title="SILO and ACCESS-S2 weather comparison at $site_name")
w2 = plot(ylabel="rainfall, mm/day")
w3 = plot(ylabel="radiation, MJ/m²/day")
w4 = plot(ylabel="vapour pressure", xlabel="date")

comparisons = (
    (panel=w1, key=:tmin, hist=:navy,       fcst=:dodgerblue, hlab="SILO minimum", flab="ACCESS-S2 median minimum"),
    (panel=w1, key=:tmax, hist=:firebrick,  fcst=:orange,     hlab="SILO maximum", flab="ACCESS-S2 median maximum"),
    (panel=w2, key=:rain, hist=:black,      fcst=:royalblue,  hlab="SILO",         flab="ACCESS-S2 median"),
    (panel=w3, key=:radn, hist=:saddlebrown,fcst=:darkorange, hlab="SILO",         flab="ACCESS-S2 median"),
    (panel=w4, key=:vapr, hist=:indigo,     fcst=:orchid,     hlab="SILO",         flab="ACCESS-S2 median"),
)

for c in comparisons
    add_comparison!(c.panel, silo_times, getfield(silo, c.key), access_times, getfield(access, c.key);
                     historical_color=c.hist, forecast_color=c.fcst,
                     historical_label=c.hlab, forecast_label=c.flab)
end

overlap_start, overlap_end = DateTime(first(access_dates)), DateTime(last(access_dates)) + Day(1)
comparison_start = max(DateTime(first(silo_dates)), overlap_start - Day(30))

for p in (w1, w2, w3, w4)
    vspan!(p, [overlap_start, overlap_end]; color=:grey, alpha=0.15, label="")
    vline!(p, [DateTime(issue_date)]; color=:black, linestyle=:dot, label="")
end

weather_plot = plot(
    w1, w2, w3, w4;
    layout=(4, 1), link=:x, size=(1100, 900),
    xlims=(comparison_start, overlap_end),
    left_margin=8Plots.mm, bottom_margin=4Plots.mm,
)

weather_path = joinpath(@__DIR__, "silo_access_weather_comparison.png")
savefig(weather_plot, weather_path)
display(weather_plot)
println("Saved weather comparison to $weather_path")

# ── Figure 2: egg trajectories and forcing ──────────────────────────────────

p1 = plot(ylabel="development",
          title="Egg development by lay date at $site_name, diapause = $diapause",
          legend=false, ylims=(0, 1))
p2 = plot(ylabel="egg mass, mg", legend=false,
          ylims=(0, ustrip(u"mg", pars.initial_egg_mass * 2.5)))
p3 = plot(ylabel="temperature, °C", legend=false)
p4 = plot(ylabel="rainfall,\nmm/day", xlabel="date", legend=false)

colors = palette(:default, length(oviposition_dates))

historical_air_temperature = degC(historical_air_values)
plot!(p3, historical_air_times, historical_air_temperature;
      color=:grey35, linewidth=1.5, label="historical air")

# Historical egg trajectories.
for (i, (lay_date, result)) in enumerate(zip(oviposition_dates, historical_egg_result))
    hasproperty(result, :trajectory) || continue
    ismissing(result.trajectory) && continue

    traj = result.trajectory
    t = dtimes(first(dates), traj.t)
    c = colors[i]

    plot!(p1, t, traj.development_fraction; color=c, linewidth=1, label=string(lay_date))
    plot!(p2, t, mg(traj.egg_mass); color=c, linewidth=1, label=string(lay_date))
    plot!(p3, t, degC(traj.temperature); color=c, linewidth=1, label=string(lay_date))
end

# Forecast egg trajectories (all chosen ensemble members, all lay dates).
for (i, lay_date) in enumerate(oviposition_dates)
    haskey(results_by_lay_date, lay_date) || continue
    c = colors[i]

    for result in results_by_lay_date[lay_date].forecast_outcomes
        hasproperty(result, :trajectory) || continue
        ismissing(result.trajectory) && continue

        traj = result.trajectory
        t = dtimes(issue_date, traj.t)

        plot!(p1, t, traj.development_fraction; color=c, linewidth=1, label="")
        plot!(p2, t, mg(traj.egg_mass); color=c, linewidth=1, label="")
        plot!(p3, t, degC(traj.temperature); color=c, linewidth=1, label="")
    end
end

# Forecast microclimate air temperature is the same across lay dates, so
# take it from the first lay date that required a forecast.
if !isempty(results_by_lay_date)
    first_forecast_result = first(values(results_by_lay_date))
    for air in first_forecast_result.forecast_air
        plot!(p3, air.times, degC(air.temperature);
              color=:steelblue, alpha=0.4, linewidth=0.8, label="")
    end
end

# Rainfall: SILO before issue date, ACCESS-S2 median from issue date onward.
historical_rain_indices = findall(silo_dates .< Date(issue_date))
bar!(p4, silo_times[historical_rain_indices], silo.rain[historical_rain_indices];
     color=:black, linewidth=0, label="SILO")
bar!(p4, access_times, ensemble_median(access.rain);
     color=:royalblue, linewidth=0, label="ACCESS-S2 median")

for p in (p1, p2, p3, p4)
    vline!(p, [DateTime(issue_date)]; color=:black, linestyle=:dash, linewidth=1.25, label="")
end
hline!(p3, [ustrip(u"°C", upper_lethal_temperature)]; color=:red, linestyle=:dash, linewidth=1.25, label="")
hline!(p3, [ustrip(u"°C", lower_lethal_temperature)]; color=:blue, linestyle=:dash, linewidth=1.25, label="")

plot_start = DateTime(minimum(oviposition_dates))
plot_end = DateTime(last(forecast_dates)) + Day(1)

egg_plot = plot(
    p1, p2, p3, p4;
    layout=(4, 1), size=(1100, 900), link=:x,
    xlims=(plot_start, plot_end),
    left_margin=8Plots.mm, bottom_margin=4Plots.mm,
)

egg_path = joinpath(@__DIR__, "egg_historical_forecast_trajectories.png")
savefig(egg_plot, egg_path)
display(egg_plot)
println("Saved egg plot to $egg_path")

# ── Figure 3: histogram of forecast hatch dates, one panel per lay date ────

lay_dates_with_forecast = [d for d in oviposition_dates if haskey(results_by_lay_date, d)]

if isempty(lay_dates_with_forecast)
    println("No lay dates had a forecast ensemble; skipping hatch-date histogram.")
else
    hatch_panels = map(lay_dates_with_forecast) do lay_date
        forecast = results_by_lay_date[lay_date]
        c = colors[findfirst(==(lay_date), oviposition_dates)]
        n_members = length(forecast.forecast_outcomes)

        hatch_days = [
            Dates.value(Date(issue_date) + Day(round(Int, ustrip(u"d", fr.hatch_time))) - issue_date)
            for fr in forecast.forecast_outcomes if fr.hatched
        ]

        p = plot(; title="lay date $lay_date ($(length(hatch_days))/$n_members hatched)",
                  xlabel="days after issue date ($issue_date)", ylabel="members",
                  legend=false, titlefontsize=9)

        if isempty(hatch_days)
            annotate!(p, 0.5, 0.5, text("no members hatched", 9))
        else
            histogram!(p, hatch_days; bins=:auto, color=c, linewidth=0)
            vline!(p, [median(hatch_days)]; color=:black, linestyle=:dash, linewidth=1.5,
                   label="median")
        end
        p
    end

    hatch_plot = plot(hatch_panels...;
        layout=(length(hatch_panels), 1),
        size=(800, 300 * length(hatch_panels)),
        left_margin=6Plots.mm, bottom_margin=4Plots.mm,
    )

    hatch_path = joinpath(@__DIR__, "hatch_date_histograms.png")
    savefig(hatch_plot, hatch_path)
    display(hatch_plot)
    println("Saved hatch-date histograms to $hatch_path")
end