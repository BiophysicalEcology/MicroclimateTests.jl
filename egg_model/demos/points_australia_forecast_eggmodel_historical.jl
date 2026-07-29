# Stage 8, HPC split egg-model stage (historical only): runs the egg model
# over the historical (SILO) leg and caches the per-point results to disk --
# nothing else. Every forecast member's egg-model job
# (points_australia_forecast_eggmodel_members.jl) needs this cached first,
# both for forecast_point_indices (which points survive the historical leg)
# and each surviving point's final_state as its forecast init.
#
# Only depends on the historical microclimate leg (spartan_01), not the
# forecast members -- safe to run in parallel with
# points_australia_forecast_microclimate.jl's forecast-member array.

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", "c:/Spatial_Data/")

using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Serialization

include(joinpath(@__DIR__, "points_australia_forecast_setup.jl"))

save_trajectory = false

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
    HardTemperatureLimit(; lower_lethal_temperature, upper_lethal_temperature),
    StagedDesiccationLimit(;
        dry_mass, early_mass_factor,
        initial_egg_mass, late_mass_factor,
        ramp_start=quiescence_windows[1][1], ramp_end=quiescence_windows[1][2],
    ),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

# cache key includes diapause since it changes nest_node/arrest params, so
# toggling it can't silently reuse a cache written under the other setting.
historical_egg_cache_path = joinpath(egg_dir, "$(historical_label())_$(diapause ? "dia" : "nodia")_egg_n$(n).jls")

if isfile(historical_egg_cache_path) && use_cache
    println("Historical egg-model results already cached at $historical_egg_cache_path -- nothing to do.")
else
    println("Loading historical SILO microclimate for $n points...")
    historical_raw = solve_batched(build_historical_model(), historical_label(), points, historical_dates,
        (; soil_moisture=fill(0.2, length(depths))), nest_node)
    historical_day_range = 1:size(historical_raw.per_point[1].soil_temperature, 1)
    historical_tspan = (0.0u"hr", length(historical_day_range) * 1.0u"hr")
    historical_forcings = [egg_nest_forcing(historical_raw.per_point[i], historical_day_range, 1, environment_pars)
                            for i in 1:n]

    # ── egg-model cache pool: one integrator per thread, reinit! per work item ──
    nworkers = min(Threads.nthreads(), n)
    build_cache() = init_egg_cache(egg_model, pars, initial_state, soil_hydraulics,
        historical_forcings[1], historical_tspan; save_trajectory)
    cache_pool = Channel{typeof(build_cache())}(nworkers)
    for _ in 1:nworkers
        put!(cache_pool, build_cache())
    end

    function run_threaded!(results, indices, f)
        isempty(indices) && return
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
    # Point 1 is solved synchronously (negligible cost against the remaining
    # n-1, threaded) purely to learn simulate_egg!'s concrete return type --
    # preallocating Vector{Any} here and letting every thread write boxed
    # results into it is what drives GC/allocator lock contention under many
    # threads (this is what "N lock conflicts" in a solve's @time output is).
    c = take!(cache_pool)
    first_result = simulate_egg!(c, initial_state, soil_hydraulics, historical_forcings[1], historical_tspan)
    put!(cache_pool, c)

    historical_egg_results = Vector{typeof(first_result)}(undef, n)
    historical_egg_results[1] = first_result
    @time run_threaded!(historical_egg_results, 2:n, (cache, i) ->
        simulate_egg!(cache, initial_state, soil_hydraulics, historical_forcings[i], historical_tspan))

    serialize(historical_egg_cache_path, historical_egg_results)
    n_hatched_historically = count(r -> r.hatched, historical_egg_results)
    n_died_historically = count(r -> r.died, historical_egg_results)
    println("$n_hatched_historically/$n hatched and $n_died_historically/$n died before $issue_date; " *
            "$(n - n_hatched_historically - n_died_historically)/$n continue into the forecast.")
    println("Saved $historical_egg_cache_path")
end
