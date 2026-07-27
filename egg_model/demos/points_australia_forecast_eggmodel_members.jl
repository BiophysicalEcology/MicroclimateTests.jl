# Stage 8, HPC split egg-model stage (per-member forecast): loads the
# historical egg-model cache (written by
# points_australia_forecast_eggmodel_historical.jl) and one or more
# ACCESS-S2 forecast members' cached microclimate (written by
# points_australia_forecast_microclimate.jl), runs the egg model over each,
# and caches each member's outcomes to disk. No aggregation/maps here --
# that's points_australia_forecast_eggmodel_aggregate.jl, run afterward once
# every job below has completed.
#
# Usage: julia points_australia_forecast_eggmodel_members.jl <member_start> [member_end]
#
# Intended as a SLURM job array over 1:n_ensembles, depending on both the
# forecast-member microclimate array and the egg-model historical job.

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

member_start, member_end = if length(ARGS) >= 2
    parse(Int, ARGS[1]), parse(Int, ARGS[2])
elseif length(ARGS) == 1
    parse(Int, ARGS[1]), parse(Int, ARGS[1])
else
    error("usage: julia points_australia_forecast_eggmodel_members.jl <member_start> [member_end]")
end
1 <= member_start <= member_end <= n_ensembles || error("member range $member_start:$member_end out of bounds 1:$n_ensembles")

historical_egg_cache_path = joinpath(egg_dir, "$(historical_label())_$(diapause ? "dia" : "nodia")_egg_n$(n).jls")
isfile(historical_egg_cache_path) || error(
    "Historical egg-model cache not found at $historical_egg_cache_path -- " *
    "run points_australia_forecast_eggmodel_historical.jl first.")
println("Loading cached historical egg-model results...")
historical_egg_results = deserialize(historical_egg_cache_path)
forecast_point_indices = [i for i in 1:n
                           if !historical_egg_results[i].hatched && !historical_egg_results[i].died]

# Everything mutated while processing members (cache_pool, nworkers) lives
# inside this function, as genuine function-locals -- Julia's top-level
# "soft scope" ambiguity (reassigning a pre-existing name inside a `for`
# loop, which is otherwise ambiguous between rebinding the outer variable
# and shadowing it) only applies to top-level script scope, not function
# bodies, so no `global` annotations are needed anywhere here.
function run_forecast_members(member_start, member_end, forecast_point_indices, historical_egg_results)
    # Built lazily on the first member actually solved below (not every
    # member in this task's range is guaranteed to need solving -- some may
    # already be cached from a previous partial run) and reused across the
    # rest of this task's members.
    cache_pool = nothing
    nworkers = min(Threads.nthreads(), length(forecast_point_indices))

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

    for member in member_start:member_end
        outcome_path = joinpath(egg_dir, "$(forecast_label(member))_$(diapause ? "dia" : "nodia")_eggout_n$(n).jls")
        if isfile(outcome_path) && use_cache
            println("[member $member] Egg-model outcomes already cached, skipping.")
            continue
        end

        println("\nACCESS-S2 member $member/$(n_ensembles): loading microclimate for $n points...")
        forecast_raw = solve_batched(build_forecast_model(member), forecast_label(member), points, forecast_dates,
            (; soil_moisture=fill(0.2, length(depths))), nest_node)   # init unused on a cache hit
        forecast_day_range = 1:size(forecast_raw.per_point[1].soil_temperature, 1)
        forecast_forcings = Dict(i => egg_nest_forcing(forecast_raw.per_point[i], forecast_day_range, 1, environment_pars)
                                  for i in forecast_point_indices)

        first_i = forecast_point_indices[1]
        if cache_pool === nothing
            build_cache() = init_egg_cache(egg_model, pars, initial_state, soil_hydraulics,
                forecast_forcings[first_i], forecast_tspan; save_trajectory)
            cache_pool = Channel{typeof(build_cache())}(nworkers)
            for _ in 1:nworkers
                put!(cache_pool, build_cache())
            end
        end

        # Point `first_i` is solved synchronously purely to learn
        # simulate_egg!'s concrete return type -- preallocating Vector{Any}
        # and letting every thread write boxed results into it is what
        # drives GC/allocator lock contention under many threads.
        c = take!(cache_pool)
        first_result = simulate_egg!(c, historical_egg_results[first_i].final_state, soil_hydraulics,
            forecast_forcings[first_i], forecast_tspan)
        put!(cache_pool, c)

        member_outcomes = Vector{Union{Nothing,typeof(first_result)}}(nothing, n)
        member_outcomes[first_i] = first_result
        remaining_indices = filter(!=(first_i), forecast_point_indices)
        run_threaded!(member_outcomes, remaining_indices, (cache, i) ->
            simulate_egg!(cache, historical_egg_results[i].final_state, soil_hydraulics,
                forecast_forcings[i], forecast_tspan))

        serialize(outcome_path, member_outcomes)
        println("[member $member] Egg-model outcomes cached.")
    end
end

if isempty(forecast_point_indices)
    println("No points continue into the forecast (all hatched or died historically) -- nothing to do.")
else
    run_forecast_members(member_start, member_end, forecast_point_indices, historical_egg_results)
end
