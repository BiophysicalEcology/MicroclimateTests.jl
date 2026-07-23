using OrdinaryDiffEqTsit5
using StaticArrays
using Unitful

# Coupled ODE state, time in hours, plain Float64 (units stripped for the
# solver): [development_fraction, egg_mass (kg), egg_water_potential (J/kg),
# chill_accumulation (hr), diapause_duration (hr), maximum_mass_achieved (kg)].
# maximum_mass_achieved ratchets up with mass, never down.
# RHS is discontinuous at every threshold crossing, so each chunk gets a fresh
# ODEProblem (bout-chunking) rather than one continuous integration.

function _unpack(u)
    development_fraction = u[1]
    egg_mass = u[2] * u"kg"
    egg_water_potential = u[3] * u"J/kg"
    arrest_state = (; chill_accumulation=u[4] * u"hr", diapause_duration=u[5] * u"hr")
    maximum_mass_achieved = u[6] * u"kg"
    (; development_fraction, egg_mass, egg_water_potential, arrest_state, maximum_mass_achieved)
end

function egg_rhs(u, p, t)
    (; egg_model, pars, soil_hydraulics, forcing) = p
    (; development_fraction, egg_mass, egg_water_potential, arrest_state, maximum_mass_achieved) = _unpack(u)
    (; environment, soil_water_potential) = forcing(t * u"hr")

    state = EggState(; development_fraction, egg_mass, egg_water_potential, maximum_mass_achieved, arrest_state)
    core_temperature = egg_temperature(egg_model.thermal_model, egg_model, state, pars, environment)
    hyd_index = hydration_index(egg_model.hydric_model, state, pars, soil_water_potential)
    halted = arrested(egg_model.arrest_model, arrest_state, development_fraction, hyd_index)

    d_dev = halted ? 0.0 : ustrip(u"hr^-1", development_rate(egg_model.development_model, core_temperature))
    (; d_mass, d_water_potential) = hydric_rate(
        egg_model.hydric_model, egg_model, state, pars, environment, soil_water_potential, soil_hydraulics,
    )
    is_cold = core_temperature < egg_model.arrest_model.cold_temperature
    d_chill = is_cold ? 1.0 : 0.0
    d_diapause_duration = in_diapause(egg_model.arrest_model, development_fraction, arrest_state) ? 1.0 : 0.0
    # ratchet: only rises, and only while mass is at the current max
    d_max_mass = egg_mass >= maximum_mass_achieved ? max(0.0, ustrip(u"kg/hr", d_mass)) : 0.0

    SVector(d_dev, ustrip(u"kg/hr", d_mass), ustrip(u"J/kg/hr", d_water_potential), d_chill, d_diapause_duration, d_max_mass)
end

# hydration_index as a function of (u[2], u[6]) directly -- desiccation_tolerance
# crossing is a condition on this, not a fixed mass threshold (m_max is dynamic).
_hydration_index_u(u, pars) =
    clamp((u[2] - ustrip(u"kg", pars.initial_egg_mass)) / (u[6] - ustrip(u"kg", pars.initial_egg_mass)), 0.0, 1.0)

# flattens a Tuple-of-Tuples via recursive tuple destructuring -- unlike
# reduce(vcat, ...)/Vector-building, this is fully unrolled/type-stable at
# compile time for any fixed-length Tuple, which is what keeps _conditions'
# return type (and hence every CallbackSet built from it) concrete regardless
# of how many quiescence windows or survival sub-models are configured.
_flatten_tuples(::Tuple{}) = ()
_flatten_tuples(t::Tuple) = (first(t)..., _flatten_tuples(Base.tail(t))...)

# DesiccationLimit's death-mass threshold(s), recursive for CombinedSurvival.
# Root-found (not just checked post-hoc), since the hydric ODE goes singular as
# egg_mass -> 0 and an explicit solver can outrun a post-hoc check within a step.
_survival_mass_thresholds(::NoSurvivalLimit) = ()
_survival_mass_thresholds(::HardTemperatureLimit) = ()
_survival_mass_thresholds(m::DesiccationLimit) = (ustrip(u"kg", m.dry_mass * (1 + m.critical_water_ratio)),)
_survival_mass_thresholds(m::CombinedSurvival) = _flatten_tuples(map(_survival_mass_thresholds, m.models))

# HardTemperatureLimit bound(s) as (u,t)->egg_temperature-relative closures,
# also root-found for the same reason. Always both bounds, even if one is
# +-Inf (never crosses, harmless) -- skipping it based on isfinite would make
# the condition count depend on parameter *values*, not just types, breaking
# the type-stability this whole rewrite is for.
function _temperature_conditions(m::HardTemperatureLimit, egg_model, pars, forcing)
    _temp(u, t) = begin
        unpacked = _unpack(u)
        (; environment) = forcing(t * u"hr")
        state = EggState(; unpacked...)
        egg_temperature(egg_model.thermal_model, egg_model, state, pars, environment)
    end
    (
        (u, t) -> ustrip(u"K", _temp(u, t) - m.lower_lethal_temperature),
        (u, t) -> ustrip(u"K", m.upper_lethal_temperature - _temp(u, t)),
    )
end
_temperature_conditions(::NoSurvivalLimit, egg_model, pars, forcing) = ()
_temperature_conditions(::DesiccationLimit, egg_model, pars, forcing) = ()
_temperature_conditions(m::CombinedSurvival, egg_model, pars, forcing) =
    _flatten_tuples(map(sub -> _temperature_conditions(sub, egg_model, pars, forcing), m.models))

# Only the soil-drier-than-egg check-valve is registered here, not the
# dry-mass floor or turgid-mass ceiling (also RHS kinks). The check-valve is
# genuinely transient (soil potential, an independent external forcing,
# crosses egg_water_potential and keeps moving past it). The floor/ceiling are
# not: mass can clamp and sit exactly on one of them for a long stretch (a dry
# spell, a long quiescence), so `u[2]-threshold` stays continuously zero
# rather than crossing once -- a ContinuousCallback root-finder can't tell
# "parked on a root" from "about to cross one" and re-fires every step,
# stalling progress. Same reasoning excludes hydration_index's clamp(...,0,1)
# corners and the ratchet's egg_mass==maximum_mass_achieved toggle -- all left
# as plain RHS clamps; Tsit5 integrates through a flat derivative region fine
# and only needs event detection for genuinely transient crossings.
_hydric_conditions(::AbstractHydricModel, pars, forcing) = ()
function _hydric_conditions(::SteadyDarcyFlux, pars, forcing)
    ((u, t) -> begin
        (; soil_water_potential) = forcing(t * u"hr")
        ustrip(u"J/kg", min(-0.001u"J/kg", soil_water_potential)) - u[3]
    end,)
end

# stage.jl's SteppedHydricStage switches conductance/wetness discontinuously
# at these development_fraction thresholds -- another undeclared RHS kink
# (same failure mode as _hydric_conditions above). Often numerically
# coincident with an arrest-window bound (e.g. both at dev=0.25) -- a
# redundant simultaneous root is harmless, CallbackSet handles it fine.
_hydric_stage_conditions(::AbstractHydricStageModel, pars, forcing) = ()
_hydric_stage_conditions(m::SteppedHydricStage, pars, forcing) =
    ((u, t) -> u[1] - m.conductance_threshold, (u, t) -> u[1] - m.wetness_threshold)

# every discontinuity in the RHS, as (u,t)->Float64 condition closures
# (root-crossing = event), as a Tuple -- fixed composition/length for a given
# egg_model type, so the CallbackSet built from it (see _callback_set) has the
# same compiled type on every call, letting solve() specialize once and reuse
# across every bout-chunk and every oviposition date instead of recompiling
# per call (a data-dependent condition subset would make solve()'s
# compilation dominate runtime -- 99%+ of wall time on a long run).
function _conditions(egg_model, pars, forcing)
    arrest = egg_model.arrest_model
    quiescence_conditions = _flatten_tuples(map(arrest.quiescence_windows) do window
        ((u, t) -> u[1] - window[1], (u, t) -> u[1] - window[2])
    end)
    survival_conditions = map(thr -> ((u, t) -> u[2] - thr), _survival_mass_thresholds(egg_model.survival_model))
    (
        (u, t) -> u[1] - 1.0,                                                # hatch
        (u, t) -> u[1] - arrest.diapause_window[1],
        (u, t) -> u[1] - arrest.diapause_window[2],
        quiescence_conditions...,
        # cold_hour_threshold IS registered: chill_accumulation keeps
        # increasing past it (only gates future diapause *eligibility*), so
        # this is a genuine one-time transient crossing.
        (u, t) -> u[4] - ustrip(u"hr", arrest.cold_hour_threshold),
        # diapause_hour_threshold is NOT registered, unlike cold_hour_threshold
        # above: diapause_duration's own accumulation is gated by comparing
        # itself to this same threshold (in_diapause requires
        # diapause_duration<=diapause_hour_threshold), so it self-terminates
        # and sits at the boundary once reached rather than crossing through
        # -- the same persistent-root pathology as the mass floor/ceiling.
        (u, t) -> _hydration_index_u(u, pars) - arrest.desiccation_tolerance,
        survival_conditions...,
        _temperature_conditions(egg_model.survival_model, egg_model, pars, forcing)...,
        _hydric_conditions(egg_model.hydric_model, pars, forcing)...,
        _hydric_stage_conditions(egg_model.hydric_stage_model, pars, forcing)...,
    )
end

# always the full fixed set of conditions, always fully armed -- built once
# per simulate_egg call (doesn't depend on u/t at all), reused unchanged
# across every chunk. Disarming a condition by checking whether its value is
# near zero at chunk start would be fragile: a condition sitting just outside
# that tolerance can still be crossed again almost immediately by genuine-but-
# tiny forward integration, producing many near-zero-duration chunks with no
# real progress. simulate_egg's mandatory grace-period time-advance after
# every chunk sidesteps this regardless of which condition is involved,
# instead of needing this function to guess a tolerance at all.
#
# recursive (not `for`/`map`) so it fully unrolls at compile time regardless
# of how many conditions there are -- same rationale as _flatten_tuples above.
@inline _fill_conditions!(out, ::Tuple{}, i, u, t) = out
@inline function _fill_conditions!(out, conditions::Tuple, i, u, t)
    out[i] = first(conditions)(u, t)
    _fill_conditions!(out, Base.tail(conditions), i + 1, u, t)
end

# One VectorContinuousCallback over all conditions, not a CallbackSet of one
# ContinuousCallback per condition. With ~15-20 conditions (arrest windows,
# quiescence windows, survival/temperature/hydric thresholds), a CallbackSet
# holds that many distinct ContinuousCallback *types* in one heterogeneous
# Tuple -- confirmed by allocation profiling that OrdinaryDiffEq's rootfinding
# falls back to a dynamically-typed path over that tuple (visible as
# Memory{Any} and individual model-struct allocations in the profile),
# reboxing captured state on every root search. A single VectorContinuousCallback
# is one concrete, homogeneous callback type regardless of condition count,
# which avoids that path entirely -- this is SciML's own documented fix for
# "many chained conditions" (see VectorContinuousCallback's docstring).
function _callback_set(conditions::Tuple)
    n = length(conditions)
    SciMLBase.VectorContinuousCallback(
        (out, u, t, integrator) -> _fill_conditions!(out, conditions, 1, u, t),
        (integrator, event_idx) -> SciMLBase.terminate!(integrator),
        n,
    )
end

# temperature isn't an ODE state (see thermal.jl), so trajectory recording
# recomputes it at each saved point from the saved state + forcing(t).
function _trajectory_point(egg_model, pars, forcing, t_hr, u)
    unpacked = _unpack(u)
    (; environment) = forcing(t_hr)
    state = EggState(; unpacked...)
    temperature = egg_temperature(egg_model.thermal_model, egg_model, state, pars, environment)
    (; t=t_hr, unpacked.development_fraction, unpacked.egg_mass, unpacked.egg_water_potential, temperature)
end

"""
    simulate_egg(egg_model, pars, initial_state, soil_hydraulics, forcing, tspan_hr; save_trajectory=false, saveat_hr=1.0)

Run the coupled development/hydric ODE from `tspan_hr[1]` to `tspan_hr[2]`
(hours), re-chunking at every arrest-condition threshold crossing. `forcing(t)`
returns `(; environment, soil_water_potential)` for time `t` (Unitful hours).
Survival (`egg_model.survival_model`) is checked at every `saveat_hr` point
regardless of `save_trajectory`, so a lethal excursion isn't missed inside a
long adaptive-step chunk. Returns `(; hatched::Bool, hatch_time, died::Bool,
death_time, death_cause::Symbol, final_state::EggState)`, plus `trajectory` (a
NamedTuple of vectors: `t`, `development_fraction`, `egg_mass`,
`egg_water_potential`, `temperature`) when `save_trajectory=true` — `nothing`
otherwise, since keeping the full history costs allocation the fast/grid-scale
path shouldn't pay for.
"""
function simulate_egg(egg_model::EggModel, pars::EggParameters, initial_state::EggState,
                      soil_hydraulics, forcing, tspan_hr; save_trajectory=false, saveat_hr=1.0)
    p = (; egg_model, pars, soil_hydraulics, forcing)
    conditions = _conditions(egg_model, pars, forcing)

    u = SVector(
        initial_state.development_fraction, ustrip(u"kg", initial_state.egg_mass),
        ustrip(u"J/kg", initial_state.egg_water_potential),
        ustrip(u"hr", initial_state.arrest_state.chill_accumulation),
        ustrip(u"hr", initial_state.arrest_state.diapause_duration),
        ustrip(u"kg", initial_state.maximum_mass_achieved),
    )
    t = ustrip(u"hr", tspan_hr[1])
    t_end = ustrip(u"hr", tspan_hr[2])

    trajectory_points = save_trajectory ? [_trajectory_point(egg_model, pars, forcing, t * u"hr", u)] : nothing
    died = false
    death_time = t * u"hr"
    death_cause = :alive

    # dtmax caps the adaptive step at the forcing's own resolution (1hr) --
    # otherwise Tsit5's error control can grow the step past a fast diurnal
    # oscillation and dense/saveat output just interpolates stale samples.
    # saveat only used when needed (survival check or trajectory) -- the
    # default fast/grid-scale path (NoSurvivalLimit) pays nothing extra.
    checking_survival = !(egg_model.survival_model isa NoSurvivalLimit)
    need_hourly_points = save_trajectory || checking_survival

    # built once -- doesn't depend on (u,t) at all (see _callback_set), so
    # every chunk (and every oviposition date, for a given egg_model) reuses
    # the exact same compiled solve() specialization.
    callback = _callback_set(conditions)

    # One integrator, built once and reused (via reinit!) across every
    # bout-chunk, instead of a fresh ODEProblem+solve() per chunk. A fresh
    # solve() re-merges/re-boxes the callback closures into a new CallbackSet
    # and rebuilds DEOptions on every call -- with ~100-200 chunks per egg and
    # 12 eggs per sweep, that dominated allocation (confirmed by profiling:
    # ~25% of bytes were callback-tuple reconstructions, not RHS evaluations,
    # which allocate a few hundred bytes each). reinit! resets the state/time
    # span on the existing integrator without rebuilding any of that.
    # need_hourly_points is fixed for the whole call, so the save behaviour
    # (saveat vs every-step) only needs configuring once, here.
    # qualified (like SciMLBase.successful_retcode/ReturnCode below) rather
    # than bare `init`/`reinit!`/`solve!` -- with the full package set
    # point_silo_deterministic.jl loads (Rasters/RasterDataSources/etc.),
    # another dependency also exports one of these generic-verb names, which
    # makes the bare name ambiguous/undefined rather than resolving to
    # SciMLBase's.
    problem = ODEProblem{false}(egg_rhs, u, (t, t_end), p)
    integrator = need_hourly_points ?
        SciMLBase.init(problem, Tsit5(); callback, saveat=saveat_hr, dtmax=1.0) :
        SciMLBase.init(problem, Tsit5(); callback, dtmax=1.0, save_everystep=false)
    # forced minimum time-advance after each chunk before the next one starts
    # -- see _callback_set's comment for why. Some conditions are self-
    # clamping accumulators (e.g. the mass floor/ceiling) that can settle at
    # their own root and stay there rather than cross it once; dt_grace backs
    # off (doubles) whenever a chunk makes near-zero real progress beyond the
    # grace nudge itself, and resets once real progress resumes -- a general
    # safety net against any such case, at the cost of a little precision
    # right at a genuine fast-but-real transition.
    dt_grace_min, dt_grace_max = 0.01, 4.0
    dt_grace = dt_grace_min

    chunk_count = 0
    while t < t_end
        chunk_count += 1
        if chunk_count % 200 == 0
            @info "simulate_egg: chunk $chunk_count, t=$(t)hr, u=$u, dt_grace=$(dt_grace)hr"
            @info "  condition values: $(map(cond -> cond(u, t), conditions))"
        end
        chunk_count > 50_000 && error(
            "simulate_egg: $(chunk_count) bout-chunks without reaching t_end=$(t_end)hr " *
            "(stuck at t=$(t)hr, u=$u, condition values=$(map(cond -> cond(u, t), conditions))) -- " *
            "likely a condition oscillating across its own threshold every chunk rather than " *
            "making forward progress.",
        )
        t_chunk_start = t
        SciMLBase.reinit!(integrator, u; t0=t, tf=t_end)
        SciMLBase.solve!(integrator)
        sol = integrator.sol
        # a chunk failing partway (e.g. dt forced below epsilon) still returns
        # a partial sol.u -- without this check that's silently treated as valid.
        if !SciMLBase.successful_retcode(sol) && sol.retcode != SciMLBase.ReturnCode.Terminated
            error("egg ODE solve failed at t=$(sol.t[end])hr: retcode=$(sol.retcode). u=$(sol.u[end])")
        end
        if need_hourly_points
            for i in 2:length(sol.t)
                t_i, u_i = sol.t[i], sol.u[i]
                point = _trajectory_point(egg_model, pars, forcing, t_i * u"hr", u_i)
                save_trajectory && push!(trajectory_points, point)
                if checking_survival
                    state_i = EggState(; _unpack(u_i)..., egg_water_potential=point.egg_water_potential)
                    if !survives(egg_model.survival_model, state_i, pars, point.temperature)
                        died, death_time, u, t = true, t_i * u"hr", u_i, t_i
                        death_cause = cause_of_death(egg_model.survival_model, state_i, pars, point.temperature)
                        break
                    end
                end
            end
            died && break
        end
        u, t = sol.u[end], sol.t[end]
        # root-finding stops a chunk *at* the crossing to within solver
        # tolerance (e.g. u[1]==0.9999999999999998, not exactly 1.0) -- a
        # strict >=1.0 here would miss that by float epsilon and spawn a
        # whole extra chunk.
        u[1] >= 1.0 - 1e-6 && break
        # nudge past the root this chunk just stopped at before starting the
        # next one -- see _callback_set for why. Adaptive backoff: a chunk
        # that barely got past its own grace nudge means we're still in a
        # stall -- widen it. Real progress (a chunk that ran meaningfully
        # longer than the nudge) resets it back down.
        chunk_duration = t - t_chunk_start
        if chunk_duration <= 2 * dt_grace
            # solve() made no real progress. This isn't always a genuine
            # zero-derivative equilibrium: e.g. development_fraction is
            # correctly frozen while halted (quiescent/diapausing), but a
            # callback condition on that frozen value alone (like the arrest
            # window's own entry boundary) then sits at a fixed near-zero
            # residual forever and keeps re-triggering solve() before any of
            # the *other*, genuinely-still-changing state (mass, water
            # potential) gets to move at all. Nudging t alone can't fix that
            # -- force one manual Euler step directly from the RHS so those
            # variables advance regardless.
            u = u + dt_grace * egg_rhs(u, p, t)
            dt_grace = min(dt_grace * 2, dt_grace_max)
        else
            dt_grace = dt_grace_min
        end
        t += dt_grace
    end

    unpacked = _unpack(u)
    hatched = !died && unpacked.development_fraction >= 1.0 - 1e-6
    # a death callback can overshoot the true zero-crossing by one step before
    # firing -- clamp for a physically sane reported final mass.
    final_state = EggState(; unpacked..., egg_mass=max(unpacked.egg_mass, 0.0u"kg"), hatched)
    trajectory = save_trajectory ? (;
        t=[p.t for p in trajectory_points],
        development_fraction=[p.development_fraction for p in trajectory_points],
        egg_mass=[p.egg_mass for p in trajectory_points],
        egg_water_potential=[p.egg_water_potential for p in trajectory_points],
        temperature=[p.temperature for p in trajectory_points],
    ) : nothing
    # hatch_time is always a concrete Quantity (not Union{Missing,Quantity}) --
    # meaningful only when hatched; check that flag, not this value, for validity.
    (; hatched, hatch_time=t * u"hr", died, death_time, death_cause, final_state, trajectory)
end
