using Unitful
using BiophysicalGeometry: Ellipsoid
using HeatExchange: HydraulicParameters

# --- development rate ---

abstract type AbstractDevelopmentModel end

# Thin wrapper around a ThermalPhysiology rate model (e.g. SharpSchoolDEBModel).
# ThermalPhysiology's own fields are bare Float64 (unitless), so the physical
# rate unit is attached here rather than assumed.
struct RateModel{M,U} <: AbstractDevelopmentModel
    tpc::M
    rate_unit::U           # e.g. 1.0u"hr^-1"
end

development_rate(m::RateModel, temperature) = m.tpc(temperature) * m.rate_unit

# --- arrest: a named-condition graph, not a single mutually-exclusive phase ---
#
# EggState.arrest_state is a NamedTuple whose fields are model-defined (not
# fixed here) so a model can track as many named conditions as it needs
# (diapause, quiescence, or a richer graph with intermediate named states).
# Moisture branches threshold `hydration_index` (the egg's own mechanistic
# read of soil moisture); temperature branches threshold `environment`.

abstract type AbstractArrestModel end

function advance_arrest end   # (model, arrest_state, development_fraction, hydration_index, environment) -> NamedTuple
function arrested end          # (model, arrest_state) -> Bool

# Diapause + an arbitrary number of moisture-limited quiescence windows,
# independent of each other (both can be active at once).
Base.@kwdef struct ProportionWindowArrest{CT,W1,W2,CH,DH,DT} <: AbstractArrestModel
    cold_temperature::CT                  # below this, an hour counts toward chill_accumulation
    diapause_window::W1
    quiescence_windows::W2                 # Tuple of (lo,hi) pairs, any length -- must be a Tuple, not a Vector, for type-stable condition-building in phases.jl
    cold_hour_threshold::CH               # chill_accumulation above this ⟹ diapause no longer entered
    diapause_hour_threshold::DH           # diapause_duration above this ⟹ diapause self-terminates
    desiccation_tolerance::DT = 0.8
end

# chill_accumulation/diapause_duration are the only state this model needs —
# both are continuous accumulators (integrated alongside development_fraction
# in phases.jl), so in_diapause/diapause_eligible are derived, not stored:
# diapause is self-terminating once diapause_duration crosses its threshold.
initial_arrest_state(::ProportionWindowArrest) = (; chill_accumulation=0.0u"hr", diapause_duration=0.0u"hr")

quiescence_active(arrest::ProportionWindowArrest, development_fraction, hydration_index) =
    any(arrest.quiescence_windows) do window
        hydration_index < arrest.desiccation_tolerance && window[1] < development_fraction < window[2]
    end

in_diapause(arrest::ProportionWindowArrest, development_fraction, arrest_state) =
    arrest.diapause_window[1] < development_fraction < arrest.diapause_window[2] &&
    arrest_state.chill_accumulation <= arrest.cold_hour_threshold &&
    arrest_state.diapause_duration <= arrest.diapause_hour_threshold

arrested(arrest::ProportionWindowArrest, arrest_state, development_fraction, hydration_index) =
    in_diapause(arrest, development_fraction, arrest_state) ||
    quiescence_active(arrest, development_fraction, hydration_index)

# Two independent rate accumulators; the second only proceeds once the first
# crosses 1 (e.g. a chilling requirement gating embryonic development).
Base.@kwdef struct DualAccumulatorArrest{C,E} <: AbstractArrestModel
    first_model::C                        # AbstractDevelopmentModel
    second_model::E                       # AbstractDevelopmentModel
end

# --- hydric stage switching: shell permeability/wetness change with development ---

abstract type AbstractHydricStageModel end

# no switching -- hydraulic_conductance/skin_wetness stay at their pars values.
struct ConstantHydricStage <: AbstractHydricStageModel end
stage_hydraulic_conductance(::ConstantHydricStage, development_fraction, pars) = pars.hydraulic_conductance
stage_skin_wetness(::ConstantHydricStage, development_fraction, pars) = pars.skin_wetness

# shell impermeable/less wetted before its threshold, boosted after — e.g. an
# embryo's shell permeability increasing partway through development.
Base.@kwdef struct SteppedHydricStage{CT,WT,DC,AC,DW,AW} <: AbstractHydricStageModel
    conductance_threshold::CT
    wetness_threshold::WT
    dormant_conductance::DC
    active_conductance::AC
    dormant_wetness::DW
    active_wetness::AW
end
stage_hydraulic_conductance(m::SteppedHydricStage, development_fraction, pars) =
    development_fraction > m.conductance_threshold ? m.active_conductance : m.dormant_conductance
stage_skin_wetness(m::SteppedHydricStage, development_fraction, pars) =
    development_fraction > m.wetness_threshold ? m.active_wetness : m.dormant_wetness

# --- hydric (liquid + vapor water exchange) ---

abstract type AbstractHydricModel end

# skips the water-balance ODE — fastest option. Quiescence still responds to
# moisture via a direct threshold on soil water potential (cheap, no integration).
Base.@kwdef struct NoHydricExchange{P} <: AbstractHydricModel
    critical_water_potential::P
end

struct SteadyDarcyFlux <: AbstractHydricModel end

# optional transient soil-moisture correction (Tracy Appendix E), not the default.
Base.@kwdef struct TransientSoilCorrection{D} <: AbstractHydricModel
    soil_diffusivity::D
end

# --- thermal (egg temperature) ---

abstract type AbstractThermalModel end

# egg temperature == soil temperature at nest depth; fast, used for grid runs.
struct SoilTemperatureEquals <: AbstractThermalModel end

# full conduction/convection/radiation/evaporation budget via HeatExchange.onelump.
struct FullHeatBudget <: AbstractThermalModel end

# --- survival (mortality, optional) ---
# survives(model, state, pars, temperature) -> Bool, so a model can key off
# temperature (heat/cold limits) or hydric state (desiccation) alike. Composable
# via CombinedSurvival, same "independent conditions" idea as arrest.
# cause_of_death(...) -> Symbol identifies which criterion failed (:alive if
# none did), so a death outcome can be reported with its actual cause.

abstract type AbstractSurvivalModel end

# no mortality tracked -- development/hydrics never stop regardless of state.
struct NoSurvivalLimit <: AbstractSurvivalModel end
survives(::NoSurvivalLimit, state, pars, temperature) = true
cause_of_death(::NoSurvivalLimit, state, pars, temperature) = :alive

# instant death outside a fixed temperature range (either bound optional).
Base.@kwdef struct HardTemperatureLimit{TL,TH} <: AbstractSurvivalModel
    lower_lethal_temperature::TL = -Inf * u"K"
    upper_lethal_temperature::TH = Inf * u"K"
end
survives(m::HardTemperatureLimit, state, pars, temperature) =
    m.lower_lethal_temperature < temperature < m.upper_lethal_temperature
cause_of_death(m::HardTemperatureLimit, state, pars, temperature) =
    temperature <= m.lower_lethal_temperature ? :cold :
    temperature >= m.upper_lethal_temperature ? :heat : :alive

# death once water content : dry mass ratio falls to or below a critical value
# (e.g. 0.6 for plague locust eggs).
Base.@kwdef struct DesiccationLimit{DM,R} <: AbstractSurvivalModel
    dry_mass::DM
    critical_water_ratio::R = 0.6
end
survives(m::DesiccationLimit, state, pars, temperature) =
    (state.egg_mass - m.dry_mass) / m.dry_mass > m.critical_water_ratio
cause_of_death(m::DesiccationLimit, state, pars, temperature) =
    survives(m, state, pars, temperature) ? :alive : :desiccation

# any number of survival criteria evaluated together -- dies if any one fails.
struct CombinedSurvival{T<:Tuple} <: AbstractSurvivalModel
    models::T
end
CombinedSurvival(models...) = CombinedSurvival(models)
survives(m::CombinedSurvival, state, pars, temperature) =
    all(sub -> survives(sub, state, pars, temperature), m.models)
function cause_of_death(m::CombinedSurvival, state, pars, temperature)
    for sub in m.models
        cause = cause_of_death(sub, state, pars, temperature)
        cause !== :alive && return cause
    end
    :alive
end

# a ThermalPhysiology.AbstractTDTModel-based cumulative thermal death time
# could be added later as another AbstractSurvivalModel (same `survives`
# protocol, new dispatch), once real time-at-temperature tolerance data exists.

# --- metabolic heat (optional) ---

abstract type AbstractMetabolicModel end

# HeatExchange.metabolic_rate(::Nothing, mass, T) = 0.0u"W" already, so this
# resolves to no heat generation with zero new code in HeatExchange.jl.
struct NoMetabolicHeat <: AbstractMetabolicModel end

# Wraps a developmental-stage-keyed rate function (destined for BiologicalScaling.jl).
Base.@kwdef struct EmpiricalStageMetabolicHeat{F} <: AbstractMetabolicModel
    stage_rate::F                         # (development_fraction, mass, temperature) -> W
end

# a DEB-computed rate would plug in here as another AbstractMetabolicModel (not implemented).

metabolic_rate_function(::NoMetabolicHeat) = nothing
metabolic_rate_function(m::EmpiricalStageMetabolicHeat) = m.stage_rate

# --- top-level model config, mirrors MicroclimateMapper's MicroModel ---

Base.@kwdef struct EggModel{D<:AbstractDevelopmentModel,A<:AbstractArrestModel,
                             H<:AbstractHydricModel,HS<:AbstractHydricStageModel,
                             T<:AbstractThermalModel,S<:AbstractSurvivalModel,
                             M<:AbstractMetabolicModel,G<:Ellipsoid}
    development_model::D
    arrest_model::A
    hydric_model::H = SteadyDarcyFlux()
    hydric_stage_model::HS = ConstantHydricStage()
    thermal_model::T = SoilTemperatureEquals()
    survival_model::S = NoSurvivalLimit()
    metabolic_model::M = NoMetabolicHeat()
    geometry::G
end

# --- egg state (time-varying, one per egg/ensemble member) ---
# Immutable: every call site constructs a fresh instance, so this stack-allocates.

Base.@kwdef struct EggState{F,M,P,HI,AS}
    development_fraction::F = 0.0
    egg_mass::M
    egg_water_potential::P
    maximum_mass_achieved::M              # running max, ratchets up (never down)
    hydration_index::HI = 0.0
    hatched::Bool = false
    arrest_state::AS                      # model-defined NamedTuple, see initial_arrest_state/advance_arrest/arrested
end

# --- egg parameters (fixed, species-level) ---

Base.@kwdef struct EggParameters{HC,SH,CF,SW,IM,MM,MX}
    hydraulic_conductance::HC             # egg shell K_e, HeatExchange's K_skin equivalent
    specific_hydration::SH
    conduction_fraction::CF               # soil-contact fraction, HeatExchange's ExternalConductionParameters equivalent
    skin_wetness::SW                      # baseline, before stage-dependent switching in phases.jl
    initial_egg_mass::IM
    minimum_egg_mass::MM                  # hard floor (dry mass)
    # hard ceiling on uptake: shell/membrane tension opposes further water gain
    # once fully turgid, a counter-force the Darcy flux term itself lacks.
    # Placeholder (2x initial mass) -- needs a real turgid-mass estimate.
    maximum_egg_mass::MX = 2 * initial_egg_mass
end

# builds a HeatExchange.HydraulicParameters snapshot for the current state.
function hydraulic_pars(state::EggState, pars::EggParameters)
    HydraulicParameters(;
        water_potential=state.egg_water_potential,
        hydraulic_conductance=pars.hydraulic_conductance,
        specific_hydration=pars.specific_hydration,
    )
end
