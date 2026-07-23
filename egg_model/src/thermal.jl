using HeatExchange
using BiophysicalGeometry: Body, Naked, Ellipsoid
using Unitful

# Both modes are algebraic, not ODE states: an egg's thermal mass is tiny
# relative to the hourly timestep, so it equilibrates effectively instantly —
# no need to integrate a core_temperature state alongside development/hydrics.
# `environment` is the `(; environment_pars, environment_vars)` NamedTuple
# HeatExchange expects.

egg_temperature(::SoilTemperatureEquals, egg_model, state, pars, environment) =
    environment.environment_vars.substrate_temperature

# egg_model.geometry only holds the initial mass -- rebuild the shape here
# from the egg's current mass (density/axis ratios fixed) so surface area
# (driving both evaporation and liquid uptake) scales with it as the egg
# swells or shrinks.
function egg_organism(egg_model::EggModel, state::EggState, pars::EggParameters)
    geometry = Ellipsoid(
        state.egg_mass, egg_model.geometry.density,
        egg_model.geometry.axis_ratio_b, egg_model.geometry.axis_ratio_c,
    )
    body = Body(geometry, Naked())
    traits = example_heat_exchange_traits(;
        shape_pars=geometry,
        conduction_pars_external=ExternalConductionParameters(; conduction_fraction=pars.conduction_fraction),
        evaporation_pars=AnimalEvaporationParameters(;
            skin_wetness=stage_skin_wetness(egg_model.hydric_stage_model, state.development_fraction, pars),
            bare_skin_fraction=1.0, insulation_fraction=0.0,
        ),
        hydraulic_pars=hydraulic_pars(state, pars),
        metabolism_pars=MetabolismParameters(; model=metabolic_rate_function(egg_model.metabolic_model)),
    )
    Organism(body, traits)
end

# full conduction/convection/radiation/evaporation budget, root-found to equilibrium.
function egg_temperature(::FullHeatBudget, egg_model, state, pars, environment)
    organism = egg_organism(egg_model, state, pars)
    solve_temperature(organism, environment).core_temperature
end
