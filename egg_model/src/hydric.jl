using HeatExchange
using BiophysicalGeometry
using Unitful

# no water balance at all -- fastest option; egg_mass/egg_water_potential don't move.
hydric_rate(::NoHydricExchange, egg_model, state, pars, environment, soil_water_potential, soil_hydraulics) =
    (d_mass=0.0u"kg/s", d_water_potential=0.0u"J/kg/s")

# direct moisture proxy, no integration: soil drier than the threshold ⟹ desiccated.
hydration_index(hydric::NoHydricExchange, state, pars, soil_water_potential) =
    soil_water_potential > hydric.critical_water_potential ? 1.0 : 0.0

# relative to the highest mass ever achieved (dynamic, ratchets up), not a
# fixed ceiling, so it tracks recovery from a dry spell relative to the egg's
# own history rather than an arbitrary absolute target.
hydration_index(::SteadyDarcyFlux, state, pars, soil_water_potential) =
    clamp((state.egg_mass - pars.initial_egg_mass) / (state.maximum_mass_achieved - pars.initial_egg_mass), 0.0, 1.0)

# soil-side unsaturated hydraulic conductivity, Campbell & Norman (2001) eq. 9.2.
# air_entry_potential is never positive -- some soil-profile sources (e.g.
# Microclimate.jl's example_soil_profile) store it as a raw magnitude (e.g.
# +0.7 J/kg) rather than the physical (negative) value, so sign-correct by
# negating rather than clamping towards zero (clamping would collapse a
# realistic ~0.5-0.7 J/kg magnitude to near-zero and send k_s -> ~0,
# switching off liquid uptake almost entirely). The tiny floor is only a
# fallback for the degenerate raw==0 case, avoiding a literal 0/0 at saturation.
function soil_hydraulic_conductivity(soil_hydraulics, soil_water_potential)
    floor_potential = -1e-6 * oneunit(soil_hydraulics.air_entry_potential)
    air_entry_potential = min(-abs(soil_hydraulics.air_entry_potential), floor_potential)
    clamped_potential = min(soil_water_potential, air_entry_potential)
    soil_hydraulics.saturated_conductivity *
        (air_entry_potential / clamped_potential)^(2 + 3 / soil_hydraulics.campbell_b)
end

# Darcy resistances in series: egg shell resistance + spherical soil-flow resistance.
function soil_liquid_flux(soil_water_potential, egg_water_potential, soil_contact_area,
                           hydraulic_conductance, soil_hydraulics)
    k_s = soil_hydraulic_conductivity(soil_hydraulics, soil_water_potential)
    shell_resistance = 1 / (soil_contact_area * hydraulic_conductance)
    soil_resistance = 1 / (sqrt(2π * soil_contact_area) * k_s)
    (soil_water_potential - egg_water_potential) / (shell_resistance + soil_resistance)
end

# vapor-phase cutaneous water loss, reusing HeatExchange's evaporation physics
# directly (same free-convection Nusselt/Sherwood correlations as the source model).
# Goes through the organism's own accessors (as heat_balance itself does) so
# ModelParameters.Param-wrapped defaults are stripped before use. Qualified with
# HeatExchange. -- Microclimate.jl also exports convection/evaporation (soil
# surface physics), so the bare names are ambiguous with both loaded.
function cutaneous_water_loss(organism, state, core_temperature, environment, air_exposed_area)
    (; environment_vars) = environment
    conv = HeatExchange.convection(;
        body=body(organism), area=air_exposed_area,
        air_temperature=environment_vars.air_temperature,
        surface_temperature=core_temperature,
        wind_speed=environment_vars.wind_speed,
        atmospheric_pressure=environment_vars.atmospheric_pressure,
        fluid=Air(),
    )
    atmos = AtmosphericConditions(environment_vars)
    HeatExchange.evaporation(
        evaporation_pars(organism), conv.mass_transfer_coefficient, atmos, air_exposed_area,
        core_temperature, environment_vars.air_temperature;
        water_potential=state.egg_water_potential,
    ).cutaneous_mass_flow
end

# combined liquid uptake + vapor loss -> (d_mass/dt, d_water_potential/dt).
function hydric_rate(::SteadyDarcyFlux, egg_model, state, pars, environment, soil_water_potential, soil_hydraulics)
    organism = egg_organism(egg_model, state, pars)
    core_temperature = egg_temperature(egg_model.thermal_model, egg_model, state, pars, environment)
    total_area = BiophysicalGeometry.total_area(body(organism))
    soil_contact_area = total_area * pars.conduction_fraction
    air_exposed_area = total_area * (1 - pars.conduction_fraction)

    conductance = stage_hydraulic_conductance(egg_model.hydric_stage_model, state.development_fraction, pars)
    # check-valve: liquid uptake is soil->egg only. Once soil dries out below
    # the egg's own water potential, block this pathway entirely rather than
    # letting the (symmetric) Darcy formula run it in reverse -- otherwise a
    # wet spell's gains keep leaking back out through this term as the soil
    # dries again, on top of ordinary evaporative loss, and the egg can never
    # net-accumulate water.
    if min(-0.001u"J/kg", soil_water_potential) < state.egg_water_potential
        conductance = zero(conductance)
    end
    m_liquid = soil_liquid_flux(
        soil_water_potential, state.egg_water_potential,
        soil_contact_area, conductance, soil_hydraulics,
    )
    m_vapor = cutaneous_water_loss(organism, state, core_temperature, environment, air_exposed_area)

    d_mass = m_liquid - m_vapor
    # dry-mass floor: egg_mass can't drop below minimum_egg_mass (the egg's dry
    # mass) -- clamp the flux itself once there, not just the denominator below,
    # or the continuous ODE state drifts through and past it (even negative).
    if state.egg_mass <= pars.minimum_egg_mass && d_mass < zero(d_mass)
        d_mass = zero(d_mass)
    end
    # turgid-mass ceiling: shell/membrane tension resists further uptake once
    # fully hydrated -- the Darcy flux term alone has no such counter-force,
    # and uptake vs. surface-area growth is a positive feedback with no other
    # brake, so egg mass would otherwise run away unboundedly.
    if state.egg_mass >= pars.maximum_egg_mass && d_mass > zero(d_mass)
        d_mass = zero(d_mass)
    end
    floored_mass = max(state.egg_mass, pars.minimum_egg_mass)
    d_water_potential = d_mass / (pars.specific_hydration * floored_mass)
    (; d_mass, d_water_potential)
end
