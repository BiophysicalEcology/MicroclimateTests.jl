# Diagnostic: print the components of hydric_rate at a specific point along
# the Feb-1 trajectory (hour ~1600, where R shows strong uptake but Julia
# stays flat) to isolate exactly where the flux magnitude diverges from R.

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model
using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Rasters, RasterDataSources, PointDataSources
using DataInterpolations
using Dates, Unitful
using Serialization

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
nest_depth = 5.0u"cm"
dates = Date(2020, 1, 1):Day(1):Date(2020, 12, 31)
soil_profile = example_soil_profile(depths)
result = deserialize(joinpath(@__DIR__, "microclimate_cache.jls"))
nest_node = nearest_node(nest_depth, depths)
day_range = 1:size(result.soil_temperature, 1)
environment_pars = example_environment_pars()
forcing = egg_nest_forcing(result, day_range, nest_node, environment_pars)
soil_hydraulics = (;
    air_entry_potential   = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b            = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)

shape_b = 0.69 / 1.82
geometry = Ellipsoid(0.0036u"g", 1000.0u"kg/m^3", 1 / shape_b, 1 / shape_b)
arrest = ProportionWindowArrest(;
    cold_temperature=u"K"(0.0u"°C"), diapause_window=(0.45, 0.46),
    quiescence_windows=[(0.25, 0.30), (0.45, 0.50)],
    cold_hour_threshold=1000.0u"hr", diapause_hour_threshold=240.0u"hr",
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
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), geometry,
)

oviposition_date = Date(2020, 2, 1)
start_hr = oviposition_offset(oviposition_date, dates)
t_probe = start_hr + 1600.0u"hr"

(; environment, soil_water_potential) = forcing(t_probe)
println("soil_water_potential = ", soil_water_potential)
println("soil_temperature = ", environment.environment_vars.substrate_temperature)

state = EggState(;
    development_fraction=0.2506, egg_mass=2.62u"mg", egg_water_potential=-1750.0u"J/kg",
    maximum_mass_achieved=2.62u"mg", arrest_state=(; chill_accumulation=0.0u"hr", diapause_duration=0.0u"hr"),
)

organism = egg_organism(egg_model, state, pars)
total_area = BiophysicalGeometry.total_area(body(organism))
soil_contact_area = total_area * pars.conduction_fraction
air_exposed_area = total_area * (1 - pars.conduction_fraction)
println("total_area = ", uconvert(u"mm^2", total_area))
println("soil_contact_area = ", uconvert(u"mm^2", soil_contact_area))

conductance = stage_hydraulic_conductance(egg_model.hydric_stage_model, state.development_fraction, pars)
println("conductance (K_e) = ", conductance)

k_s = soil_hydraulic_conductivity(soil_hydraulics, soil_water_potential)
println("k_s = ", k_s)
println("air_entry_potential = ", soil_hydraulics.air_entry_potential)
println("campbell_b = ", soil_hydraulics.campbell_b)
println("saturated_conductivity = ", soil_hydraulics.saturated_conductivity)

shell_resistance = 1 / (soil_contact_area * conductance)
soil_resistance = 1 / (sqrt(2π * soil_contact_area) * k_s)
println("shell_resistance = ", shell_resistance)
println("soil_resistance = ", soil_resistance)

m_liquid = soil_liquid_flux(soil_water_potential, state.egg_water_potential, soil_contact_area, conductance, soil_hydraulics)
println("m_liquid = ", uconvert(u"mg/hr", m_liquid))

core_temperature = egg_temperature(egg_model.thermal_model, egg_model, state, pars, environment)
m_vapor = cutaneous_water_loss(organism, state, core_temperature, environment, air_exposed_area)
println("m_vapor = ", uconvert(u"mg/hr", m_vapor))

println("net d_mass = ", uconvert(u"mg/hr", m_liquid - m_vapor))
