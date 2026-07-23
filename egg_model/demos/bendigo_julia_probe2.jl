# Diagnostic: call egg_rhs directly at the exact frozen state reported by the
# Feb-1 2020 run (t=2373.267643380025hr, u=[0.2500000886745846, 3.6e-6,
# -709.4554743596425, 0.0, 9.094322041162513, 3.5999999999999994e-6]) to see
# why d_mass is apparently exactly zero there despite soil being wetter than
# the egg and conductance nominally active (dev > conductance_threshold=0.25).

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
    cold_temperature=u"K"(0.0u"°C"), diapause_window=(0.25, 0.30),
    quiescence_windows=((0.25, 0.30), (0.45, 0.50)),
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
survival_model = CombinedSurvival(
    HardTemperatureLimit(; lower_lethal_temperature=u"K"(-5.0u"°C"), upper_lethal_temperature=u"K"(55.0u"°C")),
    DesiccationLimit(; dry_mass=0.1 * pars.initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

# exact frozen state from the Feb-1 2020 run
u = [0.2500000886745846, 3.6e-6, -709.4554743596425, 0.0, 9.094322041162513, 3.5999999999999994e-6]
t = 2373.267643380025

p = (; egg_model, pars, soil_hydraulics, forcing)
(; development_fraction, egg_mass, egg_water_potential, arrest_state, maximum_mass_achieved) = _unpack(u)
(; environment, soil_water_potential) = forcing(t * u"hr")
state = EggState(; development_fraction, egg_mass, egg_water_potential, maximum_mass_achieved, arrest_state)

hyd_index = hydration_index(egg_model.hydric_model, state, pars, soil_water_potential)
halted = arrested(egg_model.arrest_model, arrest_state, development_fraction, hyd_index)
println("development_fraction = ", development_fraction)
println("hydration_index = ", hyd_index)
println("halted (arrested) = ", halted)
println("in_diapause = ", in_diapause(egg_model.arrest_model, development_fraction, arrest_state))
println("quiescence_active = ", quiescence_active(egg_model.arrest_model, development_fraction, hyd_index))
println("soil_water_potential = ", soil_water_potential)
println("egg_water_potential = ", egg_water_potential)

conductance = stage_hydraulic_conductance(egg_model.hydric_stage_model, development_fraction, pars)
println("stage conductance (pre check-valve) = ", conductance)
println("check-valve blocks? ", min(-0.001u"J/kg", soil_water_potential) < egg_water_potential)

(; d_mass, d_water_potential) = hydric_rate(
    egg_model.hydric_model, egg_model, state, pars, environment, soil_water_potential, soil_hydraulics,
)
println("d_mass = ", uconvert(u"mg/hr", d_mass))
println("d_water_potential = ", d_water_potential)

rhs = egg_rhs(u, p, t)
println("full RHS = ", rhs)
