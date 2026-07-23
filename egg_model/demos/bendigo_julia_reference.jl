# Julia trajectory for the same Bendigo Feb-1-2020, 150-day scenario as
# bendigo_r_reference.R, for a direct numerical comparison. Reuses the cached
# microclimate result and forcing exactly as point_silo_deterministic.jl.

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

cache_path = joinpath(@__DIR__, "microclimate_cache.jls")
result = deserialize(cache_path)

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
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), geometry,
)

function run_and_save(oviposition_date, out_name; duration=150.0u"d")
    start_hr = oviposition_offset(oviposition_date, dates)
    tspan = (start_hr, start_hr + duration)
    initial_state = EggState(;
        egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
        maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
    )

    println("Running Julia egg model, Bendigo, lay date $oviposition_date, $(duration)...")
    r = simulate_egg(egg_model, pars, initial_state, soil_hydraulics, forcing, tspan; save_trajectory=true, saveat_hr=1.0)

    outcome = r.hatched ? "hatched after $(round(typeof(1.0u"d"), r.hatch_time - start_hr, digits=1))" :
              r.died ? "died of $(r.death_cause) after $(round(typeof(1.0u"d"), r.death_time - start_hr, digits=1))" :
              "did not hatch in $(duration)"
    println("outcome: $outcome")
    println("final development_fraction: $(r.final_state.development_fraction)")
    println("final egg_mass: $(uconvert(u"mg", r.final_state.egg_mass))")

    traj = r.trajectory
    hours = ustrip.(uconvert.(u"hr", traj.t .- start_hr))
    mass_mg = ustrip.(uconvert.(u"mg", traj.egg_mass))
    psi_e = ustrip.(uconvert.(u"J/kg", traj.egg_water_potential))
    temperature_C = ustrip.(uconvert.(u"°C", traj.temperature))

    open(joinpath(@__DIR__, out_name), "w") do io
        println(io, "hour,dev,mass_mg,psi_e,temperature_C")
        for i in eachindex(hours)
            println(io, "$(hours[i]),$(traj.development_fraction[i]),$(mass_mg[i]),$(psi_e[i]),$(temperature_C[i])")
        end
    end
    println("wrote ", length(hours), " rows to $out_name")
    r
end

run_and_save(Date(2020, 2, 1), "bendigo_julia_reference_output.csv")
run_and_save(Date(2020, 9, 1), "bendigo_julia_spring_output.csv")
