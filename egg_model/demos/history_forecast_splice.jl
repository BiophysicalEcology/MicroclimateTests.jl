# Stage 8: history -> forecast splice. Solves the historical microclimate +
# egg model from an oviposition date up to "now" (an ACCESS-S2 issue date),
# then continues the egg simulation into each of several ACCESS-S2 ensemble
# members' forecast microclimates, producing a per-member hatch-date
# distribution instead of one deterministic value.
#
# No combined-forcing closure is needed: simulate_egg! already accepts any
# (initial_state, soil_hydraulics, forcing, tspan) tuple, so the splice is
# just two sequential simulate_egg! calls -- historical then forecast --
# threading the first call's final_state through as the second's
# initial_state. Each call's `forcing` has its own t=0 origin (the start of
# its own microclimate solve), so tspan for the forecast leg starts at 0,
# not at the historical run's "now" offset.
#
# Microclimate-level continuity: the forecast's initial soil_moisture is
# seeded from the historical run's actual profile at "now" (not a generic
# default), matching the plan's "carry forward final soil state as init="
# idea at the level MicroVectorProblem actually supports.

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_properties_model, example_soil_hydraulic_model
using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Rasters, RasterDataSources, PointDataSources
using Dates, Unitful, Statistics
using Serialization

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))
include(joinpath(@__DIR__, "..", "src", "access_s2.jl"))

# all cached/serialized run output goes here, not directly in demos/ -- one
# gitignore entry (egg_model/demos/output/) instead of per-file patterns.
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
nest_depth = 5.0u"cm"
nest_node = nearest_node(nest_depth, depths)
environment_pars = example_environment_pars()

# Just what egg_nest_forcing (forcing.jl) actually reads, plus soil_moisture
# (kept for a possible future moisture threshold, cheap to retain -- see
# forcing.jl's comment for why everything else here was dropped: unused, or
# inert under this egg model's SoilTemperatureEquals/SteadyDarcyFlux config).
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
)

bendigo = (144.2826718, -36.7590183)

const SANDY_LOAM = (air_entry=1.5u"J/kg", b=3.1, Ksat=7.2e-4u"kg*s/m^3", field_capacity=0.21, wilting_point=0.10)
function soil_profile_from_texture(texture::NamedTuple, depths;
    bulk_density=1.3u"Mg/m^3", mineral_density=2.560u"Mg/m^3",
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
soil_profile = soil_profile_from_texture(SANDY_LOAM, depths)
soil_hydraulics = (;
    air_entry_potential    = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b             = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)

# ── egg model, identical config to point_silo_deterministic.jl ──

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
    HardTemperatureLimit(; lower_lethal_temperature=u"K"(-5.0u"°C"), upper_lethal_temperature=u"K"(52.0u"°C")),
    DesiccationLimit(; dry_mass=0.1 * pars.initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

# ── historical leg: SILO, oviposition_date through issue_date ("now") ──
# Genuinely sub-yearly (Jan-Jun) -- this is exactly what the sub-yearly-run
# fix (MicroclimateMapper.jl PR #29) unblocked.

oviposition_date = Date(2024, 5, 1)
issue_date = Date(2024, 6, 1)   # ACCESS-S2 issue date -- "now"
historical_dates = oviposition_date:Day(1):issue_date

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

println("Solving historical (SILO) microclimate: $oviposition_date to $issue_date...")
historical_problem = MicroVectorProblem(;
    model = historical_model, points=[bendigo], dates=historical_dates, soil_profile,
    init = (; soil_moisture = fill(0.2, length(depths))),
)
@time historical_output = solve(historical_problem)

historical_day_range = 1:size(historical_output.soil_temperature[point=1], 1)
historical_result = (;
    soil_temperature          = collect(historical_output.soil_temperature[point=1]),
    soil_moisture             = collect(historical_output.soil_moisture[point=1]),
    soil_water_potential      = collect(historical_output.soil_water_potential[point=1]),
    soil_thermal_conductivity = collect(historical_output.soil_thermal_conductivity[point=1]),
    soil_humidity             = collect(historical_output.soil_humidity[point=1]),
)
historical_forcing = egg_nest_forcing(historical_result, historical_day_range, nest_node, environment_pars)

# soil_moisture profile at the very last historical hour -- seeds the
# forecast's init, for continuity across the splice (microclimate-level).
now_soil_moisture = collect(historical_output.soil_moisture[point=1, Ti=lastindex(historical_day_range)])

# ── historical egg-model leg: oviposition_date -> issue_date ──

historical_tspan = (0.0u"hr", length(historical_day_range) * 1.0u"hr")
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

println("Running historical egg-model leg...")
historical_egg_result = simulate_egg(egg_model, pars, initial_state, soil_hydraulics, historical_forcing, historical_tspan)

if historical_egg_result.hatched
    println("Hatched during the historical period, before the forecast issue date -- nothing to splice.")
elseif historical_egg_result.died
    println("Died during the historical period ($(historical_egg_result.death_cause)) -- nothing to splice.")
else
    println("Still developing at issue date $issue_date (development_fraction=$(historical_egg_result.final_state.development_fraction)) -- continuing into the forecast ensemble.")

    # ── forecast leg: ACCESS-S2, issue_date + N members ──
    # Full 99-member ensemble. Each member's outcome is cached to disk
    # individually -- solving 99 separate ACCESS-S2 microclimates is the
    # dominant cost, so an interrupted run shouldn't have to redo already-
    # finished members. The egg-model integrator itself (init_egg_cache/
    # simulate_egg!) is built once from member 1's forcing and reused across
    # all members via reinit! rather than rebuilt per member (phases.jl's
    # standard cache-reuse pattern -- see points_australia.jl).

    members = 1:99
    forecast_horizon_days = 214
    forecast_dates = issue_date:Day(1):(issue_date + Day(forecast_horizon_days - 1))
    forecast_tspan = (0.0u"hr", length(forecast_dates) * 24.0u"hr")   # hours, not days

    egg_cache = nothing
    progress = 0
    report_every = max(1, length(members) ÷ 20)

    forecast_outcomes = map(members) do member
        member_cache_path = joinpath(output_dir,
            "splice_member$(member)_issue$(issue_date)_lay$(oviposition_date).jls")
        outcome = if isfile(member_cache_path)
            deserialize(member_cache_path)
        else
            forecast_model = MicroMapModel(;
                micro_model = MicroModel(;
                    depths, heights,
                    soil_properties_model = example_soil_properties_model(),
                    soil_hydraulic_model  = example_soil_hydraulic_model(),
                    snow_model            = NoSnow(),
                    config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
                ),
                dem_source              = CRUCL2,
                weather_source          = AccessS2(issue_date, member),
                surface_albedo_source   = 0.15,
                roughness_height_source = 0.004u"m",
                compute_terrain         = false,
                output_layers,
            )
            forecast_problem = MicroVectorProblem(;
                model = forecast_model, points=[bendigo], dates=forecast_dates, soil_profile,
                init = (; soil_moisture = now_soil_moisture),
            )
            forecast_output = solve(forecast_problem)
            forecast_day_range = 1:size(forecast_output.soil_temperature[point=1], 1)
            forecast_result = (;
                soil_temperature          = collect(forecast_output.soil_temperature[point=1]),
                soil_moisture             = collect(forecast_output.soil_moisture[point=1]),
                soil_water_potential      = collect(forecast_output.soil_water_potential[point=1]),
                soil_thermal_conductivity = collect(forecast_output.soil_thermal_conductivity[point=1]),
                soil_humidity             = collect(forecast_output.soil_humidity[point=1]),
            )
            forecast_forcing = egg_nest_forcing(forecast_result, forecast_day_range, nest_node, environment_pars)

            if egg_cache === nothing
                global egg_cache = init_egg_cache(egg_model, pars, historical_egg_result.final_state,
                    soil_hydraulics, forecast_forcing, forecast_tspan)
            end
            result = simulate_egg!(egg_cache, historical_egg_result.final_state,
                soil_hydraulics, forecast_forcing, forecast_tspan)
            serialize(member_cache_path, result)
            result
        end
        global progress += 1
        (progress % report_every == 0 || progress == length(members)) &&
            println("  $progress/$(length(members)) ensemble members done")
        outcome
    end

    println("\nPer-member outcomes (lay date $oviposition_date, issue date $issue_date):")
    hatch_dates = Date[]
    for (member, r) in zip(members, forecast_outcomes)
        if r.hatched
            hatch_date = issue_date + Day(round(Int, ustrip(u"d", r.hatch_time)))
            push!(hatch_dates, hatch_date)
            println("  member $member -> hatched on $hatch_date")
        elseif r.died
            println("  member $member -> died of $(r.death_cause)")
        else
            println("  member $member -> did not hatch within $forecast_horizon_days days of the forecast")
        end
    end
    if !isempty(hatch_dates)
        hatch_days_since_issue = [Dates.value(d - issue_date) for d in hatch_dates]
        println("\n$(length(hatch_dates))/$(length(members)) members hatched. ",
                "Median hatch: $(issue_date + Day(round(Int, median(hatch_days_since_issue)))), ",
                "range: $(minimum(hatch_dates)) to $(maximum(hatch_dates))")
    end
end
