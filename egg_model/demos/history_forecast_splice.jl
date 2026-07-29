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

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"#"z:/"

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))
include(joinpath(@__DIR__, "..", "src", "access_s2.jl"))

# load parameters
include(joinpath(@__DIR__, "..", "params", "chortoicetes.jl"))

# all cached/serialized run output goes here, not directly in demos/ -- one
# gitignore entry (egg_model/demos/output/) instead of per-file patterns.
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

site = geocode("Bendigo, Australia")
points = [site]

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

ensembles = 99
save_trajectory = true
diapause = true
oviposition_date = Date(2026, 4, 25)

if diapause
    nest_depth = 5.0u"cm"
    cold_hour_threshold = 30u"d"
    diapause_hour_threshold = 0.0u"d"
else
    nest_depth = 10.0u"cm"
    cold_hour_threshold = 0.0u"d"
    diapause_hour_threshold = 0.0u"d"
end
nest_node = nearest_node(nest_depth, depths)

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

environment_pars = example_environment_pars()
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
    DesiccationLimit(; dry_mass=0.1 * pars.initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

# ── historical leg: SILO, oviposition_date through issue_date ("now") ──

issue_date = Date(2026, 7, 1)   # ACCESS-S2 issue date -- "now"
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
    model = historical_model, points, dates=historical_dates, soil_profile,
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

# ── historical egg-model leg: oviposition_date -> issue_date ──

historical_tspan = (0.0u"hr", length(historical_day_range) * 1.0u"hr")
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

println("Running historical egg-model portion...")
historical_egg_result = simulate_egg(egg_model, pars, initial_state, soil_hydraulics, 
                                        historical_forcing, historical_tspan; save_trajectory)

# soil temperature and moistures profile at the very last historical hour -- seeds the
# forecast's init, for continuity across the splice (microclimate-level).
now_soil_temperature = collect(historical_output.soil_temperature[point=1, Ti=lastindex(historical_day_range)])
now_soil_moisture = collect(historical_output.soil_moisture[point=1, Ti=lastindex(historical_day_range)])

label = ""
if historical_egg_result.hatched
    println("Hatched during the historical period, before the forecast issue date -- nothing to splice.")
elseif historical_egg_result.died
    println("Died during the historical period ($(historical_egg_result.death_cause)) -- nothing to splice.")
else
    println("Still developing at issue date $issue_date (development_fraction=$(historical_egg_result.final_state.development_fraction)) -- continuing into the forecast ensemble.")

    # ── forecast leg: ACCESS-S2, issue_date + N members ──
    # doing user-specified number of the 99-member ensemble. 
    # init_egg_cache/simulate_egg!) is built once from member 1's forcing and reused across
    # all members via reinit!.

    members = 1:ensembles
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
                    config                = MicroConfig(
                                            soil_moisture_strategy = DynamicSoilMoisture(),
                                            convergence = FixedSoilTemperatureIterations(1)
                                            ),
                ),
                dem_source              = CRUCL2,
                weather_source          = AccessS2(issue_date, member),
                surface_albedo_source   = 0.15,
                roughness_height_source = 0.004u"m",
                compute_terrain         = false,
                output_layers,
            )
            forecast_problem = MicroVectorProblem(;
                model = forecast_model, points, dates=forecast_dates, soil_profile,
                init = (; soil_moisture = now_soil_moisture, soil_temperature = now_soil_temperature),
            )
            println("simulate microclimate")
            @time forecast_output = solve(forecast_problem)
            forecast_day_range = 1:size(forecast_output.soil_temperature[point=1], 1)
            forecast_result = (;
                soil_temperature          = collect(forecast_output.soil_temperature[point=1]),
                soil_moisture             = collect(forecast_output.soil_moisture[point=1]),
                soil_water_potential      = collect(forecast_output.soil_water_potential[point=1]),
                soil_thermal_conductivity = collect(forecast_output.soil_thermal_conductivity[point=1]),
                soil_humidity             = collect(forecast_output.soil_humidity[point=1]),
            )
            # # Bias-correct soil water potential at the nest depth to match the
            # # historical leg's endpoint, decaying to zero over the first day.
            # swp_bias = historical_result.soil_water_potential[end, nest_node] -
            #     forecast_result.soil_water_potential[2, nest_node]
            # temp_bias = historical_result.soil_temperature[end, nest_node] -
            #     forecast_result.soil_temperature[2, nest_node]                
            # println(swp_bias)
            # println(temp_bias)
            # blend_hours = 24
            # n_forecast = size(forecast_result.soil_water_potential, 1)
            # decay = clamp.(1.0 .- (0:n_forecast-1) ./ blend_hours, 0.0, 1.0)
            # forecast_result.soil_water_potential[:, nest_node] .+= swp_bias .* decay
            # forecast_result.soil_temperature[:, nest_node] .+= temp_bias .* decay

            forecast_forcing = egg_nest_forcing(forecast_result, forecast_day_range, nest_node, environment_pars)

            if egg_cache === nothing
                global egg_cache = init_egg_cache(egg_model, pars, historical_egg_result.final_state,
                    soil_hydraulics, forecast_forcing, forecast_tspan; save_trajectory)
            end
            result = simulate_egg!(egg_cache, historical_egg_result.final_state,
                soil_hydraulics, forecast_forcing, forecast_tspan)
            #serialize(member_cache_path, result)
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
        label = "$(length(hatch_dates))/$(length(members)) members hatched. " *
        "Median hatch: $(issue_date + Day(round(Int, median(hatch_days_since_issue)))), " *
        "range: $(minimum(hatch_dates)) to $(maximum(hatch_dates))"
        println(label)
    end
end


# trajectories for every oviposition date, overlaid on the same panels.
using Plots

p1 = plot(;
    ylabel = "development",
    title = split(site.display_name, ",")[1] *
            ", lay date: " * string(oviposition_date) *
            ", diapause = " * string(diapause),
    legend = false,
    ylims = (0, 1)
)
p2 = plot(; ylabel="egg mass", legend=false, ylims = (minimum_egg_mass, maximum_egg_mass))
p3 = plot(; ylabel="egg\ntemperature", xlabel="date", legend=false)

silo_t = ustrip.(u"s", historical_egg_result.trajectory.t)

for i in 1:ensembles
    traj = forecast_outcomes[i].trajectory
    access_t = ustrip.(u"s", traj.t) .+ last(silo_t)
    combined_t = vcat(silo_t, access_t)
    development_fraction = vcat(
        historical_egg_result.trajectory.development_fraction,
        traj.development_fraction
    )
    egg_mass = vcat(
        historical_egg_result.trajectory.egg_mass,
        traj.egg_mass
    )
    egg_water_potential = vcat(
        historical_egg_result.trajectory.egg_water_potential,
        traj.egg_water_potential
    )
    temperature = vcat(
        historical_egg_result.trajectory.temperature,
        traj.temperature
    )
    actual_times =
        DateTime(oviposition_date) .+
        Dates.Second.(round.(Int, combined_t))
    plot!(p1, actual_times, development_fraction)
    plot!(
        p2,
        actual_times,
        collect(uconvert.(u"mg", egg_mass));
        ylims = (0.0u"mg", pars.initial_egg_mass * 2.0)
    )
    plot!(
        p3,
        actual_times,
        collect(uconvert.(u"°C", temperature))
    )
end

label = if isempty(hatch_dates)
    "0/$(length(members)) members hatched.\n" *
    "Median hatch: none\n" *
    "Range: none"
else
    "$(length(hatch_dates))/$(length(members)) members hatched.\n" *
    "Median hatch: $(issue_date + Day(round(Int, median(hatch_days_since_issue)))),\n" *
    "Range: $(minimum(hatch_dates)) to $(maximum(hatch_dates))"
end

annotate!(
    p1,
    DateTime(oviposition_date),
    0.85,
    text(label, 10, :black, :left)
)

combined_plot = plot(
    p1, p2, p3;
    layout = (3, 1),
    size = (1000, 900),
    link = :x
)

savefig(
    combined_plot,
    joinpath(@__DIR__, "history_forecast_splice.png")
)

display(combined_plot)