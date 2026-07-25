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
points = [
    geocode("Bendigo, Australia"),
    geocode("Mildura, Australia"),
    geocode("Hay, NSW, Australia"),
    geocode("Ivanhoe, NSW, Australia"),
    geocode("Dubbo, NSW, Australia"),
    geocode("Bourke, NSW, Australia"),
    geocode("St George, Qld, Australia"),
    geocode("Nooyeah Downs, Qld, Australia"),
    geocode("Quilpie, Qld, Australia"),
    geocode("Roma, Qld, Australia"),
    geocode("Marree, SA, Australia"),
    geocode("Broken Hill, Australia"),
    geocode("Hawker, SA, Australia"),
    geocode("Yunta, SA, Australia"),
    geocode("Birdsville, Qld, Australia"),
]

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
    20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

ensembles = 99
save_trajectory = true
diapause = false
oviposition_date = Date(2026, 4, 25)

if diapause
    nest_depth = 5.0u"cm"
    cold_hour_threshold = 30u"d"
else
    nest_depth = 10.0u"cm"
    cold_hour_threshold = 0.0u"d"
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
    air_entry_potential=soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity=soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b=soil_profile.hydraulics.campbell_b_parameter[nest_node],
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
    initial_egg_mass, minimum_egg_mass, maximum_egg_mass,
)
survival_model = CombinedSurvival(
    HardTemperatureLimit(; lower_lethal_temperature, upper_lethal_temperature),
    DesiccationLimit(; dry_mass, critical_water_ratio),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential,
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

# ─────────────────────────────────────────────────────────────────────────────
# Run all sites together so MicroVectorProblem can parallelise over `points`.
# The microclimate model is solved once for all sites historically and once for
# all sites for each ACCESS-S2 ensemble member.
# ─────────────────────────────────────────────────────────────────────────────

issue_date = Date(2026, 7, 1)
historical_dates = oviposition_date:Day(1):issue_date
members = collect(1:ensembles)
nsites = length(points)

site_names = [split(site.display_name, ",")[1] for site in points]
site_slugs = [replace(strip(name), r"[^A-Za-z0-9]+" => "_") for name in site_names]

mkpath(output_dir)

# ── historical microclimate: all 15 sites in one solve ───────────────────────

historical_model = MicroMapModel(;
    micro_model=MicroModel(;
        depths, heights,
        soil_properties_model=example_soil_properties_model(),
        soil_hydraulic_model=example_soil_hydraulic_model(),
        snow_model=NoSnow(),
        config=MicroConfig(
            soil_moisture_strategy=DynamicSoilMoisture(),
        ),
    ),
    dem_source=CRUCL2,
    weather_source=SILO,
    surface_albedo_source=0.15,
    roughness_height_source=0.004u"m",
    compute_terrain=false,
    output_layers,
)

println("Solving historical SILO microclimate for $nsites sites: " *
        "$oviposition_date to $issue_date...")

historical_problem = MicroVectorProblem(;
    model=historical_model,
    points,                       # all sites passed together
    dates=historical_dates,
    soil_profile,
    init=(; soil_moisture=fill(0.2, length(depths))),
)

@time historical_output = solve(historical_problem)

# Run one historical egg simulation per site from the multi-site output.
historical_egg_results = Vector{Any}(undef, nsites)
historical_day_ranges = Vector{UnitRange{Int}}(undef, nsites)

for site_index in eachindex(points)
    site_name = site_names[site_index]
    println("Running historical egg model for $site_name...")

    historical_day_range =
        1:size(historical_output.soil_temperature[point=site_index], 1)
    historical_day_ranges[site_index] = historical_day_range

    historical_result = (;
        soil_temperature=collect(
            historical_output.soil_temperature[point=site_index]
        ),
        soil_moisture=collect(
            historical_output.soil_moisture[point=site_index]
        ),
        soil_water_potential=collect(
            historical_output.soil_water_potential[point=site_index]
        ),
        soil_thermal_conductivity=collect(
            historical_output.soil_thermal_conductivity[point=site_index]
        ),
        soil_humidity=collect(
            historical_output.soil_humidity[point=site_index]
        ),
    )

    historical_forcing = egg_nest_forcing(
        historical_result,
        historical_day_range,
        nest_node,
        environment_pars,
    )

    historical_tspan =
        (0.0u"hr", length(historical_day_range) * 1.0u"hr")

    site_initial_state = EggState(;
        egg_mass=pars.initial_egg_mass,
        egg_water_potential,
        maximum_mass_achieved=pars.initial_egg_mass,
        arrest_state=initial_arrest_state(arrest),
    )

    historical_egg_results[site_index] = simulate_egg(
        egg_model,
        pars,
        site_initial_state,
        soil_hydraulics,
        historical_forcing,
        historical_tspan;
        save_trajectory,
    )
end

# MicroVectorProblem currently requires init.soil_moisture and
# init.soil_temperature to each be ONE depth vector shared by all points; it
# does not accept a point × depth matrix. Preserve the all-sites threaded solve
# by averaging the site-specific historical endpoint profiles at each depth.
last_historical_hour = size(historical_output.soil_temperature, 1)

site_final_soil_temperature = [
    collect(historical_output.soil_temperature[
        point=site_index,
        Ti=last_historical_hour,
    ])
    for site_index in eachindex(points)
]

site_final_soil_moisture = [
    collect(historical_output.soil_moisture[
        point=site_index,
        Ti=last_historical_hour,
    ])
    for site_index in eachindex(points)
]

# These are vectors of length(depths), as required by MicroclimateMapper.
now_soil_temperature = reduce(+, site_final_soil_temperature) ./ nsites
now_soil_moisture = reduce(+, site_final_soil_moisture) ./ nsites

@assert now_soil_temperature isa AbstractVector
@assert now_soil_moisture isa AbstractVector
@assert length(now_soil_temperature) == length(depths)
@assert length(now_soil_moisture) == length(depths)

# Identify sites that still require a forecast egg simulation.
forecast_site_indices = [
    i for i in eachindex(points)
          if !historical_egg_results[i].hatched && !historical_egg_results[i].died
]

for site_index in eachindex(points)
    r = historical_egg_results[site_index]
    if r.hatched
        println("$(site_names[site_index]): hatched before forecast issue date.")
    elseif r.died
        println("$(site_names[site_index]): died historically ($(r.death_cause)).")
    else
        println("$(site_names[site_index]): still developing at $issue_date " *
                "(development_fraction=$(r.final_state.development_fraction)).")
    end
end

forecast_horizon_days = 214
forecast_dates = issue_date:Day(1):(issue_date+Day(forecast_horizon_days-1))
forecast_tspan = (0.0u"hr", length(forecast_dates) * 24.0u"hr")

# forecast_outcomes[site_index, member_index]. Entries remain `nothing` for
# sites that had already hatched or died during the historical leg.
forecast_outcomes = Matrix{Any}(nothing, nsites, length(members))
egg_caches = Vector{Any}(nothing, nsites)

if !isempty(forecast_site_indices)
    for (member_index, member) in pairs(members)
        println("\nACCESS-S2 member $member/$(length(members)): " *
                "solving microclimate for all $nsites sites...")

        forecast_model = MicroMapModel(;
            micro_model=MicroModel(;
                depths, heights,
                soil_properties_model=example_soil_properties_model(),
                soil_hydraulic_model=example_soil_hydraulic_model(),
                snow_model=NoSnow(),
                config=MicroConfig(
                    soil_moisture_strategy=DynamicSoilMoisture(),
                    convergence=FixedSoilTemperatureIterations(1),
                ),
            ),
            dem_source=CRUCL2,
            weather_source=AccessS2(issue_date, member),
            surface_albedo_source=0.15,
            roughness_height_source=0.004u"m",
            compute_terrain=false,
            output_layers,
        )

        forecast_problem = MicroVectorProblem(;
            model=forecast_model,
            points,                   # all 15 sites passed together
            dates=forecast_dates,
            soil_profile,
            init=(;
                soil_moisture=now_soil_moisture,
                soil_temperature=now_soil_temperature,
            ),
        )

        @time forecast_output = solve(forecast_problem)

        # The egg simulations are light relative to the microclimate solve.
        # Each site has its own cache because its historical final state differs.
        for site_index in forecast_site_indices
            forecast_day_range =
                1:size(forecast_output.soil_temperature[point=site_index], 1)

            forecast_result = (;
                soil_temperature=collect(
                    forecast_output.soil_temperature[point=site_index]
                ),
                soil_moisture=collect(
                    forecast_output.soil_moisture[point=site_index]
                ),
                soil_water_potential=collect(
                    forecast_output.soil_water_potential[point=site_index]
                ),
                soil_thermal_conductivity=collect(
                    forecast_output.soil_thermal_conductivity[point=site_index]
                ),
                soil_humidity=collect(
                    forecast_output.soil_humidity[point=site_index]
                ),
            )

            forecast_forcing = egg_nest_forcing(
                forecast_result,
                forecast_day_range,
                nest_node,
                environment_pars,
            )

            historical_egg_result = historical_egg_results[site_index]

            if egg_caches[site_index] === nothing
                egg_caches[site_index] = init_egg_cache(
                    egg_model,
                    pars,
                    historical_egg_result.final_state,
                    soil_hydraulics,
                    forecast_forcing,
                    forecast_tspan;
                    save_trajectory,
                )
            end

            forecast_outcomes[site_index, member_index] = simulate_egg!(
                egg_caches[site_index],
                historical_egg_result.final_state,
                soil_hydraulics,
                forecast_forcing,
                forecast_tspan,
            )
        end
    end
end

# ── summarise and plot each site ─────────────────────────────────────────────

using Plots

for site_index in eachindex(points)
    site = points[site_index]
    site_name = site_names[site_index]
    site_slug = site_slugs[site_index]
    historical_egg_result = historical_egg_results[site_index]

    site_forecast_outcomes = [
        forecast_outcomes[site_index, member_index]
        for member_index in eachindex(members)
        if forecast_outcomes[site_index, member_index] !== nothing
    ]

    hatch_dates = Date[]

    if historical_egg_result.hatched
        label = "Hatched during the historical period, before $issue_date."
    elseif historical_egg_result.died
        label = "Died during the historical period:\n" *
                "$(historical_egg_result.death_cause)."
    else
        println("\nPer-member outcomes for $site_name:")

        for (member, result) in zip(members, site_forecast_outcomes)
            if result.hatched
                hatch_date = issue_date +
                             Day(round(Int, ustrip(u"d", result.hatch_time)))
                push!(hatch_dates, hatch_date)
                println("  member $member -> hatched on $hatch_date")
            elseif result.died
                println("  member $member -> died of $(result.death_cause)")
            else
                println("  member $member -> did not hatch within " *
                        "$forecast_horizon_days forecast days")
            end
        end

        label = if isempty(hatch_dates)
            "0/$(length(members)) members hatched.\n" *
            "Median hatch: none\n" *
            "Range: none"
        else
            hatch_days_since_issue = [
                Dates.value(date - issue_date) for date in hatch_dates
            ]
            median_hatch_date = issue_date +
                                Day(round(Int, median(hatch_days_since_issue)))

            "$(length(hatch_dates))/$(length(members)) members hatched.\n" *
            "Median hatch: $median_hatch_date\n" *
            "Range: $(minimum(hatch_dates)) to $(maximum(hatch_dates))"
        end
    end

    println("$site_name: " * replace(label, '\n' => ' '))

    p1 = plot(;
        ylabel="development",
        title="$site_name, lay date: $oviposition_date, diapause = $diapause",
        legend=false,
        ylims=(0, 1),
    )
    p2 = plot(;
        ylabel="egg mass",
        legend=false,
        ylims=(0.0u"mg", pars.initial_egg_mass * 2.0),
    )
    p3 = plot(;
        ylabel="egg\ntemperature",
        xlabel="date",
        legend=false,
    )

    silo_t = ustrip.(u"s", historical_egg_result.trajectory.t)

    if isempty(site_forecast_outcomes)
        actual_times = DateTime(oviposition_date) .+
                       Dates.Second.(round.(Int, silo_t))

        plot!(
            p1,
            actual_times,
            historical_egg_result.trajectory.development_fraction,
        )
        plot!(
            p2,
            actual_times,
            collect(uconvert.(u"mg", historical_egg_result.trajectory.egg_mass)),
        )
        plot!(
            p3,
            actual_times,
            collect(uconvert.(u"°C", historical_egg_result.trajectory.temperature)),
        )
    else
        for result in site_forecast_outcomes
            traj = result.trajectory
            access_t = ustrip.(u"s", traj.t) .+ last(silo_t)
            combined_t = vcat(silo_t, access_t)

            development_fraction = vcat(
                historical_egg_result.trajectory.development_fraction,
                traj.development_fraction,
            )
            egg_mass = vcat(
                historical_egg_result.trajectory.egg_mass,
                traj.egg_mass,
            )
            temperature = vcat(
                historical_egg_result.trajectory.temperature,
                traj.temperature,
            )

            actual_times = DateTime(oviposition_date) .+
                           Dates.Second.(round.(Int, combined_t))

            plot!(p1, actual_times, development_fraction)
            plot!(p2, actual_times, collect(uconvert.(u"mg", egg_mass)))
            plot!(p3, actual_times, collect(uconvert.(u"°C", temperature)))
        end
    end

    # Place text just inside the left plot boundary so it is not clipped.
    annotate!(
        p1,
        DateTime(oviposition_date) + Day(2),
        0.85,
        text(label, 10, :black, :left),
    )

    combined_plot = plot(
        p1, p2, p3;
        layout=(3, 1),
        size=(1000, 900),
        link=:x,
    )

    plot_path = joinpath(
        output_dir,
        "history_forecast_splice_$(site_slug).png",
    )
    savefig(combined_plot, plot_path)
    println("Saved $plot_path")
    display(combined_plot)
end
