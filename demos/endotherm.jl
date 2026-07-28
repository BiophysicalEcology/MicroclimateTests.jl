# Endotherm thermoregulation driven by MicroclimateMapper microclimates and
# BiophysicalBehaviour's endotherm thermoregulation model.
#
# Swap WEATHER_SOURCE/DEM_SOURCE below (e.g. CRUCL2, TerraClimate, SILO, BARRA)
# to run the same organism at the same point against a different climate
# dataset — MicroclimateMapper resolves any of them onto the same interface.

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model
using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using FluidProperties
using Rasters, RasterDataSources
using Dates, Unitful, UnitfulMoles, Statistics
using Plots

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

# ── Data source (swap here) ────────────────────────────────────────────────
const WEATHER_SOURCE = CRUCL2     # or TerraClimate, SILO, BARRA, ...
const DEM_SOURCE      = CRUCL2
const YEAR             = 2000      # ignored by CRUCL2's climatology; used by TerraClimate etc.

location = "Palm Springs, CA"
points = [geocode(location)]
site_name = points[1].display_name
dates  = Date(YEAR, 1, 1):Day(1):Date(YEAR, 12, 31)

months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
hours  = collect(0.0:1.0:23.0)

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.1, 2.0]u"m"    # ground node + 2 m reference height

# ── Step 1: Microclimate ───────────────────────────────────────────────────
# Extra layers (beyond air_temperature/relative_humidity/wind_speed used by the
# steady-state loop below) are needed later for the transient model's `EnvironmentForcing`.
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:air_temperature, :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:diffuse_fraction, :scalar),
    LayerSpec(:pressure, :scalar),
    LayerSpec(:zenith_angle, :solar),
)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
    ),
    dem_source              = DEM_SOURCE,
    weather_source          = WEATHER_SOURCE,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
    output_layers,
)

println("Solving microclimate...")
cache  = init(MicroVectorProblem(;
    model, points, dates,
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths))),
))
output = solve!(cache)

elevation = terrain(cache).elevation[point=1]

T_air_series = collect(output.air_temperature[point=1, height=1])    # 2 m reference height
T_sky_series   = collect(output.sky_temperature[point=1])
T_ground_series = collect(output.soil_temperature[point=1, depth=1])
rh_series    = collect(output.relative_humidity[point=1, height=1])
wind_series  = collect(output.wind_speed[point=1, height=1])
solar_radiation_series = collect(output.global_radiation[point=1])
zenith_series = collect(output.zenith_angle[point=1])
soil_k_series  = collect(output.soil_thermal_conductivity[point=1, depth=1])
diffuse_series = collect(output.diffuse_fraction[point=1])

# Step count comes from the solved output, not `length(dates)`: monthly-
# climatology sources (e.g. CRUCL2) solve one representative day per month
# (12 × 24 = 288 steps) regardless of the requested date range.
nsteps = length(T_air_series)
ndays  = nsteps ÷ 24

# ── Step 2: Set up endotherm ──────────────────────────────────────────────
shape_pars       = example_shape_pars(mass = 1.0u"kg", axis_ratio_b = 3.0, axis_ratio_c = 3.0)
insulation_pars  = example_insulation_pars(;
                    insulation_depth_dorsal  = 2.0u"mm",
                    insulation_depth_ventral = 2.0u"mm",
                    )
radiation_pars   = example_radiation_pars()
metabolic_heat_flow = metabolic_rate(McKechnieWolf(), shape_pars.mass)
metabolism_pars = example_metabolism_pars(; core_temperature = (38.0 + 273.15)u"K", q10 = 2, metabolic_heat_flow)
evaporation_pars = example_evaporation_pars()
respiration_pars = example_respiration_pars()

conduction_pars_internal = example_conduction_pars_internal()
fat = FatLayer(conduction_pars_internal.fat_fraction, conduction_pars_internal.fat_density)

mean_insulation_depth = insulation_pars.dorsal.depth * (1 - radiation_pars.ventral_fraction) +
                        insulation_pars.ventral.depth * radiation_pars.ventral_fraction
mean_fibre_diameter   = insulation_pars.dorsal.diameter * (1 - radiation_pars.ventral_fraction) +
                        insulation_pars.ventral.diameter * radiation_pars.ventral_fraction
mean_fibre_density    = insulation_pars.dorsal.density * (1 - radiation_pars.ventral_fraction) +
                        insulation_pars.ventral.density * radiation_pars.ventral_fraction
fur      = FibrousLayer(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

physiology_traits = HeatExchangeTraits(
    shape_pars,
    insulation_pars,
    example_conduction_pars_external(),
    conduction_pars_internal,
    radiation_pars,
    ConvectionParameters(),
    evaporation_pars,
    example_hydraulic_pars(),
    respiration_pars,
    metabolism_pars,
    example_metabolic_rate_options(),
)

core_temperature_ref = metabolism_pars.core_temperature

thermoregulation_limits = ThermoregulationLimits(;
    control          = RuleBasedSequentialControl(; mode = CorePantingSweatingFirst(), tolerance = 0.005, max_iterations = 200),
    minimum_heat_flow = metabolism_pars.metabolic_heat_flow,
    insulation    = InsulationLimits(;
        dorsal  = SteppedParameter(; current   = insulation_pars.dorsal.depth,
                                     reference = insulation_pars.dorsal.depth,
                                     max       = insulation_pars.dorsal.depth, step = 0.0),
        ventral = SteppedParameter(; current   = insulation_pars.ventral.depth,
                                     reference = insulation_pars.ventral.depth,
                                     max       = insulation_pars.ventral.depth, step = 0.0),
    ),
    axis_ratio_factor = SteppedParameter(; current = shape_pars.axis_ratio_b, max = 5.0, step = 0.1),
    flesh_conductivity   = SteppedParameter(; current = conduction_pars_internal.flesh_conductivity,
                                       max = 2.8u"W/m/K", step = 0.1u"W/m/K"),
    core_temperature     = SteppedParameter(; current = core_temperature_ref, reference = core_temperature_ref,
                                       max = core_temperature_ref + 5.0u"K", step = 0.1u"K"),
    panting = PantingLimits(;
        pant       = SteppedParameter(; current = respiration_pars.pant, max = 15.0, step = 0.5),
        cost       = 0.0u"W",
        multiplier = 1.0,
        core_temperature_ref,
    ),
    skin_wetness = SteppedParameter(; current = evaporation_pars.skin_wetness,
                                       max = 1.0, step = 0.05),
)

behavioral_traits = BehavioralTraits(;
    thermoregulation = thermoregulation_limits,
    activity_period  = Diurnal(),
)
organism_traits = OrganismTraits(Endotherm(), physiology_traits, behavioral_traits)
organism        = Organism(geometry, organism_traits)

environment_pars = example_environment_pars(; elevation)

# ── Step 3: Thermoregulation loop ─────────────────────────────────────────
# Warm-start: carry skin/insulation temperature forward between hours for faster convergence.
skin_temperature       = core_temperature_ref - 3.0u"K"
insulation_temperature = T_air_series[1]
metabolic_heat_flow    = 0.0u"W"

endo_results = Vector{Any}(undef, nsteps)
println("Running endotherm thermoregulation loop...")

let skin_temperature = skin_temperature, insulation_temperature = insulation_temperature
    for step in 1:nsteps
        environment_vars = EnvironmentalVars(;
            air_temperature        = T_air_series[step],
            sky_temperature        = T_sky_series[step],
            ground_temperature     = T_ground_series[step],
            substrate_temperature  = T_ground_series[step],
            relative_humidity      = rh_series[step],
            wind_speed              = wind_series[step],
            atmospheric_pressure    = atmospheric_pressure(elevation),
            zenith_angle            = zenith_series[step],
            substrate_conductivity  = soil_k_series[step],
            global_radiation        = solar_radiation_series[step],
            diffuse_fraction        = diffuse_series[step],
            shade                   = 0.0,   # fully exposed
        )

        out = thermoregulate(
            organism,
            (; environment_pars, environment_vars),
            (; metabolic_heat_flow, skin_temperature, insulation_temperature),
        )

        endo_results[step]     = out
        skin_temperature       = out.thermoregulation.skin_temperature
        insulation_temperature = out.thermoregulation.insulation_temperature
    end
end

# ── Extract outputs ───────────────────────────────────────────────────────
T_air_C      = uconvert.(u"°C", T_air_series[1:nsteps])
T_core_C     = [uconvert(u"°C",   endo_results[i].thermoregulation.core_temperature) for i in 1:nsteps]
metabolic_heat_flow_W = [endo_results[i].energy_flows.metabolic_heat_flow for i in 1:nsteps]
m_evap_gh    = [endo_results[i].mass_flows.m_evap                 for i in 1:nsteps]
m_resp_gh    = [endo_results[i].mass_flows.respiration_mass_flow  for i in 1:nsteps]
m_sweat_gh   = [endo_results[i].mass_flows.m_sweat                for i in 1:nsteps]
axis_ratio_b = [endo_results[i].thermoregulation.axis_ratio_b                      for i in 1:nsteps]
skin_wetness = [endo_results[i].thermoregulation.skin_wetness                      for i in 1:nsteps]
pant         = [endo_results[i].thermoregulation.pant                              for i in 1:nsteps]

month_ranges    = [(m-1)*24+1 : m*24 for m in 1:ndays]
month_Ta        = [T_air_C[r]       for r in month_ranges]
month_Tc        = [T_core_C[r]      for r in month_ranges]
month_metabolic_heat_flow = [metabolic_heat_flow_W[r] for r in month_ranges]
month_evap      = [m_evap_gh[r]     for r in month_ranges]
month_resp      = [m_resp_gh[r]     for r in month_ranges]
month_sweat     = [m_sweat_gh[r]    for r in month_ranges]
month_axis_b    = [axis_ratio_b[r]  for r in month_ranges]
month_wetness   = [skin_wetness[r]  for r in month_ranges]
month_pant      = [pant[r]          for r in month_ranges]

println("\n── Annual metabolic summary ──")
println("  Mean metabolic_heat_flow: $(round(typeof(metabolic_heat_flow_W[1]), mean(metabolic_heat_flow_W); digits=2))")
println("  Max  metabolic_heat_flow: $(round(typeof(metabolic_heat_flow_W[1]), maximum(metabolic_heat_flow_W); digits=2))")
println("  Min  metabolic_heat_flow: $(round(typeof(metabolic_heat_flow_W[1]), minimum(metabolic_heat_flow_W); digits=2))")
println("\n── Annual water loss summary ──")
println("  Mean total evap: $(round(typeof(m_evap_gh[1]), mean(m_evap_gh); digits=3))")
println("  Mean resp loss:  $(round(typeof(m_evap_gh[1]), mean(m_resp_gh); digits=3))")
println("  Mean cutaneous:  $(round(typeof(m_evap_gh[1]), mean(m_sweat_gh); digits=3))")

# ── Fig. 1 – Core temperature by month (4×3 grid) ────────────────────────
panels_Tc = map(1:ndays) do m
    p = plot(hours, month_Tc[m];
        lw = 2, color = :red, label = "",
        title = months[mod1(m, 12)], titlefontsize = 9)
    plot!(p, hours, month_Ta[m];
        lw = 1, color = :steelblue, linestyle = :dash, label = "")
    p
end

display(plot(panels_Tc...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Core temperature, $site_name\n" *
                 "(red = T_core, blue dashed = T_air)"))

# ── Fig. 1b – Metabolic heat flow by month (4×3 grid) ─────────────────────
panels_mhf_ss = map(1:ndays) do m
    plot(hours, month_metabolic_heat_flow[m];
        lw = 2, color = :firebrick, label = "",
        title = months[mod1(m, 12)], titlefontsize = 9)
end
display(plot(panels_mhf_ss...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Metabolic heat flow, $site_name"))

# ── Fig. 1c – Water loss rates by month (4×3 grid) ─────────────────────────
panels_evap = map(1:ndays) do m
    p = plot(hours, month_evap[m];
        lw = 2, color = :teal, label = "total",
        title = months[mod1(m, 12)], titlefontsize = 9,
        legend = m == 1 ? :topright : false)
    plot!(p, hours, month_resp[m];  lw = 1, color = :orange, linestyle = :dash, label = "respiratory")
    plot!(p, hours, month_sweat[m]; lw = 1, color = :purple, linestyle = :dot,  label = "cutaneous")
    p
end
display(plot(panels_evap...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Evaporative water loss, $site_name"))

# ── Fig. 1d – Posture (axis ratio b) by month (4×3 grid) ───────────────────
panels_axis_b = map(1:ndays) do m
    plot(hours, month_axis_b[m];
        lw = 2, color = :goldenrod, label = "",
        title = months[mod1(m, 12)], ylabel = "axis ratio b", titlefontsize = 9)
end
display(plot(panels_axis_b...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Posture, $site_name"))

# ── Fig. 1e – Pant rate by month (4×3 grid) ────────────────────────────────
panels_pant = map(1:ndays) do m
    plot(hours, month_pant[m];
        lw = 2, color = :crimson, label = "",
        title = months[mod1(m, 12)], ylabel = "pant rate", titlefontsize = 9)
end
display(plot(panels_pant...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Pant rate, $site_name"))

# ── Fig. 2 – Annual heatmaps (core temperature and metabolic heat flow) ───
tc_matrix = zeros(Float64, 24, ndays)
mhf_matrix = zeros(Float64, 24, ndays)
for m in 1:ndays
    tc_matrix[:, m]  = ustrip.(u"°C", month_Tc[m])
    mhf_matrix[:, m] = ustrip.(u"W", month_metabolic_heat_flow[m])
end

p_tc = heatmap(1:ndays, hours, tc_matrix;
    color = cgrad(:RdYlBu, rev = true),
    colorbar_title = "°C",
    title = "Core temperature (°C)", ylabel = "hour", xlabel = "day of year")

p_mhf = heatmap(1:ndays, hours, mhf_matrix;
    color = :heat,
    colorbar_title = "W",
    title = "Metabolic heat flow (W)", ylabel = "hour", xlabel = "day of year")

display(plot(p_tc, p_mhf; layout = (2, 1), size = (900, 600), left_margin = 6Plots.mm))

# ── Fig. 3 – Annual heatmap of total water loss ───────────────────────────
evap_matrix = zeros(Float64, 24, ndays)
for m in 1:ndays
    evap_matrix[:, m] = ustrip.(u"g/hr", month_evap[m])
end

p_evap = heatmap(1:ndays, hours, evap_matrix;
    color = :Blues,
    colorbar_title = "g/hr",
    title = "Total evaporative water loss (g/hr)", ylabel = "hour", xlabel = "day of year")

display(plot(p_evap; size = (900, 350), left_margin = 6Plots.mm))

# ── Fig. 4 – Annual heatmaps (posture, skin wetness, panting) ─────────────
sb_matrix  = zeros(Float64, 24, ndays)
sw_matrix  = zeros(Float64, 24, ndays)
pt_matrix  = zeros(Float64, 24, ndays)
for m in 1:ndays
    sb_matrix[:, m] = month_axis_b[m]
    sw_matrix[:, m] = month_wetness[m]
    pt_matrix[:, m] = month_pant[m]
end

p_sb = heatmap(1:ndays, hours, sb_matrix;
    color = :YlOrBr, colorbar_title = "",
    title = "Posture (axis ratio b)", ylabel = "hour", xlabel = "day of year")
p_sw = heatmap(1:ndays, hours, sw_matrix;
    color = :Blues, colorbar_title = "",
    title = "Skin wetness", ylabel = "hour", xlabel = "day of year")
p_pt = heatmap(1:ndays, hours, pt_matrix;
    color = :Reds, colorbar_title = "",
    title = "Pant rate", ylabel = "hour", xlabel = "day of year")

display(plot(p_sb, p_sw, p_pt; layout = (3, 1), size = (900, 800), left_margin = 6Plots.mm))

# =============================================================================
# Transient (thermal-mass) endotherm model — active/resting cycling
#
# The steady-state loop above assumes body temperature reaches equilibrium
# instantly every hour. Real animals have thermal mass, so it takes time to
# heat up or cool down. `simulate_endotherm_activity_cycle` instead tracks
# core temperature continuously (an ODE): the animal is "active" (elevated
# metabolic rate + wind, e.g. foraging) until it gets too hot, then rests
# (baseline metabolic rate) until it cools back down, and repeats.
# =============================================================================

pressure_series = collect(output.pressure[point=1])

# Resting uses ground-level (sheltered) conditions; active uses the 2 m
# reference height, with wind speed bumped up to at least the animal's
# movement speed — same "height choice" idea as
# BiophysicalBehaviour.jl's examples/endotherm_transient.jl.
T_air_local_series = collect(output.air_temperature[point=1, height=1])
rh_local_series     = collect(output.relative_humidity[point=1, height=1])
wind_local_series   = collect(output.wind_speed[point=1, height=1])
movement_speed      = 3.0u"m/s"   # e.g. walking/foraging pace

# Builds one day's `EnvironmentForcing` from column slices - `shade` fixes how
# exposed the animal is (0 = full sun, close to 1 = sheltered).
# Parameters are named to match `EnvironmentalVarsVec`'s fields, so the constructor
# call below can pun them directly.
function day_forcing(times_day, air_temperature, relative_humidity, wind_speed, ground_temperature,
                      substrate_conductivity, sky_temperature, global_radiation, diffuse_fraction,
                      atmospheric_pressure, zenith_angle, shade)
    n = length(times_day)
    EnvironmentForcing(times_day, HeatExchange.EnvironmentalVarsVec(;
        air_temperature, sky_temperature,
        ground_temperature, substrate_temperature = ground_temperature,
        relative_humidity, wind_speed,
        atmospheric_pressure, zenith_angle,
        substrate_conductivity, global_radiation,
        diffuse_fraction, shade = fill(shade, n),
    ))
end

resting_metabolic_heat_flow = metabolism_pars.metabolic_heat_flow
active_metabolic_heat_flow  = resting_metabolic_heat_flow * 3.0   # e.g. sustained foraging
active_temperature_max      = core_temperature_ref + 2.0u"K"      # hyperthermia tolerance before resting
resume_temperature          = core_temperature_ref                # cools back to baseline before going active again

println("Running transient (active/resting cycle) endotherm model for each representative day...")
transient_results = Vector{Any}(undef, ndays)

@time for m in 1:ndays
    day_range = (m-1)*24+1 : m*24
    times     = (0:23) .* 1.0u"hr" .|> u"s"

    resting_forcing = day_forcing(times,
        T_air_local_series[day_range], rh_local_series[day_range], wind_local_series[day_range],
        T_ground_series[day_range], soil_k_series[day_range], T_sky_series[day_range],
        solar_radiation_series[day_range], diffuse_series[day_range], pressure_series[day_range],
        zenith_series[day_range], 0.9)
    active_forcing = day_forcing(times,
        T_air_series[day_range], rh_series[day_range],
        max.(ustrip.(u"m/s", movement_speed), ustrip.(u"m/s", wind_series[day_range])) .* u"m/s",
        T_ground_series[day_range], soil_k_series[day_range], T_sky_series[day_range],
        solar_radiation_series[day_range], diffuse_series[day_range], pressure_series[day_range],
        zenith_series[day_range], 0.0)

    transient_results[m] = simulate_endotherm_activity_cycle(
        times, core_temperature_ref, organism, environment_pars, active_forcing, resting_forcing;
        active_metabolic_heat_flow, resting_metabolic_heat_flow, active_temperature_max, resume_temperature,
    )
end

# Fraction of simulated time spent active, weighting by how long each accepted
# solver step actually lasted (a plain step count would over-weight bouts
# where the solver takes many small steps).
function active_time_fraction(t, state)
    total = 0.0u"s"
    for i in eachindex(t)
        state[i] isa Active || continue
        seg_end = i < lastindex(t) ? t[i + 1] : t[end]
        total += seg_end - t[i]
    end
    return total / (t[end] - t[1])
end

active_fractions = [active_time_fraction(r.t, r.state) for r in transient_results]
transient_min = minimum(minimum(ustrip.(u"°C", r.core_temperature)) for r in transient_results)
transient_max = maximum(maximum(ustrip.(u"°C", r.core_temperature)) for r in transient_results)

println("\n── Transient active/resting summary ──")
println("  Mean fraction of time active: $(round(100 * mean(active_fractions); digits=1))%")
println("  Core temperature range (transient): $(round(transient_min; digits=1)) to $(round(transient_max; digits=1)) °C")
println("  Core temperature range (steady-state): $(round(typeof(T_core_C[1]), minimum(T_core_C); digits=1)) to $(round(typeof(T_core_C[1]), maximum(T_core_C); digits=1)) °C")

# ── Fig. 5 – Core temperature by month: steady-state vs transient overlay ──
panels_transient = map(1:ndays) do m
    p = plot(hours, month_Tc[m];
        lw = 2, color = :red, label = "steady-state",
        title = months[mod1(m, 12)], ylabel = "°C", titlefontsize = 9,
        legend = m == 1 ? :topright : false)
    plot!(p, ustrip.(u"hr", transient_results[m].t), ustrip.(u"°C", transient_results[m].core_temperature);
        lw = 2, color = :darkgreen, label = "transient (active/resting)")
    p
end

display(plot(panels_transient...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Core temperature — steady-state vs transient, $site_name\n" *
                 "(red = steady-state, green = transient thermal-mass model)"))

# ── Fig. 5b – Metabolic heat flow by month: steady-state vs transient ─────
# Transient metabolic rate is a known input (active vs resting), not solved -
# plotted as a step function switching on `state`.
panels_mhf_trans = map(1:ndays) do m
    p = plot(hours, month_metabolic_heat_flow[m];
        lw = 2, color = :red, label = "steady-state",
        title = months[mod1(m, 12)], ylabel = "W", titlefontsize = 9,
        legend = m == 1 ? :topright : false)
    t_h = ustrip.(u"hr", transient_results[m].t)
    mhf_trans = [s isa Active ? ustrip(u"W", active_metabolic_heat_flow) : ustrip(u"W", resting_metabolic_heat_flow)
                 for s in transient_results[m].state]
    plot!(p, t_h, mhf_trans; lw = 2, color = :darkgreen, label = "transient", seriestype = :steppost)
    p
end
display(plot(panels_mhf_trans...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Metabolic heat flow — steady-state vs transient, $site_name"))

# ── Fig. 5c/5d – Posture and pant rate: steady-state vs transient. Both
#    effectors are held fixed at the organism's baseline value throughout
#    the transient model (only core temperature is solved) - shown as a flat
#    reference line rather than a solved trajectory.
panels_axis_b_trans = map(1:ndays) do m
    p = plot(hours, month_axis_b[m];
        lw = 2, color = :red, label = "steady-state",
        title = months[mod1(m, 12)], ylabel = "axis ratio b", titlefontsize = 9,
        legend = m == 1 ? :topright : false)
    hline!(p, [shape_pars.axis_ratio_b]; lw = 2, color = :darkgreen, label = "transient (fixed)")
    p
end
display(plot(panels_axis_b_trans...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Posture — steady-state vs transient, $site_name"))

panels_pant_trans = map(1:ndays) do m
    p = plot(hours, month_pant[m];
        lw = 2, color = :red, label = "steady-state",
        title = months[mod1(m, 12)], ylabel = "pant rate", titlefontsize = 9,
        legend = m == 1 ? :topright : false)
    hline!(p, [respiration_pars.pant]; lw = 2, color = :darkgreen, label = "transient (fixed)")
    p
end
display(plot(panels_pant_trans...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Pant rate — steady-state vs transient, $site_name"))

# =============================================================================
# Grid simulation — the same models applied across a raster instead of one
# point (parallels demos/ectotherm.jl's grid sections): a baseline with no
# thermoregulatory adjustment, the steady-state effector-adjustment loop, and
# the transient active/resting model, each run per pixel.
# =============================================================================

# ── Disk caching for the raster microclimate solve ──────────────────────────
# Set `USE_DISK_CACHE = true` to write the solved grid to disk and reload it
# on a later run instead of re-solving (see demos/ectotherm.jl for details).
const USE_DISK_CACHE = false
const CACHE_DIR = joinpath(@__DIR__, "cache")
isdir(CACHE_DIR) || mkpath(CACHE_DIR)

# Native unit for each `output_layers` entry above - Zarr can't carry Unitful
# types, so values are stripped before writing and reattached after reading.
const GRID_LAYER_UNITS = (
    soil_temperature          = u"K",
    soil_thermal_conductivity = u"W/m/K",
    air_temperature            = u"K",
    relative_humidity          = NoUnits,
    wind_speed                 = u"m/s",
    global_radiation           = u"W/m^2",
    sky_temperature            = u"K",
    diffuse_fraction           = NoUnits,
    pressure                   = u"Pa",
    zenith_angle               = u"°",
)

strip_grid_units(stack) = RasterStack(NamedTuple{propertynames(stack)}(
    map(n -> ustrip.(GRID_LAYER_UNITS[n], stack[n]), propertynames(stack))))
reattach_grid_units(stack) = RasterStack(NamedTuple{propertynames(stack)}(
    map(n -> stack[n] .* GRID_LAYER_UNITS[n], propertynames(stack))))

function solve_grid_cached(problem, cache_name)
    path = joinpath(CACHE_DIR, cache_name * ".zarr")
    if USE_DISK_CACHE && ispath(path)
        println("Loading cached grid microclimate: $cache_name")
        return reattach_grid_units(RasterStack(path))
    end
    println("Solving grid microclimate: $cache_name")
    cache  = init(problem)
    output = solve!(cache)
    cache  = nothing
    if USE_DISK_CACHE
        println("Caching grid microclimate to disk: $path")
        ispath(path) && rm(path; recursive = true)
        write(path, strip_grid_units(output))
    end
    GC.gc()
    return output
end

# ── Grid microclimate ────────────────────────────────────────────────────────
# One solve covers all three grid sections below: unlike demos/ectotherm.jl,
# nothing here needs a separate shaded/unshaded run - "sheltered" vs "exposed"
# is expressed via the `shade`/height choice already used at the point level.
grid_site     = geocode(location; buffer = 2.0)
raster_area   = grid_site.extent
grid_template = load_template(WEATHER_SOURCE, grid_site)

@time grid_output = solve_grid_cached(MicroRasterProblem(;
    model, area = raster_area, dates,
    template     = grid_template,
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
), "endotherm_grid")

grid_shape = grid_output.soil_temperature[depth=1]   # (X, Y, Ti) template - shape/dims only
out_dims   = dims(grid_shape)
nx, ny, grid_nsteps = size(grid_shape)
grid_ndays = grid_nsteps ÷ 24

_nan_layer(u) = Raster(fill(NaN * u, nx, ny, grid_nsteps), out_dims)

# No-data pixels (ocean, edge-of-coverage cells) come back as `missing` or `NaN`.
_no_data(v) = ismissing(v) || isnan(ustrip(v))

# Thread-safe progress reporting, logged periodically with a simple ETA.
function grid_progress_reporter(label, total; report_every = max(1, total ÷ 20))
    counter = Threads.Atomic{Int}(0)
    t_start = time()
    return function ()
        n = Threads.atomic_add!(counter, 1) + 1
        if n % report_every == 0 || n == total
            elapsed = time() - t_start
            eta = elapsed / n * (total - n)
            @info "$label: $n / $total pixels ($(round(Int, 100n / total))%) — ETA $(round(Int, eta))s"
        end
        return nothing
    end
end

# ── Grid section 1: baseline (no thermoregulatory adjustment) ──────────────
# `solve_metabolic_rate` fixes core temperature and the organism's default
# (as-built) effectors, solving only for the metabolic rate needed to balance
# the heat budget - the endotherm parallel of ectotherm.jl's "operative
# temperature" grid: what physiological cost the environment imposes with no
# compensating response.
grid_baseline_result = RasterStack((;
    metabolic_heat_flow    = _nan_layer(u"W"),
    evaporative_water_loss = _nan_layer(u"g/hr"),
    skin_temperature        = _nan_layer(u"K"),
))

function solve_baseline_grid!(grid_baseline_result, grid_output, organism, environment_pars,
                               core_temperature_ref, elevation, nx::Int, ny::Int, grid_nsteps::Int)
    report_progress = grid_progress_reporter("Baseline grid", nx * ny)
    Threads.@threads :static for xy in CartesianIndices((nx, ny))
        x, y = Tuple(xy)

        T_air_col = collect(grid_output.air_temperature[X(x), Y(y), height=2])
        if !_no_data(T_air_col[1])
            rh_col   = collect(grid_output.relative_humidity[X(x), Y(y), height=2])
            wind_col = collect(grid_output.wind_speed[X(x), Y(y), height=2])

            skin_temperature       = core_temperature_ref - 3.0u"K"
            insulation_temperature = u"K"(10.0u"°C")
            @inbounds for step in 1:grid_nsteps
                environment_vars = example_environment_vars(;
                    air_temperature      = T_air_col[step],
                    relative_humidity    = rh_col[step],
                    wind_speed           = wind_col[step],
                    atmospheric_pressure = atmospheric_pressure(elevation),
                    global_radiation     = 0.0u"W/m^2",   # sheltered, as at the point level above
                    zenith_angle         = 20.0u"°",
                )
                # solve_metabolic_rate takes skin/insulation temperature as separate positional
                # arguments (not a bundled `init` NamedTuple like `thermoregulate` below).
                out = solve_metabolic_rate(
                    organism, (; environment_pars, environment_vars), skin_temperature, insulation_temperature,
                )

                grid_baseline_result.metabolic_heat_flow[x, y, step]    = out.energy_flows.metabolic_heat_flow
                grid_baseline_result.evaporative_water_loss[x, y, step] = out.mass_flows.m_evap
                grid_baseline_result.skin_temperature[x, y, step]       = out.thermoregulation.skin_temperature

                skin_temperature       = out.thermoregulation.skin_temperature
                insulation_temperature = out.thermoregulation.insulation_temperature
            end
        end

        report_progress()
    end
    return grid_baseline_result
end

println("Solving baseline (no thermoregulation) grid...")
@time solve_baseline_grid!(grid_baseline_result, grid_output, organism, environment_pars,
    core_temperature_ref, elevation, nx, ny, grid_nsteps)

# ── Grid section 2: WITH steady-state thermoregulation ──────────────────────
# The same `thermoregulate` effector-adjustment loop used at the point level
# (Step 3 above), applied per pixel.
grid_thermoreg_result = RasterStack((;
    core_temperature       = _nan_layer(u"K"),
    metabolic_heat_flow    = _nan_layer(u"W"),
    evaporative_water_loss = _nan_layer(u"g/hr"),
    axis_ratio_b            = _nan_layer(NoUnits),
    skin_wetness             = _nan_layer(NoUnits),
    pant                     = _nan_layer(NoUnits),
))

function solve_thermoreg_grid!(grid_thermoreg_result, grid_output, organism, environment_pars,
                                core_temperature_ref, elevation, nx::Int, ny::Int, grid_nsteps::Int)
    report_progress = grid_progress_reporter("Thermoregulation grid", nx * ny)
    Threads.@threads :static for xy in CartesianIndices((nx, ny))
        x, y = Tuple(xy)

        T_air_col = collect(grid_output.air_temperature[X(x), Y(y), height=2])
        if !_no_data(T_air_col[1])
            rh_col   = collect(grid_output.relative_humidity[X(x), Y(y), height=2])
            wind_col = collect(grid_output.wind_speed[X(x), Y(y), height=2])

            skin_temperature       = core_temperature_ref - 3.0u"K"
            insulation_temperature = u"K"(10.0u"°C")
            metabolic_heat_flow    = 0.0u"W"
            @inbounds for step in 1:grid_nsteps
                environment_vars = example_environment_vars(;
                    air_temperature      = T_air_col[step],
                    relative_humidity    = rh_col[step],
                    wind_speed           = wind_col[step],
                    atmospheric_pressure = atmospheric_pressure(elevation),
                    global_radiation     = 0.0u"W/m^2",
                    zenith_angle         = 20.0u"°",
                )
                out = thermoregulate(organism, (; environment_pars, environment_vars),
                    (; metabolic_heat_flow, skin_temperature, insulation_temperature))

                grid_thermoreg_result.core_temperature[x, y, step]       = out.thermoregulation.core_temperature
                grid_thermoreg_result.metabolic_heat_flow[x, y, step]    = out.energy_flows.metabolic_heat_flow
                grid_thermoreg_result.evaporative_water_loss[x, y, step] = out.mass_flows.m_evap
                grid_thermoreg_result.axis_ratio_b[x, y, step]           = out.thermoregulation.axis_ratio_b
                grid_thermoreg_result.skin_wetness[x, y, step]           = out.thermoregulation.skin_wetness
                grid_thermoreg_result.pant[x, y, step]                   = out.thermoregulation.pant

                skin_temperature       = out.thermoregulation.skin_temperature
                insulation_temperature = out.thermoregulation.insulation_temperature
            end
        end

        report_progress()
    end
    return grid_thermoreg_result
end

println("Solving grid thermoregulation...")
@time solve_thermoreg_grid!(grid_thermoreg_result, grid_output, organism, environment_pars,
    core_temperature_ref, elevation, nx, ny, grid_nsteps)

# ── Grid section 3: WITH transient (active/resting) thermoregulation ────────
# The same `simulate_endotherm_activity_cycle` model used in the point-level
# transient section above, applied per pixel per representative day. Most
# expensive step in this file - each pixel-day runs a full adaptive-step ODE
# solve. Only per-day summaries are kept (not the full sub-daily trajectory).
day_out_dims = (out_dims[1], out_dims[2], Ti(1:grid_ndays))
_nan_day_layer(u) = Raster(fill(NaN * u, nx, ny, grid_ndays), day_out_dims)

grid_transient_result = RasterStack((;
    core_temperature_max = _nan_day_layer(u"K"),
    core_temperature_min = _nan_day_layer(u"K"),
    active_fraction        = _nan_day_layer(NoUnits),
))

function solve_transient_grid!(grid_transient_result, grid_output, organism, environment_pars,
                                core_temperature_ref, active_metabolic_heat_flow, resting_metabolic_heat_flow,
                                active_temperature_max, resume_temperature, movement_speed,
                                nx::Int, ny::Int, grid_ndays::Int)
    report_progress = grid_progress_reporter("Transient grid", nx * ny)
    Threads.@threads :static for xy in CartesianIndices((nx, ny))
        x, y = Tuple(xy)

        T_air_ref_col = collect(grid_output.air_temperature[X(x), Y(y), height=2])
        if !_no_data(T_air_ref_col[1])
            rh_ref_col      = collect(grid_output.relative_humidity[X(x), Y(y), height=2])
            wind_ref_col    = collect(grid_output.wind_speed[X(x), Y(y), height=2])
            T_air_local_col = collect(grid_output.air_temperature[X(x), Y(y), height=1])
            rh_local_col    = collect(grid_output.relative_humidity[X(x), Y(y), height=1])
            wind_local_col  = collect(grid_output.wind_speed[X(x), Y(y), height=1])
            soil_T_col      = collect(grid_output.soil_temperature[X(x), Y(y), depth=1])
            soil_k_col      = collect(grid_output.soil_thermal_conductivity[X(x), Y(y), depth=1])
            sky_T_col       = collect(grid_output.sky_temperature[X(x), Y(y)])
            solar_col       = collect(grid_output.global_radiation[X(x), Y(y)])
            diffuse_col     = collect(grid_output.diffuse_fraction[X(x), Y(y)])
            pressure_col    = collect(grid_output.pressure[X(x), Y(y)])
            zenith_col      = collect(grid_output.zenith_angle[X(x), Y(y)])

            for m in 1:grid_ndays
                day_range = (m-1)*24+1 : m*24
                times     = (0:23) .* 1.0u"hr" .|> u"s"

                resting_forcing = day_forcing(times,
                    T_air_local_col[day_range], rh_local_col[day_range], wind_local_col[day_range],
                    soil_T_col[day_range], soil_k_col[day_range], sky_T_col[day_range],
                    solar_col[day_range], diffuse_col[day_range], pressure_col[day_range],
                    zenith_col[day_range], 0.9)
                active_forcing = day_forcing(times,
                    T_air_ref_col[day_range], rh_ref_col[day_range],
                    max.(ustrip.(u"m/s", movement_speed), ustrip.(u"m/s", wind_ref_col[day_range])) .* u"m/s",
                    soil_T_col[day_range], soil_k_col[day_range], sky_T_col[day_range],
                    solar_col[day_range], diffuse_col[day_range], pressure_col[day_range],
                    zenith_col[day_range], 0.0)

                sol = simulate_endotherm_activity_cycle(
                    times, core_temperature_ref, organism, environment_pars, active_forcing, resting_forcing;
                    active_metabolic_heat_flow, resting_metabolic_heat_flow, active_temperature_max, resume_temperature,
                )

                grid_transient_result.core_temperature_max[x, y, m] = maximum(sol.core_temperature)
                grid_transient_result.core_temperature_min[x, y, m] = minimum(sol.core_temperature)
                grid_transient_result.active_fraction[x, y, m]      = ustrip(active_time_fraction(sol.t, sol.state))
            end
        end

        report_progress()
    end
    return grid_transient_result
end

println("Solving grid transient thermoregulation...")
@time solve_transient_grid!(grid_transient_result, grid_output, organism, environment_pars,
    core_temperature_ref, active_metabolic_heat_flow, resting_metabolic_heat_flow,
    active_temperature_max, resume_temperature, movement_speed, nx, ny, grid_ndays)

# ── Plots: baseline vs steady-state vs transient, across the raster ─────────
nanmax(A; dims) = (m = dropdims(maximum(ifelse.(isnan.(A), -Inf, A); dims); dims); ifelse.(isinf.(m), NaN, m))
nanmin(A; dims) = (m = dropdims(minimum(ifelse.(isnan.(A),  Inf, A); dims); dims); ifelse.(isinf.(m), NaN, m))

Tb_palette  = cgrad([:blue, :lightblue, :orange, :red, :purple])
MHF_palette = cgrad([:black, :orange, :red])

EWL_palette = cgrad([:white, :teal, :blue])

# ── Baseline: energy and water cost of the environment with no thermoregulatory response ──
mean_MHF_baseline = dropdims(mean(ustrip.(u"W", grid_baseline_result.metabolic_heat_flow); dims = Ti); dims = Ti) .* u"W"
max_MHF_baseline  = nanmax(ustrip.(u"W", grid_baseline_result.metabolic_heat_flow); dims = Ti) .* u"W"
mean_EWL_baseline = dropdims(mean(ustrip.(u"g/hr", grid_baseline_result.evaporative_water_loss); dims = Ti); dims = Ti) .* u"g/hr"
max_EWL_baseline  = nanmax(ustrip.(u"g/hr", grid_baseline_result.evaporative_water_loss); dims = Ti) .* u"g/hr"

display(plot(
    plot(mean_MHF_baseline; title = "Mean metabolic heat flow", seriestype = :heatmap, color = MHF_palette, yflip = false),
    plot(max_MHF_baseline;  title = "Max metabolic heat flow",  seriestype = :heatmap, color = MHF_palette, yflip = false),
    plot(mean_EWL_baseline; title = "Mean water loss",          seriestype = :heatmap, color = EWL_palette, yflip = false),
    plot(max_EWL_baseline;  title = "Max water loss",           seriestype = :heatmap, color = EWL_palette, yflip = false);
    layout = (2, 2), size = (1000, 800), left_margin = 6Plots.mm,
    plot_title = "Energy and water cost with no thermoregulation"))

# ── Steady-state: energy/water cost and thermoregulatory effector use ─────
mean_MHF_reg = dropdims(mean(ustrip.(u"W", grid_thermoreg_result.metabolic_heat_flow); dims = Ti); dims = Ti) .* u"W"
max_MHF_reg  = nanmax(ustrip.(u"W", grid_thermoreg_result.metabolic_heat_flow); dims = Ti) .* u"W"
mean_EWL_reg = dropdims(mean(ustrip.(u"g/hr", grid_thermoreg_result.evaporative_water_loss); dims = Ti); dims = Ti) .* u"g/hr"
max_EWL_reg  = nanmax(ustrip.(u"g/hr", grid_thermoreg_result.evaporative_water_loss); dims = Ti) .* u"g/hr"

display(plot(
    plot(mean_MHF_reg; title = "Mean metabolic heat flow", seriestype = :heatmap, color = MHF_palette, yflip = false),
    plot(max_MHF_reg;  title = "Max metabolic heat flow",  seriestype = :heatmap, color = MHF_palette, yflip = false),
    plot(mean_EWL_reg; title = "Mean water loss",          seriestype = :heatmap, color = EWL_palette, yflip = false),
    plot(max_EWL_reg;  title = "Max water loss",           seriestype = :heatmap, color = EWL_palette, yflip = false);
    layout = (2, 2), size = (1000, 800), left_margin = 6Plots.mm,
    plot_title = "Energy and water cost — steady-state thermoregulation"))

mean_axis_b_reg = dropdims(mean(grid_thermoreg_result.axis_ratio_b; dims = Ti); dims = Ti)
mean_wetness_reg = dropdims(mean(grid_thermoreg_result.skin_wetness; dims = Ti); dims = Ti)
mean_pant_reg = dropdims(mean(grid_thermoreg_result.pant; dims = Ti); dims = Ti)

display(plot(
    plot(mean_axis_b_reg;  title = "Mean posture (axis ratio b)", seriestype = :heatmap, color = :YlOrBr, yflip = false),
    plot(mean_wetness_reg; title = "Mean skin wetness",            seriestype = :heatmap, color = :Blues,  yflip = false),
    plot(mean_pant_reg;    title = "Mean pant rate",                seriestype = :heatmap, color = :Reds,   yflip = false);
    layout = (3, 1), size = (900, 900), left_margin = 6Plots.mm,
    plot_title = "Thermoregulatory effector use — steady-state (annual mean)"))

max_Tc_reg = nanmax(ustrip.(u"°C", grid_thermoreg_result.core_temperature); dims = Ti) .* u"°C"
min_Tc_reg = nanmin(ustrip.(u"°C", grid_thermoreg_result.core_temperature); dims = Ti) .* u"°C"
max_Tc_transient = nanmax(ustrip.(u"°C", grid_transient_result.core_temperature_max); dims = Ti) .* u"°C"
min_Tc_transient = nanmin(ustrip.(u"°C", grid_transient_result.core_temperature_min); dims = Ti) .* u"°C"

display(plot(
    plot(max_Tc_reg;        title = "Annual max Tc — steady-state", seriestype = :heatmap, color = Tb_palette, yflip = false),
    plot(max_Tc_transient;  title = "Annual max Tc — transient",    seriestype = :heatmap, color = Tb_palette, yflip = false),
    plot(min_Tc_reg;        title = "Annual min Tc — steady-state", seriestype = :heatmap, color = Tb_palette, yflip = false),
    plot(min_Tc_transient;  title = "Annual min Tc — transient",    seriestype = :heatmap, color = Tb_palette, yflip = false);
    layout = (2, 2), size = (1000, 800), left_margin = 6Plots.mm,
    plot_title = "Steady-state vs transient core temperature extremes"))

# Fraction of the year spent active (transient model), per pixel.
active_fraction_grid = dropdims(mean(grid_transient_result.active_fraction; dims = Ti); dims = Ti) .* 100
display(plot(active_fraction_grid;
    title = "Time spent active (% of simulated hours, transient model)",
    seriestype = :heatmap, color = cgrad([:steelblue, :orange, :firebrick]), yflip = false))

# Transient metabolic rate isn't tracked per step (it's a fixed input per phase,
# not solved) - implied mean cost = time-weighted average of the two fixed rates.
mean_active_fraction = dropdims(mean(grid_transient_result.active_fraction; dims = Ti); dims = Ti)
mean_MHF_transient = mean_active_fraction .* ustrip(u"W", active_metabolic_heat_flow) .+
    (1 .- mean_active_fraction) .* ustrip(u"W", resting_metabolic_heat_flow)
display(plot(mean_MHF_transient .* u"W";
    title = "Implied mean metabolic heat flow (transient model)",
    seriestype = :heatmap, color = MHF_palette, yflip = false))
