# Ectotherm thermoregulation driven by MicroclimateMapper microclimates and
# BiophysicalBehaviour's ectotherm thermoregulation model.
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
using Rasters, RasterDataSources
using ZarrDatasets   # activates Rasters.jl's Zarr read/write support (used for disk-caching grid solves)
using Dates, Unitful, Statistics
using Plots

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

# ── Data source (swap here) ────────────────────────────────────────────────
const WEATHER_SOURCE = CRUCL2     # or TerraClimate, SILO, BARRA, ...
const DEM_SOURCE     = CRUCL2
const YEAR           = 2000      # ignored by CRUCL2's climatology; used by TerraClimate etc.

# place and time
location = "Palm Springs, CA"
points = [geocode(location)]
site_name = points[1].display_name
dates = Date(YEAR, 1, 1):Day(1):Date(YEAR, 12, 31)

# just for plotting
months     = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
hours      = collect(0.0:1.0:23.0)
hour_edges = collect(0.0:1.0:24.0)   # bin edges for hourly_state_fraction below

# Microclimate solves the full depth/height suite (needed for a physically correct soil/air
# profile) - the steady-state organism's retreat range is pinned to a single node below
# (depth_min=depth_max, height_min=height_max) matching the transient driver's single fixed
# burrow depth and climb height, so the two drivers use exactly the same retreat site. Foraging
# stays at the true surface/ground node (depth_foraging/height_foraging default), independent
# of the retreat range.
depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 0.5, 1.0]u"m"

# The one depth/height site each animal is allowed to use, as physical values (snapped to the
# nearest node in `depths`/`heights` above).
burrow_depth = 10.0u"cm"
climb_height = 1.0u"m"

minimum_shade = 0.0
maximum_shade = 0.9

# ── Step 1: Microclimate at 0% and 90% shade ──────────────────────────────
# Extra output layers (beyond MicroclimateMapper's defaults) needed by
# BiophysicalBehaviour's AvailableEnvironments/interpolate_environment.
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
    LayerSpec(:air_temperature, :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:diffuse_fraction, :scalar),
    LayerSpec(:reference_temperature, :scalar),
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

problem_kwargs = (;
    model, points, dates,
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths))),
)

println("Solving microclimate (0% shade)...")
low_shade_cache   = init(MicroVectorProblem(; problem_kwargs...))
low_shade_output  = solve!(low_shade_cache)
println("Solving microclimate (90% shade)...")
high_shade_output = solve(MicroVectorProblem(; problem_kwargs..., data = (; shade = maximum_shade)))

elevation = terrain(low_shade_cache).elevation[point=1]

# ── Adapter: MicroclimateMapper's point-mode RasterStack → the plain
#    Matrix/Vector shape BiophysicalBehaviour's AvailableEnvironments expects.
function point_environment(output)
    profile = (;
        air_temperature   = collect(output.air_temperature[point=1]),
        relative_humidity = collect(output.relative_humidity[point=1]),
        wind_speed        = collect(output.wind_speed[point=1]),
    )
    return (;
        pressure                  = collect(output.pressure[point=1]),
        reference_temperature     = collect(output.reference_temperature[point=1]),
        global_radiation          = collect(output.global_radiation[point=1]),
        diffuse_fraction          = collect(output.diffuse_fraction[point=1]),
        sky_temperature           = collect(output.sky_temperature[point=1]),
        soil_temperature          = collect(output.soil_temperature[point=1]),
        soil_humidity             = collect(output.soil_humidity[point=1]),
        soil_thermal_conductivity = collect(output.soil_thermal_conductivity[point=1]),
        profile,
        solar_radiation = (; zenith_angle = collect(output.zenith_angle[point=1])),
    )
end

low_shade_result  = point_environment(low_shade_output)
high_shade_result = point_environment(high_shade_output)

# Step count comes from the solved output, not `length(dates)`: monthly-
# climatology sources (e.g. CRUCL2) solve one representative day per month
# (12 × 24 = 288 steps) regardless of the requested date range.
nsteps = length(low_shade_result.global_radiation)
ndays  = nsteps ÷ 24

available_environments = AvailableEnvironments(
    low_shade_result, high_shade_result, minimum_shade, maximum_shade, depths, heights
)

# ── Step 2: Organism ──────────────────────────────────────────────────────
shape_pars = DesertIguana(0.07u"kg", 1000.0u"kg/m^3")
body       = Body(shape_pars, Naked())

organism_traits = example_ectotherm_organism_traits(
    activity_period         = Diurnal(),
    
    target_temperature      = u"K"(38.5u"°C"),
    active_temperature_min  = u"K"(38.0u"°C"),
    active_temperature_max  = u"K"(43.0u"°C"),
    basking_temperature_min = u"K"(37.8u"°C"),
    emerge_temperature_min  = u"K"(15.0u"°C"),
    escape_temperature_min = u"K"(3.0u"°C"),
    escape_temperature_max = u"K"(44.0u"°C"),
    
    can_climb               = false,
    can_retreat_underground = true,
    can_seek_shade          = true,
    can_solar_orient        = true,
    can_press_to_ground     = true,
    can_change_absorptivity = true,
    
    burrow_shade_mode       = AdaptiveBurrowShade(), # MinShadeOnly()
    depths                  = depths,
    heights                 = heights,
    
    depth_foraging          = depths[1],      # surface
    height_foraging         = heights[1],     # ground level
    depth_min               = burrow_depth,   # pins the retreat search to the transient driver's fixed burrow depth
    depth_max               = burrow_depth,
    height_min              = climb_height,   # pins the climb search to the transient driver's fixed climb height
    height_max              = climb_height,
    shade_min               = minimum_shade,
    shade_max               = maximum_shade,
    
    absorptivity_min        = 0.6,
    absorptivity_max        = 0.8,
    absorptivity_step       = 0.01,
    
    can_pant                = false,
    pant_max                = 1.0,
    
    heat_exchange = example_ectotherm_heat_exchange_traits(;
        shape_pars,
        conduction_pars_external = example_ectotherm_conduction_pars_external(
            conduction_fraction = 0.1),
        evaporation_pars = example_ectotherm_evaporation_pars(
            eye_fraction = 0.0003, skin_wetness = 0.001),
        radiation_pars = example_ectotherm_radiation_pars(
            body_absorptivity_dorsal  = 0.8,
            body_absorptivity_ventral = 0.8,
            solar_orientation         = Intermediate(),
            body_emissivity_dorsal    = 0.95,
            body_emissivity_ventral   = 0.95),
        respiration_pars = example_ectotherm_respiration_pars(mouth_fraction = 0.0),
    ),
)

organism = Organism(body, organism_traits)
limits   = thermoregulation(organism)
env_pars = example_environment_pars(;
    elevation,
    ground_albedo = 0.15,
)

# ── Step 3: Thermoregulation loop ─────────────────────────────────────────
println("Running thermoregulation loop...")
results            = NamedTuple[]
previous_depth     = limits.depth.reference
activity_commenced = false

let previous_depth = previous_depth, activity_commenced = activity_commenced
    for step in 1:nsteps
        if (step - 1) % 24 == 0
            activity_commenced = false
        end
        out = thermoregulate(
            organism, available_environments, limits, env_pars, step, previous_depth;
            activity_commenced,
        )
        previous_depth     = out.depth_node
        activity_commenced = activity_commenced || out.state isa Active || out.state isa Basking
        push!(results, out)
    end
end

# ── Extract outputs ───────────────────────────────────────────────────────
T_body   = [r.core_temperature for r in results]
state    = [r.state            for r in results]
height   = [r.height           for r in results]
depth_nd = [r.depth_node       for r in results]
shade    = [r.shade            for r in results]
T_air    = low_shade_result.profile.air_temperature[:, 1]
T_sub    = low_shade_result.soil_temperature[:, 1]

T_body_C = ustrip.(u"°C", T_body)
act      = [s isa Active ? 2 : s isa Basking ? 1 : 0 for s in state]

_ground_ht = ustrip(u"cm", heights[1])
pos_cm = [depth_nd[i] > 1 ?
    -ustrip(u"cm", depths[depth_nd[i]]) :
    (h = ustrip(u"cm", height[i]); h > _ground_ht ? h : 0.0)
    for i in 1:nsteps]

month_ranges = [(m-1)*24+1 : m*24 for m in 1:ndays]
month_Tb     = [T_body_C[r]              for r in month_ranges]
month_Ta     = [ustrip.(u"°C", T_air[r]) for r in month_ranges]
month_Ts     = [ustrip.(u"°C", T_sub[r]) for r in month_ranges]
month_pos    = [pos_cm[r]               for r in month_ranges]
month_shade  = [shade[r]                for r in month_ranges]
month_state  = [state[r]                for r in month_ranges]

T_active_min_C = ustrip(u"°C", limits.active_temperature_min)
T_active_max_C = ustrip(u"°C", limits.active_temperature_max)

println("\n── Annual activity summary ──")
println("  Resting=$(sum(act.==0)), Basking=$(sum(act.==1)), Active=$(sum(act.==2))")

# Fraction of each hour bin spent in state `S`, from a piecewise-constant
# (t, state) series — state[i] holds from t[i] until t[i+1] (the last segment
# extends to `t_end`). For the steady-state hourly loop (one value per hour,
# `t = hours`, `t_end = 24.0`) every bin aligns exactly with a sample, so this
# reduces to a plain 0-or-1 indicator; for the transient model's adaptive
# event series, a bout boundary can fall mid-hour, so the fraction is
# continuously graded. Shared by the steady-state activity heatmaps (Fig. 2)
# and their transient counterparts below.
function hourly_state_fraction(t::AbstractVector{<:Real}, state::AbstractVector, t_end::Real,
                                ::Type{S}, hour_edges::AbstractVector{<:Real}) where S
    nbins = length(hour_edges) - 1
    total = zeros(Float64, nbins)
    for i in eachindex(t)
        state[i] isa S || continue
        seg_start = t[i]
        seg_end   = i < lastindex(t) ? t[i + 1] : t_end
        for b in 1:nbins
            lo, hi = hour_edges[b], hour_edges[b + 1]
            overlap = min(seg_end, hi) - max(seg_start, lo)
            overlap > 0 && (total[b] += overlap)
        end
    end
    return total ./ diff(hour_edges)
end

# Linear interpolation of a continuous transient series (e.g. core
# temperature) onto fixed query hours, for building an hourly heatmap
# comparable to the steady-state one.
function interp_at_hours(t_hours, y, query_hours)
    map(query_hours) do h
        h <= t_hours[1]   && return y[1]
        h >= t_hours[end] && return y[end]
        j = clamp(searchsortedlast(t_hours, h), 1, length(t_hours) - 1)
        frac = (h - t_hours[j]) / (t_hours[j + 1] - t_hours[j])
        return y[j] + frac * (y[j + 1] - y[j])
    end
end

# Zero-order-hold lookup of a categorical series onto query hours (state only
# changes at bout boundaries, so holding the last value is exact).
function state_at_hours(t_hours, state, query_hours)
    idx = clamp.(searchsortedlast.(Ref(t_hours), query_hours), 1, length(state))
    return state[idx]
end

# ── Fig. 1 – Body temperature by month (4×3 grid) ────────────────────────
panels_Tb = map(1:ndays) do m
    p = plot(hours, month_Tb[m];
        lw = 2, color = :red, label = "",
        title = months[mod1(m, 12)], ylabel = "°C", ylim = (0, 50), titlefontsize = 9)
    plot!(p, hours, month_Ta[m]; lw = 1, color = :steelblue, linestyle = :dash, label = "")
    plot!(p, hours, month_Ts[m]; lw = 1, color = :grey, linestyle = :dash, label = "")
    hline!(p, [T_active_min_C, T_active_max_C];
        color = :orange, linestyle = :dash, lw = 1, label = "")
    p
end

display(plot(panels_Tb...; layout = (ceil(Int, ndays/3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Body temperature, $site_name\n" *
                 "(red = Tb, blue dashed = T_air, orange dashed = active range)"))

# ── Fig. 2 – Annual heatmaps (body temp, active %, basking %) ─────────────
tb_matrix  = zeros(Float64, 24, ndays)
for m in 1:ndays
    tb_matrix[:, m] = month_Tb[m]
end

p1 = heatmap(1:ndays, hours, tb_matrix;
    color = cgrad(:RdYlBu, rev = true), clims = (0, 50),
    colorbar_title = "°C",
    title = "Body temperature (°C)", ylabel = "hour", xlabel = "day of year")

active_matrix  = zeros(Float64, 24, ndays)
basking_matrix = zeros(Float64, 24, ndays)
for m in 1:ndays
    active_matrix[:,  m] = hourly_state_fraction(hours, month_state[m], 24.0, Active,  hour_edges) .* 100
    basking_matrix[:, m] = hourly_state_fraction(hours, month_state[m], 24.0, Basking, hour_edges) .* 100
end

p2a = heatmap(1:ndays, hours, active_matrix;
    color = :Reds, clims = (0, 100),
    colorbar_title = "%", title = "Active/foraging (%)", ylabel = "hour", xlabel = "day of year")

p2b = heatmap(1:ndays, hours, basking_matrix;
    color = :Oranges, clims = (0, 100),
    colorbar_title = "%", title = "Basking (%)", ylabel = "hour", xlabel = "day of year")

display(plot(p1, p2a, p2b; layout = (3, 1), size = (900, 900), left_margin = 6Plots.mm))

# ── Fig. 3 – Annual heatmaps (shade and position) ─────────────────────────
# Scaled to the actually-accessible burrow depth/climb height, not the full microclimate
# depth/height suite (most of which the animal can never reach - see burrow_depth/
# climb_height above).
depth_cm_max  = ustrip(u"cm", burrow_depth)
height_cm_max = ustrip(u"cm", climb_height)
pos_clims      = (-depth_cm_max, height_cm_max)
total_range    = depth_cm_max + height_cm_max
norm(v)        = (depth_cm_max + v) / total_range
pos_shallowest = norm(-depth_cm_max)
pos_surface    = norm(0.0)
pos_lowest_ht  = norm(height_cm_max)
pos_cmap = cgrad(
    [:saddlebrown, :saddlebrown, :limegreen, :limegreen, :skyblue, :steelblue],
    [0.0,
     (pos_shallowest + pos_surface) / 2,
     pos_surface - 0.001,
     pos_surface + 0.001,
     (pos_surface + pos_lowest_ht) / 2,
     1.0],
)

shade_matrix = zeros(Float64, 24, ndays)
pos_matrix   = zeros(Float64, 24, ndays)
for m in 1:ndays
    shade_matrix[:, m] = month_shade[m] .* 100
    pos_matrix[:,   m] = month_pos[m]
end

p_shade = heatmap(1:ndays, hours, shade_matrix;
    color = :Greens, clims = (0, 100),
    colorbar_title = "%", title = "Shade selection (%)", ylabel = "hour", xlabel = "day of year")

p_pos = heatmap(1:ndays, hours, pos_matrix;
    color = pos_cmap, clims = pos_clims,
    colorbar_title = "cm (+ above, − below)",
    title = "Position (cm above/below ground)", ylabel = "hour", xlabel = "day of year")

display(plot(p_shade, p_pos; layout = (2, 1), size = (900, 600), left_margin = 6Plots.mm))

# =============================================================================
# Transient (thermal-mass-aware) body temperature — one-lump ODE model, run
# per representative day and contrasted with the discrete loop above, which
# assumes Tb reaches operative temperature instantly each hourly step.
# =============================================================================

# `simulate_onelump`/`simulate_transient_behavior` take `organism` directly, sourcing
# physics from its traits — same organism `thermoregulate()` uses above.

# One day's columns → a continuous-time forcing at the given height node (default 1, ground
# level). `climb_height_node` gives the ClimbPhase forcing below.
function diurnal_forcing(result, day_range, shade_fraction, height_node = 1)
    n = length(day_range)
    times = (0:n-1) .* 1.0u"hr" .|> u"s"
    EnvironmentForcing(times, EnvironmentalVarsVec(;
        air_temperature        = result.profile.air_temperature[day_range, height_node],
        sky_temperature        = result.sky_temperature[day_range],
        ground_temperature     = result.soil_temperature[day_range, 1],
        substrate_temperature  = result.soil_temperature[day_range, 1],
        relative_humidity      = result.profile.relative_humidity[day_range, height_node],
        wind_speed             = result.profile.wind_speed[day_range, height_node],
        atmospheric_pressure   = result.pressure[day_range],
        zenith_angle           = result.solar_radiation.zenith_angle[day_range],
        substrate_conductivity = result.soil_thermal_conductivity[day_range, 1],
        global_radiation       = result.global_radiation[day_range],
        diffuse_fraction       = result.diffuse_fraction[day_range],
        shade                  = fill(shade_fraction, n),
    ))
end

# BurrowPhase forcing at `burrow_depth_node` — mirrors `interpolate_environment`'s BELOWGROUND
# branch (thermoregulation.jl): no radiation, still air, everything set to soil temperature.
nearest_node(value, nodes) = argmin(abs.(nodes .- value))
burrow_depth_node = nearest_node(burrow_depth, depths)
climb_height_node = nearest_node(climb_height, heights)
function underground_forcing(result, day_range, depth_node)
    n = length(day_range)
    times = (0:n-1) .* 1.0u"hr" .|> u"s"
    soil_T = result.soil_temperature[day_range, depth_node]
    EnvironmentForcing(times, EnvironmentalVarsVec(;
        air_temperature        = soil_T,
        sky_temperature        = soil_T,
        ground_temperature     = soil_T,
        substrate_temperature  = soil_T,
        relative_humidity      = result.soil_humidity[day_range, depth_node],
        wind_speed              = fill(0.01u"m/s", n),
        atmospheric_pressure    = result.pressure[day_range],
        zenith_angle            = result.solar_radiation.zenith_angle[day_range],
        substrate_conductivity  = result.soil_thermal_conductivity[day_range, depth_node],
        global_radiation        = fill(0.0u"W/m^2", n),
        diffuse_fraction        = result.diffuse_fraction[day_range],
        shade                   = fill(1.0, n),
    ))
end

println("Running transient (one-lump) body temperature model for each representative day...")
transient_open   = Vector{Any}(undef, ndays)   # non-thermoregulating baseline (full sun)
transient_thermo = Vector{Any}(undef, ndays)   # thermoregulating (shuttles sun ⇄ shade)

# Movement, not just day/night, costs energy: active foraging and basking (postural
# adjustment, moving into/out of sun) scale metabolic heat production above resting baseline;
# unlisted phases (sleep/cool/climb/burrow/refuge) default to 1.0 (no scaling).
 metabolic_multipliers = (; forage = 1.0, bask = 1.0)

@time for m in 1:ndays
    day_range = (m-1)*24+1 : m*24
    times     = (0:23) .* 1.0u"hr" .|> u"s"
    sun_forcing        = diurnal_forcing(low_shade_result,  day_range, minimum_shade)
    shade_forcing      = diurnal_forcing(low_shade_result, day_range, maximum_shade)
    climb_forcing      = diurnal_forcing(low_shade_result,  day_range, minimum_shade, climb_height_node)
    underground_result = underground_forcing(low_shade_result, day_range, burrow_depth_node)
    core_temperature_init = low_shade_result.soil_temperature[day_range[1], burrow_depth_node]  # starts in BurrowPhase at midnight

    transient_open[m] = simulate_onelump(
        times, core_temperature_init, organism, env_pars, sun_forcing;
        posture = Intermediate(),
    )
    transient_thermo[m] = simulate_transient_behavior(
        times, core_temperature_init, organism, env_pars, sun_forcing, shade_forcing, limits;
        climb_forcing, underground_forcing = underground_result, metabolic_multipliers,
        cool_resume_margin = 2.0u"K", cool_resume_offset = 1.0u"K",
    )
    t_end_h  = ustrip(u"hr", transient_thermo[m].t[end])
    t_target_h = ustrip(u"hr", times[end])
    n_points = length(transient_thermo[m].t)
    if t_end_h < t_target_h - 0.1
        @warn "month $m ($(months[mod1(m,12)])): simulate_transient_behavior stopped at $(round(t_end_h; digits=2)) h (of $t_target_h h) after $n_points points — likely hit max_bouts"
    end
    n_climb  = count(p -> p isa ClimbPhase,  transient_thermo[m].phase)
    n_burrow = count(p -> p isa BurrowPhase, transient_thermo[m].phase)
    n_refuge = count(p -> p isa RefugePhase, transient_thermo[m].phase)
    (n_climb > 0 || n_burrow > 0 || n_refuge > 0) &&
        println("  month $m ($(months[mod1(m,12)])): climb=$n_climb, burrow=$n_burrow, refuge=$n_refuge accepted steps")
end

# ── Diagnostic: is CoolPhase exiting on the temperature-based resume signal
#    (which cool_resume_margin/cool_resume_offset control), or on sunset
#    (activity_signal) instead? If CoolPhase bouts consistently end right at
#    sunset, the shade site never cools the animal down to the resume
#    threshold before nightfall does the job instead — in that case
#    cool_resume_margin/cool_resume_offset have nothing left to affect.
function coolphase_end_hours(result)
    t_h    = ustrip.(u"hr", result.t)
    phases = result.phase
    ends   = Float64[]
    for i in eachindex(phases)
        if phases[i] isa CoolPhase && (i == lastindex(phases) || !(phases[i + 1] isa CoolPhase))
            push!(ends, t_h[i])
        end
    end
    return ends
end

println("\n── CoolPhase exit timing vs sunset (diagnostic) ──")
for m in 1:ndays
    zenith_day    = ustrip.(u"°", low_shade_result.solar_radiation.zenith_angle[(m-1)*24+1 : m*24])
    daytime_hours = hours[zenith_day .< 90]
    isempty(daytime_hours) && continue
    sunset_hr = maximum(daytime_hours)
    ends = coolphase_end_hours(transient_thermo[m])
    isempty(ends) && continue
    near_sunset = count(e -> abs(e - sunset_hr) < 0.5, ends)
    println("  month $m ($(months[mod1(m,12)])): $(length(ends)) CoolPhase bout(s), " *
            "$near_sunset ending within 0.5h of sunset ($(round(sunset_hr; digits=1))h)")
end

# ── Diagnostic: what does each CoolPhase bout hand off to? ForagePhase means
#    the resume signal won (cool_resume_margin/cool_resume_offset controlled
#    it); RefugePhase means the escape signal won instead (those two
#    parameters never got a chance to matter for that bout).
function coolphase_next_phases(result)
    phases = result.phase
    nexts  = Symbol[]
    for i in eachindex(phases)
        if phases[i] isa CoolPhase && (i == lastindex(phases) || !(phases[i + 1] isa CoolPhase))
            if i == lastindex(phases)
                push!(nexts, :end_of_sim)
            else
                next = phases[i + 1]
                push!(nexts, next isa RefugePhase ? :refuge : next isa ForagePhase ? :forage : :other)
            end
        end
    end
    return nexts
end

println("\n── CoolPhase exit reason (diagnostic) ──")
for m in 1:ndays
    nexts = coolphase_next_phases(transient_thermo[m])
    isempty(nexts) && continue
    n_forage = count(==(:forage), nexts)
    n_refuge = count(==(:refuge), nexts)
    n_other  = count(x -> x ∉ (:forage, :refuge), nexts)
    println("  month $m ($(months[mod1(m,12)])): $(length(nexts)) bout(s) — " *
            "resumed to forage=$n_forage, escalated to refuge=$n_refuge, other=$n_other")
end

# `simulate_onelump` saves only at `times` (hourly) — evaluate its dense ODE
# interpolant on a fine grid to see the sub-hourly dynamics it actually
# resolves. `simulate_transient_behavior` has no such restriction; its `.t`
# already carries every accepted adaptive solver step.
fine_hours   = range(0.0, 24.0; length = 288)
fine_times_s = fine_hours .* 3600.0
open_fine(m) = [ustrip(u"°C", transient_open[m].solution(t) * u"K") for t in fine_times_s]

transient_min = minimum(minimum(ustrip.(u"°C", r.core_temperature)) for r in transient_thermo)
transient_max = maximum(maximum(ustrip.(u"°C", r.core_temperature)) for r in transient_thermo)
open_min      = minimum(minimum(ustrip.(u"°C", r.core_temperature)) for r in transient_open)
open_max      = maximum(maximum(ustrip.(u"°C", r.core_temperature)) for r in transient_open)

println("\n── Transient vs discrete body temperature comparison ──")
println("  Discrete (instantaneous-equilibrium per step): $(round(minimum(T_body_C); digits=1)) to $(round(maximum(T_body_C); digits=1)) °C")
println("  Transient, thermoregulating (one-lump):        $(round(transient_min; digits=1)) to $(round(transient_max; digits=1)) °C")
println("  Transient, non-thermoregulating (open, sun):   $(round(open_min; digits=1)) to $(round(open_max; digits=1)) °C")

# ── Fig. 1 (transient) – Body temperature by month: steady-state and
#    transient overlaid on the same axes (4×3-ish grid) ────────────────────
panels_transient = map(1:ndays) do m
    p = plot(hours, month_Tb[m];
        lw = 2, color = :red, label = "steady-state",
        title = months[mod1(m, 12)], ylabel = "°C", ylim = (0, 55), titlefontsize = 9,
        legend = m == 1 ? :topright : false)
    plot!(p, ustrip.(u"hr", transient_thermo[m].t), ustrip.(u"°C", transient_thermo[m].core_temperature);
        lw = 2, color = :darkgreen, label = "transient (thermoreg)")
    plot!(p, fine_hours, open_fine(m);
        lw = 1, color = :grey, linestyle = :dot, label = "transient (open, no thermoreg)")
    hline!(p, [T_active_min_C, T_active_max_C];
        color = :orange, linestyle = :dash, lw = 1, label = "")
    p
end

display(plot(panels_transient...; layout = (ceil(Int, ndays / 3), 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Body temperature — steady-state vs transient (one-lump) models, $site_name"))

# ── Fig. 2 (transient) – Tb heatmap + active/basking-fraction heatmaps, side
#    by side with steady-state (cf. Fig. 2 above). Body temperature is
#    linearly interpolated onto the hourly grid; activity fractions are
#    integrated exactly over each hour bin (`hourly_state_fraction`), so a
#    bout that starts or ends mid-hour shows up as a graded value instead of
#    the steady-state loop's all-or-nothing per hour. ───────────────────────
tb_matrix_trans      = zeros(Float64, 24, ndays)
active_matrix_trans  = zeros(Float64, 24, ndays)
basking_matrix_trans = zeros(Float64, 24, ndays)
for m in 1:ndays
    t_h = ustrip.(u"hr", transient_thermo[m].t)
    y_C = ustrip.(u"°C", transient_thermo[m].core_temperature)
    tb_matrix_trans[:, m]      = interp_at_hours(t_h, y_C, hours)
    active_matrix_trans[:, m]  = hourly_state_fraction(t_h, transient_thermo[m].state, 24.0, Active,  hour_edges) .* 100
    basking_matrix_trans[:, m] = hourly_state_fraction(t_h, transient_thermo[m].state, 24.0, Basking, hour_edges) .* 100
end

p1_trans = heatmap(1:ndays, hours, tb_matrix_trans;
    color = cgrad(:RdYlBu, rev = true), clims = (0, 50),
    colorbar_title = "°C",
    title = "Body temperature (°C) — transient", ylabel = "hour", xlabel = "day of year")

p2a_trans = heatmap(1:ndays, hours, active_matrix_trans;
    color = :Reds, clims = (0, 100),
    colorbar_title = "%", title = "Active/foraging (%) — transient", ylabel = "hour", xlabel = "day of year")

p2b_trans = heatmap(1:ndays, hours, basking_matrix_trans;
    color = :Oranges, clims = (0, 100),
    colorbar_title = "%", title = "Basking (%) — transient", ylabel = "hour", xlabel = "day of year")

display(plot(p1, p1_trans, p2a, p2a_trans, p2b, p2b_trans; layout = (3, 2), size = (1400, 1200), left_margin = 6Plots.mm,
    plot_title = "Steady-state (left) vs transient (right)"))

# ── Fig. 3 (transient) – Shade selection and position, side by side with
#    steady-state (cf. Fig. 3 above). With climb/burrow enabled, `phase` gives
#    a transient equivalent of the steady-state position panel — one climb
#    height and one burrow depth, matching `heights`/`depths` above. Shade is
#    the fixed shade fraction each phase's forcing was built with (not a
#    continuously-chosen fraction, since each phase is one fixed site). ──────
phase_shade_frac(::SleepPhase) = maximum_shade
phase_shade_frac(::CoolPhase) = maximum_shade
phase_shade_frac(::BurrowPhase) = 1.0
phase_shade_frac(::RefugePhase) = 1.0
phase_shade_frac(::TransientBehavioralPhase) = minimum_shade  # Bask/Forage/Climb

# height_cm_max/depth_cm_max already computed above for the steady-state position panel.
phase_pos_cm(::ClimbPhase) = height_cm_max
phase_pos_cm(::BurrowPhase) = -depth_cm_max
phase_pos_cm(::RefugePhase) = -depth_cm_max
phase_pos_cm(::TransientBehavioralPhase) = 0.0

shade_matrix_trans = zeros(Float64, 24, ndays)
pos_matrix_trans   = zeros(Float64, 24, ndays)
for m in 1:ndays
    t_h          = ustrip.(u"hr", transient_thermo[m].t)
    phases_at_h  = state_at_hours(t_h, transient_thermo[m].phase, hours)
    shade_matrix_trans[:, m] = phase_shade_frac.(phases_at_h) .* 100
    pos_matrix_trans[:, m]   = phase_pos_cm.(phases_at_h)
end

p_shade_trans = heatmap(1:ndays, hours, shade_matrix_trans;
    color = :Greens, clims = (0, 100),
    colorbar_title = "%", title = "Shade selection (%) — transient (sun/shade phase)",
    ylabel = "hour", xlabel = "day of year")

p_pos_trans = heatmap(1:ndays, hours, pos_matrix_trans;
    color = pos_cmap, clims = pos_clims,
    colorbar_title = "cm (+ above, − below)",
    title = "Position (cm above/below ground) — transient", ylabel = "hour", xlabel = "day of year")

display(plot(p_shade, p_shade_trans, p_pos, p_pos_trans; layout = (2, 2), size = (1400, 800), left_margin = 6Plots.mm,
    plot_title = "Shade selection and position — steady-state (left) vs transient (right)"))

# =============================================================================
# Grid simulation — operative body temperature across a raster, with NO
# behavioural thermoregulation (cf. crucl2.jl's gridded microclimate example).
# This is the classic "operative temperature" map: the body temperature an
# ectotherm would equilibrate to if it couldn't behaviourally escape (no
# shade-seeking, no burrowing) — fixed full-sun exposure, fixed posture.
# =============================================================================

# ── Disk caching for the raster microclimate solves ────────────────────────
# The grid solves are the most memory-hungry step (a full multi-layer
# RasterStack per shade level, held alongside the derived output grids) —
# large enough areas can OOM even with the memory-freeing calls below. Set
# `USE_DISK_CACHE = true` to write each solved grid to disk and reload it on
# a later run instead of re-solving — trades disk space for peak RAM, and
# also skips the (slow) raster solve entirely on reruns of just the
# downstream plotting/analysis code.
const USE_DISK_CACHE = false
const CACHE_DIR = joinpath(@__DIR__, "cache")
isdir(CACHE_DIR) || mkpath(CACHE_DIR)

# Native (not "canonical" — see MicroclimateMapper's `canonical_unit`) unit
# for each layer in `output_layers` above. Zarr can't carry Unitful types, so
# every value is stripped to a plain Float64 before writing and reattached
# in the SAME unit (no lossy conversion) after reading back.
const GRID_LAYER_UNITS = (
    soil_temperature          = u"K",
    soil_moisture              = NoUnits,
    soil_thermal_conductivity = u"W/m/K",
    soil_humidity              = NoUnits,
    air_temperature            = u"K",
    relative_humidity          = NoUnits,
    wind_speed                 = u"m/s",
    global_radiation           = u"W/m^2",
    sky_temperature            = u"K",
    diffuse_fraction           = NoUnits,
    reference_temperature      = u"K",
    pressure                   = u"Pa",
    zenith_angle               = u"°",
)

strip_grid_units(stack) = RasterStack(NamedTuple{propertynames(stack)}(
    map(n -> ustrip.(GRID_LAYER_UNITS[n], stack[n]), propertynames(stack))))
reattach_grid_units(stack) = RasterStack(NamedTuple{propertynames(stack)}(
    map(n -> stack[n] .* GRID_LAYER_UNITS[n], propertynames(stack))))

# Solve a raster microclimate problem, or load a previously-cached copy from
# `CACHE_DIR/<cache_name>.zarr` if `USE_DISK_CACHE` and that file exists.
function solve_grid_cached(problem, cache_name)
    path = joinpath(CACHE_DIR, cache_name * ".zarr")
    if USE_DISK_CACHE && ispath(path)
        println("Loading cached grid microclimate: $cache_name")
        return reattach_grid_units(RasterStack(path))
    end
    println("Solving grid microclimate: $cache_name")
    cache  = init(problem)
    output = solve!(cache)
    cache  = nothing   # drop the per-worker scratch pool before the (large) write/GC below
    if USE_DISK_CACHE
        println("Caching grid microclimate to disk: $path")
        ispath(path) && rm(path; recursive = true)
        write(path, strip_grid_units(output))
    end
    GC.gc()
    return output
end

# ── Grid microclimate (reuses `model`/`dates` from the point run above) ───
# `geocode` returns a `GeocodeResult` with a `.extent` (buffered bounding box)
# alongside the point — use it directly as the raster area, as in crucl2.jl.
grid_site     = geocode(location; buffer = 2.0)
raster_area   = grid_site.extent
grid_template = load_template(WEATHER_SOURCE, grid_site)

@time grid_output = solve_grid_cached(MicroRasterProblem(;
    model, area = raster_area, dates,
    template     = grid_template,
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
), "grid_0pct_shade")

# ── Preallocate the output RasterStack ─────────────────────────────────────
# One layer per heat-budget output of interest, all sharing the grid's
# (X, Y, Ti) dims. NaN-filled so no-data pixels stay visibly distinct from a
# converged zero.
grid_shape = grid_output.soil_temperature[depth=1]   # (X, Y, Ti) template — shape/dims only
out_dims   = dims(grid_shape)
nx, ny, grid_nsteps = size(grid_shape)

_nan_layer(u) = Raster(fill(NaN * u, nx, ny, grid_nsteps), out_dims)

# Thread-safe progress reporting for the per-pixel grid loops below: an
# atomic counter shared across threads, logged periodically (not every
# pixel, to avoid interleaved output) with a simple elapsed-time ETA.
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

grid_result = RasterStack((;
    core_temperature       = _nan_layer(u"K"),
    metabolic_heat_flow    = _nan_layer(u"W"),
    evaporative_water_loss = _nan_layer(u"g/hr"),
    heat_balance_residual  = _nan_layer(u"W"),
))

# No-data pixels (e.g. ocean, or edge-of-coverage cells) come back as either
# `missing` or `NaN` depending on source/mask — `ismissing` short-circuits
# before `isnan` so this never has to compare `missing` in a boolean context.
_no_data(v) = ismissing(v) || isnan(ustrip(v))

# ── Per-pixel operative temperature (no thermoregulation loop) ────────────
# Fixed conditions since there's no behavioural loop to choose them: full sun
# (shade = 0), posture/absorptivity as already fixed in `organism`'s traits.
# `env_pars` is reused from the point run above — `elevation` isn't actually
# read by `heat_balance` (atmospheric pressure comes from `env_vars`, taken
# per-pixel from the grid below), so sharing it across pixels is safe.
fixed_shade = 0.0

# Wrapped in a function (rather than run as a top-level loop) so every
# variable it touches — `organism`, `env_pars`, array element types, etc. —
# is a concrete, specialised argument/local instead of an untyped global.
# Global-scope code can't be type-inferred in Julia, which both blocks
# compiler optimisation and forces heap-boxing of otherwise-cheap values;
# wrapping the hot loop in a function is the single biggest lever for both
# allocation and speed here. `:static` scheduling gives each thread a fixed,
# contiguous share of pixels (a good fit — the per-pixel workload is uniform),
# and each thread only ever writes to its own `grid_result.<layer>[x, y, :]`
# columns, so no locks/atomics/channels are needed for the writes themselves.
function solve_operative_grid!(grid_result, grid_output, organism, env_pars,
                                fixed_shade::Float64, nx::Int, ny::Int, grid_nsteps::Int)
    report_progress = grid_progress_reporter("Operative-temperature grid", nx * ny)
    Threads.@threads :static for xy in CartesianIndices((nx, ny))
        x, y = Tuple(xy)

        air_T_col = collect(grid_output.air_temperature[X(x), Y(y), height=1])
        if !_no_data(air_T_col[1])
            rh_col          = collect(grid_output.relative_humidity[X(x), Y(y), height=1])
            wind_col        = collect(grid_output.wind_speed[X(x), Y(y), height=1])
            sky_T_col       = collect(grid_output.sky_temperature[X(x), Y(y)])
            ground_T_col    = collect(grid_output.soil_temperature[X(x), Y(y), depth=1])
            substrate_k_col = collect(grid_output.soil_thermal_conductivity[X(x), Y(y), depth=1])
            pressure_col    = collect(grid_output.pressure[X(x), Y(y)])
            zenith_col      = collect(grid_output.zenith_angle[X(x), Y(y)])
            rad_col         = collect(grid_output.global_radiation[X(x), Y(y)])
            diffuse_col     = collect(grid_output.diffuse_fraction[X(x), Y(y)])

            @inbounds for step in 1:grid_nsteps
                env_vars = EnvironmentalVars(;
                    air_temperature        = air_T_col[step],
                    sky_temperature        = sky_T_col[step],
                    ground_temperature     = ground_T_col[step],
                    substrate_temperature  = ground_T_col[step],
                    relative_humidity      = clamp(rh_col[step], 0.0, 1.0),
                    wind_speed              = wind_col[step],
                    atmospheric_pressure    = pressure_col[step],
                    zenith_angle            = zenith_col[step],
                    substrate_conductivity  = substrate_k_col[step],
                    global_radiation        = rad_col[step],
                    diffuse_fraction        = diffuse_col[step],
                    shade                   = fixed_shade,
                )

                core_temperature = solve_body_temperature(organism, env_vars, env_pars)
                hb = heat_balance(core_temperature, organism, (; environment_pars = env_pars, environment_vars = env_vars))

                grid_result.core_temperature[x, y, step]       = hb.core_temperature
                grid_result.metabolic_heat_flow[x, y, step]    = hb.energy_balance.metabolic_heat_flow
                grid_result.evaporative_water_loss[x, y, step] = uconvert(u"g/hr",
                    hb.mass_balance.cutaneous_mass + hb.mass_balance.eye_mass + hb.mass_balance.respiration_mass)
                grid_result.heat_balance_residual[x, y, step]  = hb.heat_balance
            end
        end

        report_progress()
    end
    return grid_result
end

println("Solving operative body temperature grid...")
@time solve_operative_grid!(grid_result, grid_output, organism, env_pars, fixed_shade, nx, ny, grid_nsteps)

# ── Plots: annual max/min and per-month max/min, all three grid variables ──
# NaN-aware max/min along a dimension (skip no-data pixels), operating on
# already-stripped (unitless) data — reduction can't propagate `-Inf`/`Inf`
# through Unitful arithmetic the way it can through plain Float64.
nanmax(A; dims) = (m = dropdims(maximum(ifelse.(isnan.(A), -Inf, A); dims); dims); ifelse.(isinf.(m), NaN, m))
nanmin(A; dims) = (m = dropdims(minimum(ifelse.(isnan.(A),  Inf, A); dims); dims); ifelse.(isinf.(m), NaN, m))
nanextrema(A)   = (vals = filter(!isnan, vec(parent(A))); isempty(vals) ? (0.0, 1.0) : (minimum(vals), maximum(vals)))

# One annual heatmap (max or min) for a stripped (X, Y, Ti) grid.
function annual_extreme_plot(A_stripped, reducer, label, unit, palette)
    m = reducer(A_stripped; dims = Ti) .* unit
    display(plot(m; title = label, seriestype = :heatmap, color = palette, yflip = false))
end

# One representative-day-per-month small-multiples grid (max or min), all
# panels sharing a single colour scale so months are visually comparable.
function monthly_extreme_plot(A_stripped, reducer, ndays, label, unit, palette, month_names)
    maps = [reducer(A_stripped[Ti((m-1)*24+1 : m*24)]; dims = Ti) for m in 1:ndays]
    clim_lo, clim_hi = nanextrema(reduce(vcat, vec.(parent.(maps))))
    panels = [plot(maps[m] .* unit; title = month_names[mod1(m, 12)], seriestype = :heatmap,
                   color = palette, clims = (clim_lo, clim_hi) .* unit, colorbar = false,
                   axis = false, ticks = false, titlefontsize = 9)
              for m in 1:ndays]
    display(plot(panels...; layout = (ceil(Int, ndays / 4), 4),
                 size = (1200, 260 * ceil(Int, ndays / 4)), plot_title = label))
end

grid_ndays  = grid_nsteps ÷ 24
month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

Tb_C  = ustrip.(u"°C",   grid_result.core_temperature)
MHF_W = ustrip.(u"W",    grid_result.metabolic_heat_flow)
EWL_g = ustrip.(u"g/hr", grid_result.evaporative_water_loss)

Tb_palette  = cgrad([:blue, :lightblue, :orange, :red, :purple])
MHF_palette = cgrad([:black, :orange, :red])
EWL_palette = cgrad([:white, :teal, :blue])

# Kept as a named variable (not just plotted inline) — reused below for the
# thermoregulation cooling-benefit comparison plot.
max_Tb = nanmax(Tb_C; dims = Ti) .* u"°C"

for (A, unit, palette, name) in (
    (Tb_C,  u"°C",   Tb_palette,  "operative body temperature"),
    (MHF_W, u"W",    MHF_palette, "metabolic heat flow"),
    (EWL_g, u"g/hr", EWL_palette, "evaporative water loss"),
)
    annual_extreme_plot(A, nanmax, "Annual max $name (no thermoregulation)", unit, palette)
    annual_extreme_plot(A, nanmin, "Annual min $name (no thermoregulation)", unit, palette)
    monthly_extreme_plot(A, nanmax, grid_ndays, "Monthly max $name (no thermoregulation)", unit, palette, month_names)
    monthly_extreme_plot(A, nanmin, grid_ndays, "Monthly min $name (no thermoregulation)", unit, palette, month_names)
end

# `grid_result` (4 full (X, Y, Ti) layers) and the stripped copies made for
# plotting are only needed for the plots just displayed — `max_Tb` above
# already captured the one summary value reused later. Drop the rest before
# the thermoregulation section's fresh 90%-shade raster solve.
grid_result = nothing
Tb_C = MHF_W = EWL_g = nothing
GC.gc()

# =============================================================================
# Grid simulation WITH thermoregulation — the same behavioural loop as the
# point run above (shade-seeking, burrowing, basking), applied per pixel
# across the raster. Contrast against the operative-temperature grid above:
# thermoregulation should narrow the annual Tb range.
# =============================================================================

# 0% shade is identical to the operative-temperature grid's microclimate run
# above (same model/area/dates/template/soil_profile/init, no shade override)
# — reuse it instead of re-solving the whole raster a second time.
grid_low_output = grid_output

@time grid_high_output = solve_grid_cached(MicroRasterProblem(;
    model, area = raster_area, dates,
    template     = grid_template,
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
    data         = (; shade = maximum_shade),
), "grid_90pct_shade")

# ── Adapter: one pixel's columns from a grid RasterStack → the plain
#    Matrix/Vector shape AvailableEnvironments expects (cf. `point_environment`).
function grid_pixel_environment(output, x, y)
    profile = (;
        air_temperature   = collect(output.air_temperature[X(x), Y(y)]),
        relative_humidity = collect(output.relative_humidity[X(x), Y(y)]),
        wind_speed        = collect(output.wind_speed[X(x), Y(y)]),
    )
    return (;
        pressure                  = collect(output.pressure[X(x), Y(y)]),
        reference_temperature     = collect(output.reference_temperature[X(x), Y(y)]),
        global_radiation          = collect(output.global_radiation[X(x), Y(y)]),
        diffuse_fraction          = collect(output.diffuse_fraction[X(x), Y(y)]),
        sky_temperature           = collect(output.sky_temperature[X(x), Y(y)]),
        soil_temperature          = collect(output.soil_temperature[X(x), Y(y)]),
        soil_humidity             = collect(output.soil_humidity[X(x), Y(y)]),
        soil_thermal_conductivity = collect(output.soil_thermal_conductivity[X(x), Y(y)]),
        profile,
        solar_radiation = (; zenith_angle = collect(output.zenith_angle[X(x), Y(y)])),
    )
end

# ── Preallocate the thermoregulated output RasterStack ─────────────────────
_nan_layer_plain() = Raster(fill(NaN, nx, ny, grid_nsteps), out_dims)

grid_thermoreg_result = RasterStack((;
    core_temperature = _nan_layer(u"K"),
    state            = _nan_layer_plain(),   # 0 = resting, 1 = basking, 2 = active
    shade            = _nan_layer_plain(),
))

# ── Per-pixel thermoregulation loop ────────────────────────────────────────
# `organism`/`limits`/`env_pars` are immutable and shared read-only across
# threads; `previous_depth`/`activity_commenced` are pixel-local state, so
# each thread's column is fully independent — same lock-free write pattern
# as the operative-temperature grid above. Wrapped in a function for the same
# reason as `solve_operative_grid!` above: global-scope code can't be
# specialised/inferred by the compiler, so every access to `organism`,
# `limits`, etc. from a top-level loop is a slow, boxing dynamic lookup.
function solve_thermoreg_grid!(grid_thermoreg_result, grid_low_output, grid_high_output,
                                organism, limits, env_pars,
                                minimum_shade::Float64, maximum_shade::Float64,
                                depths, heights, nx::Int, ny::Int, grid_nsteps::Int)
    report_progress = grid_progress_reporter("Thermoregulation grid", nx * ny)
    Threads.@threads :static for xy in CartesianIndices((nx, ny))
        x, y = Tuple(xy)

        low_env = grid_pixel_environment(grid_low_output, x, y)
        if !_no_data(low_env.profile.air_temperature[1, 1])
            high_env = grid_pixel_environment(grid_high_output, x, y)

            available_environments_px = AvailableEnvironments(
                low_env, high_env, minimum_shade, maximum_shade, depths, heights)

            previous_depth     = limits.depth.reference
            activity_commenced = false

            for step in 1:grid_nsteps
                if (step - 1) % 24 == 0
                    activity_commenced = false
                end
                out = thermoregulate(
                    organism, available_environments_px, limits, env_pars, step, previous_depth;
                    activity_commenced,
                )
                previous_depth     = out.depth_node
                activity_commenced = activity_commenced || out.state isa Active || out.state isa Basking

                grid_thermoreg_result.core_temperature[x, y, step] = out.core_temperature
                grid_thermoreg_result.state[x, y, step] = out.state isa Active ? 2.0 : out.state isa Basking ? 1.0 : 0.0
                grid_thermoreg_result.shade[x, y, step] = out.shade
            end
        end

        report_progress()
    end
    return grid_thermoreg_result
end

println("Solving grid thermoregulation...")
solve_thermoreg_grid!(grid_thermoreg_result, grid_low_output, grid_high_output,
                       organism, limits, env_pars, minimum_shade, maximum_shade,
                       depths, heights, nx, ny, grid_nsteps)

# ── Plots: annual max/min thermoregulated body temperature ─────────────────
Tb_C_reg    = ustrip.(u"°C", grid_thermoreg_result.core_temperature)
Tb_C_reg_hi = ifelse.(isnan.(Tb_C_reg), -Inf, Tb_C_reg)
Tb_C_reg_lo = ifelse.(isnan.(Tb_C_reg),  Inf, Tb_C_reg)
max_Tb_reg = dropdims(maximum(Tb_C_reg_hi; dims = Ti); dims = Ti)
min_Tb_reg = dropdims(minimum(Tb_C_reg_lo; dims = Ti); dims = Ti)
max_Tb_reg = ifelse.(isinf.(max_Tb_reg), NaN, max_Tb_reg) .* u"°C"
min_Tb_reg = ifelse.(isinf.(min_Tb_reg), NaN, min_Tb_reg) .* u"°C"

display(plot(max_Tb_reg;
    title = "Annual max body temperature (WITH thermoregulation)",
    seriestype = :heatmap, color = Tb_palette, yflip = false))

display(plot(min_Tb_reg;
    title = "Annual min body temperature (WITH thermoregulation)",
    seriestype = :heatmap, color = Tb_palette, yflip = false))

# Thermoregulation's cooling benefit: how much lower is annual max Tb than
# the operative (behaviour-free) temperature computed above?
Tb_buffering = max_Tb .- max_Tb_reg
display(plot(Tb_buffering;
    title = "Thermoregulation cooling benefit (max Tb reduction, °C)",
    seriestype = :heatmap, color = cgrad([:white, :blue]), yflip = false))

# Fraction of the year spent active, per pixel.
active_fraction = dropdims(
    mean(grid_thermoreg_result.state .== 2.0; dims = Ti); dims = Ti) .* 100

display(plot(active_fraction;
    title = "Time spent active (% of simulated hours)",
    seriestype = :heatmap, color = cgrad([:steelblue, :orange, :firebrick]), yflip = false))

# =============================================================================
# Grid simulation WITH transient (thermal-mass-aware) thermoregulation — the
# same `simulate_transient_behavior` model used at the point level above,
# applied per pixel across the raster. By far the most expensive step in this
# file: each pixel runs a full adaptive-step ODE solve (several bouts, each
# its own `solve()` call) per representative day, vs. one cheap discrete
# hourly step for the steady-state grid above. Only per-day summary stats
# (max/min core temperature, active fraction) are kept per pixel — not the
# full sub-daily trajectory — to keep memory bounded.
# =============================================================================

day_out_dims = (out_dims[1], out_dims[2], Ti(1:grid_ndays))
_nan_day_layer(u) = Raster(fill(NaN * u, nx, ny, grid_ndays), day_out_dims)

grid_transient_result = RasterStack((;
    core_temperature_max = _nan_day_layer(u"K"),
    core_temperature_min = _nan_day_layer(u"K"),
    active_fraction       = _nan_day_layer(NoUnits),
))

# `grid_low_output`/`grid_high_output` (0%/90% shade grids) are reused from the
# steady-state grid section above — no re-solve needed.
function solve_transient_grid!(grid_transient_result, grid_low_output, grid_high_output,
                                organism, limits, env_pars,
                                minimum_shade::Float64, maximum_shade::Float64,
                                climb_height_node::Int, burrow_depth_node::Int,
                                metabolic_multipliers::NamedTuple, hour_edges::Vector{Float64},
                                nx::Int, ny::Int, grid_ndays::Int)
    report_progress = grid_progress_reporter("Transient grid", nx * ny)
    Threads.@threads :static for xy in CartesianIndices((nx, ny))
        x, y = Tuple(xy)

        low_env = grid_pixel_environment(grid_low_output, x, y)
        if !_no_data(low_env.profile.air_temperature[1, 1])
            high_env = grid_pixel_environment(grid_high_output, x, y)

            for m in 1:grid_ndays
                day_range = (m-1)*24+1 : m*24
                times     = (0:23) .* 1.0u"hr" .|> u"s"
                sun_forcing         = diurnal_forcing(low_env,  day_range, minimum_shade)
                shade_forcing       = diurnal_forcing(high_env, day_range, maximum_shade)
                climb_forcing       = diurnal_forcing(low_env,  day_range, minimum_shade, climb_height_node)
                underground_result  = underground_forcing(high_env, day_range, burrow_depth_node)
                core_temperature_init = high_env.soil_temperature[day_range[1], burrow_depth_node]

                sol = simulate_transient_behavior(
                    times, core_temperature_init, organism, env_pars, sun_forcing, shade_forcing, limits;
                    climb_forcing, underground_forcing = underground_result, metabolic_multipliers,
                )

                grid_transient_result.core_temperature_max[x, y, m] = maximum(sol.core_temperature)
                grid_transient_result.core_temperature_min[x, y, m] = minimum(sol.core_temperature)
                grid_transient_result.active_fraction[x, y, m] = mean(hourly_state_fraction(
                    ustrip.(u"hr", sol.t), sol.state, 24.0, Active, hour_edges))
            end
        end

        report_progress()
    end
    return grid_transient_result
end

println("Solving grid transient thermoregulation...")
@time solve_transient_grid!(grid_transient_result, grid_low_output, grid_high_output,
    organism, limits, env_pars, minimum_shade, maximum_shade,
    climb_height_node, burrow_depth_node, metabolic_multipliers, hour_edges, nx, ny, grid_ndays)

# ── Plots: annual max/min transient body temperature ───────────────────────
max_Tb_transient = nanmax(ustrip.(u"°C", grid_transient_result.core_temperature_max); dims = Ti) .* u"°C"
min_Tb_transient = nanmin(ustrip.(u"°C", grid_transient_result.core_temperature_min); dims = Ti) .* u"°C"
active_fraction_transient = dropdims(
    mean(grid_transient_result.active_fraction; dims = Ti); dims = Ti) .* 100

display(plot(
    plot(max_Tb_reg;         title = "Annual max Tb — steady-state",  seriestype = :heatmap, color = Tb_palette, yflip = false),
    plot(max_Tb_transient;   title = "Annual max Tb — transient",     seriestype = :heatmap, color = Tb_palette, yflip = false),
    plot(min_Tb_reg;         title = "Annual min Tb — steady-state",  seriestype = :heatmap, color = Tb_palette, yflip = false),
    plot(min_Tb_transient;   title = "Annual min Tb — transient",     seriestype = :heatmap, color = Tb_palette, yflip = false);
    layout = (2, 2), size = (1000, 800), left_margin = 6Plots.mm,
    plot_title = "Steady-state vs transient body temperature extremes"))

# Does modelling thermal mass narrow or widen the simulated extremes relative
# to the instantaneous-equilibrium steady-state assumption?
display(plot(
    plot(max_Tb_transient .- max_Tb_reg;
        title = "Max Tb: transient − steady-state", seriestype = :heatmap,
        color = cgrad([:blue, :white, :red]), yflip = false),
    plot(min_Tb_transient .- min_Tb_reg;
        title = "Min Tb: transient − steady-state", seriestype = :heatmap,
        color = cgrad([:blue, :white, :red]), yflip = false);
    layout = (1, 2), size = (1000, 400), left_margin = 6Plots.mm,
    plot_title = "Transient − steady-state (°C)"))

display(plot(
    plot(active_fraction;             title = "Time active (%) — steady-state", seriestype = :heatmap,
        color = cgrad([:steelblue, :orange, :firebrick]), yflip = false),
    plot(active_fraction_transient;   title = "Time active (%) — transient",    seriestype = :heatmap,
        color = cgrad([:steelblue, :orange, :firebrick]), yflip = false);
    layout = (1, 2), size = (1000, 400), left_margin = 6Plots.mm,
    plot_title = "Fraction of time active — steady-state vs transient"))
