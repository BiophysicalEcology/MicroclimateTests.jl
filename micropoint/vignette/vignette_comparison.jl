# vignette_comparison.jl — reproduces micropoint's own published vignette
# forest example (see run_micropoint_vignette.R) with an equivalent
# Microclimate.jl MultilayerCanopy run driven from the *exact same* hourly
# forcing, and compares a reasonable subset of outputs + solve-time
# benchmark. See README.md for the full parameter-mapping/known-mismatch
# writeup — this is a cross-implementation physics sanity check, not a
# regression test with a pass/fail threshold.
#
# depths uses Microclimate.jl's own DEFAULT_DEPTHS, not micropoint's
# geometricCpp node placement -- geometricCpp's ~2.4mm first node forced the
# soil ODE solver into a slow, marginally-accurate regime near the surface,
# which fed an overheated ground_temperature into the canopy's longwave
# exchange every hour (day and night alike). That was the actual cause of
# both the 40x slowdown and the persistent leaf-warmer-than-air offset seen
# in earlier versions of this script -- not cloud_cover, not the longwave
# input pathway, not the canopy's own leaf parameters (all independently
# ruled out first). DEFAULT_DEPTHS (1.25cm first node) still resolves the
# near-surface diurnal wave (a uniform grid does not) without that cost.
#
# Run once run_micropoint_vignette.R has produced outputs/:
#   cd c:/git/BiophysicalEcologyEnv
#   julia --project=. c:/git/MicroclimateTests.jl/micropoint/vignette/vignette_comparison.jl

using Microclimate, Unitful
using FluidProperties: wet_air_properties
using CSV, DataFrames, Statistics, Printf, Plots

gr()

outdir = joinpath(@__DIR__, "outputs")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  READ micropoint's outputs                                               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

climdata = DataFrame(CSV.File(joinpath(outdir, "climdata.csv")))
paii_df  = DataFrame(CSV.File(joinpath(outdir, "paii.csv")))
soilc    = DataFrame(CSV.File(joinpath(outdir, "soilc.csv")))
params   = DataFrame(CSV.File(joinpath(outdir, "params.csv")))[1, :]
timing_r = DataFrame(CSV.File(joinpath(outdir, "timing.csv")))[1, :]
soil_depths_df = DataFrame(CSV.File(joinpath(outdir, "soil_depths.csv")))

read_matrix(name) = Matrix(DataFrame(CSV.File(joinpath(outdir, name * ".csv"))))

r_tair      = read_matrix("tair")       # °C, [8760 x 20], col 1 = ground, col 20 = canopy top
r_tleaf     = read_matrix("tleaf")      # °C
r_relhum    = read_matrix("relhum")     # %
r_windspeed = read_matrix("windspeed")  # m/s
r_Rdirdown  = read_matrix("Rdirdown")   # W/m^2, beam-normal (not horizontal)
r_Rdifdown  = read_matrix("Rdifdown")   # W/m^2
r_Rswup     = read_matrix("Rswup")      # W/m^2
r_Rlwdown   = read_matrix("Rlwdown")    # W/m^2
r_Rlwup     = read_matrix("Rlwup")      # W/m^2
r_tsoil     = read_matrix("tsoil")      # °C, [8760 x 16], col 1 = surface, col 16 = deepest
r_theta     = read_matrix("theta")      # m^3/m^3, [8760 x 16], col 1 = surface, col 16 = deepest

nhours = nrow(climdata)
n_layers = 20

# micropoint's canopy-layer columns run bottom (1) to top (20) — Microclimate.jl's
# run top (1) to bottom (20). Flip R's columns once here so every subsequent
# comparison lines up index-for-index against the Julia output.
flip(m) = m[:, end:-1:1]
r_tair, r_tleaf, r_relhum, r_windspeed = flip(r_tair), flip(r_tleaf), flip(r_relhum), flip(r_windspeed)
r_Rdirdown, r_Rdifdown, r_Rswup, r_Rlwdown, r_Rlwup =
    flip(r_Rdirdown), flip(r_Rdifdown), flip(r_Rswup), flip(r_Rlwdown), flip(r_Rlwup)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  BUILD an equivalent Microclimate.jl MultilayerCanopy run                ║
# ╚══════════════════════════════════════════════════════════════════════════╝

lat_deg, lon_deg = params.lat, params.long
canopy_height = params.canopy_height * 1.0u"m"

# 20 evenly-spaced canopy layers (micropoint's RunModelFull doesn't return
# per-layer heights to verify against — this is a documented assumption, see
# README.md), plus one height above canopy top (matching climdata's own
# weatherhgt_adjust(zout=27) — the height micropoint's own forcing is valid at).
canopy_heights = collect(1:n_layers) .* (ustrip(u"m", canopy_height) / n_layers) .* u"m"
heights = vcat(canopy_heights, [27.0u"m"])

# paii/Lfrac are bottom-to-top in R; Microclimate.jl's plant_area_index is
# top-to-bottom — reverse here, once, at translation time.
plant_area_index = reverse(paii_df.paii) .* 1.0

site = Site(;
    latitude = lat_deg * 1.0u"°",
    longitude = lon_deg * 1.0u"°",
    elevation = 0.0u"m",
    slope = 0.0u"°",
    aspect = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24),
    sky_view_fraction = 1.0,
    albedo = params.gref,
    roughness_height = 0.004u"m",
    atmospheric_pressure = mean(climdata.pres) * 1.0u"kPa" |> p -> uconvert(u"Pa", p),
)

# Soil: Clay loam (createsoilc defaults) translated via the already-verified
# psi_e (m -> J/kg, x9.81) and Ksat (kg*s/m^3, exact unit match) conversions
# from the sibling micropoint/config.jl bare-ground comparison. Smax/Smin/n/
# Vq/Vm/Vo/Mc/rho have no Microclimate.jl equivalent (different retention-
# curve formulation) -- not translated, see README.md.
# Microclimate.jl's own default depth grid, not micropoint's geometricCpp
# nodes (first step ~2.4mm) -- that fine a near-surface spacing forces the
# soil ODE solver into pathologically small steps (444s vs ~10s for a
# comparable canopy_site_check.jl run). DEFAULT_DEPTHS (first step 1.25cm)
# still resolves the near-surface diurnal wave (a uniform grid does not,
# confirmed earlier) without that cost -- same grid the rest of the codebase
# already validates against.
depths = Microclimate.DEFAULT_DEPTHS  # not exported

psi_e_Jkg = mean(soilc.psi_e) * 9.81 * -1.0u"J/kg"
Ksat = mean(soilc.Ksat) * 1.0u"kg*s/m^3"
campbell_b = mean(soilc.b) * 1.0

soil_profile = example_soil_profile(depths;
    hydraulics = example_campbell_hydraulic_profile(depths;
        air_entry_water_potential = fill(psi_e_Jkg, length(depths)),
        saturated_hydraulic_conductivity = fill(Ksat, length(depths)),
        campbell_b_parameter = fill(campbell_b, length(depths)),
    ),
)
soil_hydraulic_model = example_soil_hydraulic_model()

canopy_model = MultilayerCanopy(;
    canopy_height, plant_area_index, woody_area_fraction = 0.0,
    shortwave_model = TwoStreamRadiation(; leaf_reflectance = 0.4, leaf_transmittance = 0.2),  # vegp$lref/ltra
    interception_model = LayeredRainInterception(; leaf_water_storage_capacity = 0.2u"kg/m^2"),
    leaf_parameters = LeafParameters(;
        leaf_length = 0.12u"m", leaf_width = 0.05u"m",
        leaf_emissivity = 0.97, leaf_angle_distribution_parameter = 1.6,
    ),
    convergence_model = PicardCanopyConvergence(;
    convergence=IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=200), relaxation=0.7),
)

environment_daily = DailyTimeseries(;
    shade = fill(0.0, 365),   # ignored by MultilayerCanopy (own radiative transfer)
    soil_wetness = fill(0.0, 365),
    surface_emissivity = fill(params.groundem, 365),
    cloud_emissivity = fill(params.groundem, 365),
    rainfall = [sum(climdata.precip[(d-1)*24+1:d*24]) for d in 1:365] .* 1.0u"kg/m^2",
    deep_soil_temperature = fill((mean(climdata.temp) + 273.15) * 1.0u"K", 365),
    leaf_area_index = fill(sum(plant_area_index), 365),
)

ref_temp_K = (climdata.temp .+ 273.15) .* 1.0u"K"
ref_pres_Pa = uconvert.(u"Pa", climdata.pres .* 1.0u"kPa")
ref_humidity = clamp.(climdata.relhum, 0.0, 100.0) ./ 100.0

# precompute_longwave_sky (canopy_longwave!'s sky term) ignores
# environment_instant.longwave_radiation entirely -- it always derives its
# own clear-sky estimate (CampbellNormanAtmosphericRadiation) blended with a
# synthetic cloud term via cloud_cover. A flat placeholder cloud_cover would
# silently overstate sky longwave on real clear nights (suppressing
# radiative cooling) and understate it on real cloudy ones. Invert that same
# blend against climdata's real measured lwdown instead, so the canopy sees
# the actual sky condition each hour.
const SIGMA_SB = 5.670374419e-8u"W/m^2/K^4"
cloud_em = params.groundem
vapour_pressure = [wet_air_properties(ref_temp_K[i], ref_humidity[i], ref_pres_Pa[i]; vapour_pressure_equation = GoffGratch()).vapour_pressure
                    for i in 1:nhours]
atmospheric_lw = [Microclimate.atmospheric_radiation(CampbellNormanAtmosphericRadiation(), vapour_pressure[i], ref_temp_K[i]) for i in 1:nhours]
cloud_lw = SIGMA_SB .* cloud_em .* (ref_temp_K .- 2.0u"K") .^ 4
measured_lw = climdata.lwdown .* 1.0u"W/m^2"
cloud_cover_derived = clamp.(ustrip.(NoUnits, (measured_lw .- atmospheric_lw) ./ (cloud_lw .- atmospheric_lw)), 0.0, 1.0)

# environment_hourly = HourlyTimeseries(;
#     pressure = ref_pres_Pa,
#     reference_temperature = ref_temp_K,
#     reference_humidity = ref_humidity.* 0.0 .+ 1.0,
#     reference_wind_speed = clamp.(climdata.windspeed, 0.1, Inf) .* 1.0u"m/s",
#     global_radiation = climdata.swdown .* 0.0u"W/m^2",
#     longwave_radiation = 5.687e-8 .* ustrip.(u"K", ref_temp_K).^4 .* 1.0u"W/m^2",
#     cloud_cover = cloud_cover_derived .* 0.0 .+ 1.0,
#     rainfall = climdata.precip .* 0.0u"kg/m^2",
#     zenith_angle = nothing,            # computed internally from site/day/hour
# )
environment_hourly = HourlyTimeseries(;
    pressure = ref_pres_Pa,
    reference_temperature = ref_temp_K,
    reference_humidity = ref_humidity,
    reference_wind_speed = clamp.(climdata.windspeed, 0.1, Inf) .* 1.0u"m/s",
    global_radiation = climdata.swdown .* 1.0u"W/m^2",
    longwave_radiation = climdata.lwdown .* 1.0u"W/m^2",
    cloud_cover = cloud_cover_derived,
    rainfall = climdata.precip .* 1.0u"kg/m^2",
    zenith_angle = nothing,            # computed internally from site/day/hour
)
radiation = RadiationModel()

model = MicroModel(;
    hours = 0.0:1.0:23.0,
    depths, heights,
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model,
    radiation,
    canopy_model,
    config = MicroConfig(; rainfall_schedule = HourlyRainfall(), soil_moisture_strategy = DynamicSoilMoisture(;
        #moisture_tolerance = microinput[:IM]u"kg/m^2/s",
        #moisture_max_iterations = 36,
        #moisture_timestep = 50u"s",
        )),
)

inputs = MicroInputs(;
    site, soil_profile,
    environment_minmax = nothing,
    environment_daily, environment_hourly,
    initial_soil_temperature = fill(u"K"((mean(climdata.temp) + 273.15) * 1.0u"K"), length(depths)),
    initial_soil_moisture = fill(0.3, length(depths)),
)

problem = MicroProblem(model, inputs; days = collect(1.0:365.0), time_mode = ConsecutiveDayMode())

println("Solving Microclimate.jl (MultilayerCanopy, $n_layers layers, $nhours hours)...")
julia_elapsed = @elapsed (jl_out = Microclimate.solve(problem))
@printf("Microclimate.jl solve: %.2f s\n", julia_elapsed)
@printf("micropoint RunModelFull: %.2f s\n", timing_r.elapsed_s)
@printf("%s is %.1fx faster on this run\n",
    julia_elapsed <= timing_r.elapsed_s ? "Microclimate.jl" : "micropoint",
    max(timing_r.elapsed_s / julia_elapsed, julia_elapsed / timing_r.elapsed_s))

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  DIAGNOSE: is the leaf-air gap explained by a nonzero energy-balance      ║
# ║  residual (a solver bug), or is net_balance ≈ 0 and the gap comes from    ║
# ║  the inputs (absorbed_radiation) feeding a correctly-solved leaf?         ║
# ╚══════════════════════════════════════════════════════════════════════════╝

jl_net = ustrip.(u"W", jl_out.canopy.net_balance)
jl_absorbed_rad = ustrip.(u"W/m^2", jl_out.canopy.absorbed_radiation)
jl_gap = ustrip.(u"K", jl_out.canopy.leaf_temperature .- jl_out.canopy.air_temperature)

worst_net_idx = argmax(abs.(jl_net))
worst_gap_idx = argmax(abs.(jl_gap))
println("\n── Diagnostic: energy-balance residual vs. leaf-air gap ──")
println("Max |net_balance| across all hours/layers: ", jl_net[worst_net_idx], " W",
    " at step=", worst_net_idx[1], " layer=", worst_net_idx[2],
    " (gap there = ", round(jl_gap[worst_net_idx]; digits=2), " K)")
println("Max |leaf - air| gap across all hours/layers: ", round(jl_gap[worst_gap_idx]; digits=2), " K",
    " at step=", worst_gap_idx[1], " layer=", worst_gap_idx[2],
    " (net_balance there = ", jl_net[worst_gap_idx], " W,",
    " absorbed_radiation = ", round(jl_absorbed_rad[worst_gap_idx]; digits=1), " W/m^2)")

println("\n-- Top 10 largest |leaf - air| gaps: step, layer, gap_K, net_balance_W, absorbed_radiation_Wm2 --")
gap_order = sortperm(vec(abs.(jl_gap)); rev=true)[1:10]
for lin in gap_order
    step, layer = Tuple(CartesianIndices(jl_gap)[lin])
    @printf("  step=%5d layer=%2d gap=%6.2f K  net=%8.3f W  absorbed_radiation=%7.1f W/m^2  (%s)\n",
        step, layer, jl_gap[step, layer], jl_net[step, layer], jl_absorbed_rad[step, layer], climdata.obs_time[step])
end

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  COMPARE                                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

struct Stat; r::Float64; rmse::Float64; bias::Float64; end
function compute_stats(a, b)
    m = isfinite.(a) .& isfinite.(b)
    count(m) < 2 && return Stat(NaN, NaN, NaN)
    Stat(cor(a[m], b[m]), sqrt(mean((b[m] .- a[m]) .^ 2)), mean(b[m] .- a[m]))
end
fmt(s::Stat) = @sprintf("r=%+.3f  RMSE=%6.3f  bias=%+7.4f", s.r, s.rmse, s.bias)

jl_air  = ustrip.(u"°C", jl_out.canopy.air_temperature)
jl_leaf = ustrip.(u"°C", jl_out.canopy.leaf_temperature)
jl_rh   = jl_out.canopy.relative_humidity .* 100.0
jl_wind = ustrip.(u"m/s", jl_out.canopy.wind_speed)

println("\n── Canopy top vs. bottom (Microclimate.jl vs. micropoint, hourly, direct) ──")
# Both r_* and jl_* now follow the same top(1)-to-bottom(n_layers) convention
# after the flip() above -- same column index means the same physical layer.
for (label, jl_col, r_col) in [("top", 1, 1), ("bottom", n_layers, n_layers)]
    println("  Layer: $label")
    for (name, jlm, rm) in [
        ("Air temp (°C)",  jl_air,  r_tair),
        ("Leaf temp (°C)", jl_leaf, r_tleaf),
        ("Rel. humidity (%)", jl_rh, r_relhum),
        ("Wind speed (m/s)", jl_wind, r_windspeed),
    ]
        s = compute_stats(rm[:, r_col], jlm[:, jl_col])
        println("    " * rpad(name, 20) * fmt(s))
    end
end

println("\n── Downward longwave: canopy top vs. near-ground (direct) ──")
# r_Rlwdown is per-CANOPY-LAYER, not literally at the ground -- column n_layers
# (post-flip) is the lowest canopy layer, the closest available proxy to
# ground-level incoming longwave, not an exact match.
jl_downLW_all = ustrip.(u"W/m^2", jl_out.canopy.boundary_downward_longwave)
jl_ground_incoming_lw = ustrip.(u"W/m^2", jl_out.canopy.ground_absorbed_longwave) ./ params.groundem
for (label, jlv, rv) in [
    ("canopy top",  jl_downLW_all[:, 1], r_Rlwdown[:, 1]),
    ("near-ground (ground vs. lowest layer)", jl_ground_incoming_lw, r_Rlwdown[:, n_layers]),
]
    s = compute_stats(rv, jlv)
    println("  " * rpad(label, 40) * fmt(s))
end

println("\n── Soil temperature: surface vs. deepest (direct) ──")
jl_soilT = ustrip.(u"°C", jl_out.soil_temperature)
for (label, jl_col, r_col) in [("surface", 1, 1), ("deepest", size(jl_soilT, 2), size(r_tsoil, 2))]
    s = compute_stats(r_tsoil[:, r_col], jl_soilT[:, jl_col])
    println("  " * rpad(label, 10) * fmt(s))
end

println("\n── Soil moisture: surface vs. deepest (direct) ──")
jl_soilM = jl_out.soil_moisture   # fraction, m^3/m^3
for (label, jl_col, r_col) in [("surface", 1, 1), ("deepest", size(jl_soilM, 2), size(r_theta, 2))]
    s = compute_stats(r_theta[:, r_col], jl_soilM[:, jl_col])
    println("  " * rpad(label, 10) * fmt(s))
end

# ── Height-profile plot at a representative daytime hour (peak swdown) ───────
hour_idx = argmax(climdata.swdown)
heights_m = ustrip.(u"m", canopy_heights)   # 1.25, 2.5, ..., 25.0 -- shared axis, both flipped to top-to-bottom above
heights_desc = reverse(heights_m)           # top-to-bottom, matching column 1 = top after flip()

jl_col_downSW = ustrip.(u"W/m^2", jl_out.canopy.boundary_downward_shortwave[hour_idx, 1:n_layers])
jl_col_upSW   = ustrip.(u"W/m^2", jl_out.canopy.boundary_upward_shortwave[hour_idx, 1:n_layers])
jl_col_downLW = ustrip.(u"W/m^2", jl_out.canopy.boundary_downward_longwave[hour_idx, 1:n_layers])
jl_col_upLW   = ustrip.(u"W/m^2", jl_out.canopy.boundary_upward_longwave[hour_idx, 1:n_layers])
r_downSW = r_Rdirdown[hour_idx, :] .+ r_Rdifdown[hour_idx, :]
r_upSW   = r_Rswup[hour_idx, :]
r_downLW = r_Rlwdown[hour_idx, :]
r_upLW   = r_Rlwup[hour_idx, :]

fig = plot(layout = (2, 3), size = (1300, 700), dpi = 120,
           left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
           plot_title = "Profile at hour $(hour_idx) ($(climdata.obs_time[hour_idx])) — " *
               "Microclimate.jl ($(round(julia_elapsed, digits=1))s) vs micropoint ($(round(timing_r.elapsed_s, digits=1))s)")

plot!(fig[1], jl_air[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
      title = "Air temp", xlabel = "°C", ylabel = "Height (m)")
plot!(fig[1], r_tair[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

plot!(fig[2], jl_leaf[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
      title = "Leaf temp", xlabel = "°C")
plot!(fig[2], r_tleaf[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

plot!(fig[3], jl_rh[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
      title = "Relative humidity", xlabel = "%")
plot!(fig[3], r_relhum[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

plot!(fig[4], jl_col_downSW, heights_desc; label = "Julia down", color = :tomato, marker = :circle,
      title = "Shortwave", xlabel = "W/m²", ylabel = "Height (m)")
plot!(fig[4], r_downSW, heights_desc; label = "micropoint down", color = :steelblue, marker = :circle)
plot!(fig[4], jl_col_upSW, heights_desc; label = "Julia up", color = :tomato, ls = :dash, marker = :diamond)
plot!(fig[4], r_upSW, heights_desc; label = "micropoint up", color = :steelblue, ls = :dash, marker = :diamond)

plot!(fig[5], jl_col_downLW, heights_desc; label = "Julia down", color = :tomato, marker = :circle,
      title = "Longwave", xlabel = "W/m²")
plot!(fig[5], r_downLW, heights_desc; label = "micropoint down", color = :steelblue, marker = :circle)
plot!(fig[5], jl_col_upLW, heights_desc; label = "Julia up", color = :tomato, ls = :dash, marker = :diamond)
plot!(fig[5], r_upLW, heights_desc; label = "micropoint up", color = :steelblue, ls = :dash, marker = :diamond)

plot!(fig[6], jl_wind[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
      title = "Wind speed", xlabel = "m/s")
plot!(fig[6], r_windspeed[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

mkpath(joinpath(@__DIR__, "outputs"))
savefig(fig, joinpath(outdir, "vignette_profile_comparison.png"))
display(fig)
println("\nProfile plot saved to $(joinpath(outdir, "vignette_profile_comparison.png"))")

# ── Time-series plot: canopy top air/leaf temp over a representative week ────
week_range = 1:175200#hour_idx-84:hour_idx+84  #1:175200  # ~1 week centred on the profile hour
week_range = clamp.(week_range, 1, nhours)

fig2 = plot(size = (1100, 400), dpi = 120, title = "Canopy-top air/leaf temperature, centred on hour $hour_idx",
            xlabel = "Hour", ylabel = "°C")
plot!(fig2, week_range, jl_air[week_range, 1]; label = "Julia air (top)", color = :tomato)
plot!(fig2, week_range, r_tair[week_range, 1]; label = "micropoint air (top)", color = :steelblue)
plot!(fig2, week_range, jl_leaf[week_range, 1]; label = "Julia leaf (top)", color = :tomato, ls = :dash)
plot!(fig2, week_range, r_tleaf[week_range, 1]; label = "micropoint leaf (top)", color = :steelblue, ls = :dash)
savefig(fig2, joinpath(outdir, "vignette_timeseries_comparison.png"))
display(fig2)
println("Time-series plot saved to $(joinpath(outdir, "vignette_timeseries_comparison.png"))")

# ── Diagnostic: which hours drive the wind-speed stats blowup? ───────────────
# r≈-0.03 with RMSE in the hundreds/thousands (m/s!) at both top AND bottom
# points at a shared handful of outlier hours, not a systematic unit error —
# print the worst offenders (date + both raw values) to characterize it.
println("\n── Wind-speed outlier diagnostic (top layer, largest |Julia - micropoint|) ──")
diffs = abs.(jl_wind[:, 1] .- r_windspeed[:, 1])
worst = sortperm(diffs; rev = true)[1:10]
for i in worst
    @printf("  hour %5d (%s): Julia=%10.3f m/s  micropoint=%10.3f m/s  |diff|=%10.3f\n",
        i, climdata.obs_time[i], jl_wind[i, 1], r_windspeed[i, 1], diffs[i])
end

# ── Time-series diagnostic: wind speed / humidity / radiation at top+bottom ──
fig3 = plot(layout = (2, 2), size = (1300, 700), dpi = 120,
            left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
            plot_title = "Diagnostic time series, centred on hour $hour_idx")

plot!(fig3[1], week_range, jl_wind[week_range, 1]; label = "Julia (top)", color = :tomato,
      title = "Wind speed", ylabel = "m/s")
plot!(fig3[1], week_range, r_windspeed[week_range, 1]; label = "micropoint (top)", color = :steelblue)
plot!(fig3[1], week_range, jl_wind[week_range, n_layers]; label = "Julia (bottom)", color = :tomato, ls = :dash)
plot!(fig3[1], week_range, r_windspeed[week_range, n_layers]; label = "micropoint (bottom)", color = :steelblue, ls = :dash)

plot!(fig3[2], week_range, jl_rh[week_range, 1]; label = "Julia (top)", color = :tomato,
      title = "Relative humidity", ylabel = "%")
plot!(fig3[2], week_range, r_relhum[week_range, 1]; label = "micropoint (top)", color = :steelblue)
plot!(fig3[2], week_range, jl_rh[week_range, n_layers]; label = "Julia (bottom)", color = :tomato, ls = :dash)
plot!(fig3[2], week_range, r_relhum[week_range, n_layers]; label = "micropoint (bottom)", color = :steelblue, ls = :dash)

jl_downSW_top = ustrip.(u"W/m^2", jl_out.canopy.boundary_downward_shortwave[week_range, 1])
r_downSW_top  = (r_Rdirdown[week_range, 1] .+ r_Rdifdown[week_range, 1])
plot!(fig3[3], week_range, jl_downSW_top; label = "Julia", color = :tomato,
      title = "Downward shortwave (canopy top)", ylabel = "W/m²")
plot!(fig3[3], week_range, r_downSW_top; label = "micropoint", color = :steelblue)

jl_downLW_top = ustrip.(u"W/m^2", jl_out.canopy.boundary_downward_longwave[week_range, 1])
plot!(fig3[4], week_range, jl_downLW_top; label = "Julia", color = :tomato,
      title = "Downward longwave (canopy top)", ylabel = "W/m²", xlabel = "Hour")
plot!(fig3[4], week_range, r_Rlwdown[week_range, 1]; label = "micropoint", color = :steelblue)

savefig(fig3, joinpath(outdir, "vignette_diagnostic_timeseries.png"))
display(fig3)
println("Diagnostic time-series plot saved to $(joinpath(outdir, "vignette_diagnostic_timeseries.png"))")

# ── Time series: soil temperature + soil moisture, surface & deepest, full year ─
fig4 = plot(layout = (2, 1), size = (1200, 700), dpi = 120,
            left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
            plot_title = "Soil temperature and moisture, full year")

plot!(fig4[1], 1:nhours, jl_soilT[:, 1]; label = "Julia (surface)", color = :tomato,
      title = "Soil temperature", ylabel = "°C")
plot!(fig4[1], 1:nhours, r_tsoil[:, 1]; label = "micropoint (surface)", color = :steelblue)
plot!(fig4[1], 1:nhours, jl_soilT[:, end]; label = "Julia (deepest)", color = :tomato, ls = :dash)
plot!(fig4[1], 1:nhours, r_tsoil[:, end]; label = "micropoint (deepest)", color = :steelblue, ls = :dash)

plot!(fig4[2], 1:nhours, jl_soilM[:, 1]; label = "Julia (surface)", color = :tomato,
      title = "Soil moisture", ylabel = "m³/m³", xlabel = "Hour")
plot!(fig4[2], 1:nhours, r_theta[:, 1]; label = "micropoint (surface)", color = :steelblue)
plot!(fig4[2], 1:nhours, jl_soilM[:, end]; label = "Julia (deepest)", color = :tomato, ls = :dash)
plot!(fig4[2], 1:nhours, r_theta[:, end]; label = "micropoint (deepest)", color = :steelblue, ls = :dash)

savefig(fig4, joinpath(outdir, "vignette_soil_timeseries.png"))
display(fig4)
println("Soil time-series plot saved to $(joinpath(outdir, "vignette_soil_timeseries.png"))")
