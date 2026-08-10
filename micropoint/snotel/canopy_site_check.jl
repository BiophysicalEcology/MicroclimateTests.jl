# canopy_site_check.jl — self-contained sanity check of MultilayerCanopy leaf
# vs. air temperature, driven by real GRIDMET forcing at one SNOTEL site via
# the already-working comparisons/scan_snotel pipeline (same machinery
# comparison.jl uses for the bare-ground micropoint comparison). Also runs
# micropoint (bare ground -- it has no vegetated mode with a Microclimate.jl
# equivalent, see micropoint/README.md) on the same above-canopy forcing, so
# the soil-temperature side gets the same cross-implementation check
# comparison.jl already does, as a sanity check that the canopy addition
# doesn't corrupt the (independently validated) soil physics.
#
# Point of this script: the micropoint/vignette/ comparison shows a
# persistent leaf-warmer-than-air offset (even at night) that doesn't
# reproduce in isolated single-hour tests; this checks whether it reproduces
# under a real, independently-validated weather pipeline instead of the
# vignette's bespoke R-forcing translation -- if not, the bug is likely
# specific to that translation, not the canopy energy balance itself.
#
# Run with:
#   cd c:/git/BiophysicalEcologyEnv
#   julia --project=. c:/git/MicroclimateTests.jl/micropoint/canopy_site_check.jl

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots, Serialization
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using GeoInterface: Wrappers as GIW

gr()

const SNOTEL_DIR = joinpath(@__DIR__, "..", "comparisons", "scan_snotel")
include(joinpath(SNOTEL_DIR, "utils.jl"))
include(joinpath(SNOTEL_DIR, "config.jl"))
include(joinpath(SNOTEL_DIR, "pipeline.jl"))
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "write_micropoint_inputs.jl"))

outdir = joinpath(@__DIR__, "outputs")
mkpath(outdir)

# ── Canopy-enabled MicroModel (mirrors build_micro_model() from pipeline.jl,
# adding canopy_model and swapping in canopy-spanning heights -- 10 canopy
# layers, matching what the user's own ad-hoc test needed) ──────────────────
canopy_height = 1.0u"m"
canopy_heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"  # 10 canopy layers + reference height above canopy
n_layers = length(canopy_heights) - 1
heights_desc = reverse(ustrip.(u"m", canopy_heights[1:n_layers]))  # top(1)-to-bottom(n_layers), matching output column order

canopy_model = MultilayerCanopy(;
    # plant_area_index = 3.0 over a 1m canopy (3.0 PAI/m) made micropoint's
    # vegetated RunModelFull NaN out under calm/stable nights -- ~19x denser
    # per metre than its own forest example (PAI=4 over 25m); 1.0 keeps a
    # modest, still-plausible low canopy without that numerical fragility.
    canopy_height, plant_area_index = 1.0,
    leaf_temperature_solver = RootFindLeafTemperature(),
    convergence = IterationToleranceConvergence(; tolerance = 0.01u"K", max_iterations_per_day = 30),
)

soil_properties_model = CampbelldeVriesSoilProperties(;
    de_vries_shape_factor, recirculation_power, return_flow_threshold,
)
soil_profile = SoilProfile(;
    bulk_density    = fill(bulk_density, length(depths)),
    mineral_density = fill(mineral_density, length(depths)),
    mineral_conductivity, mineral_heat_capacity,
    hydraulics = CampbellHydraulicProfile(;
        air_entry_water_potential        = fill(air_entry_potential, length(depths)),
        saturated_hydraulic_conductivity = fill(sat_hydraulic_cond, length(depths)),
        campbell_b_parameter             = fill(campbell_b, length(depths)),
        root_density,
    ),
)
soil_hydraulic_model = CampbellSoilHydraulics(;
    root_resistance, stomatal_closure_potential, leaf_resistance,
    stomatal_stability_parameter, root_radius,
    infiltration_algorithm = infiltration_algorithm_choice,
)
snow_model = SnowModel(;
    snow_temperature_threshold = snow_temp_threshold, snow_density = snow_density_init,
    snow_melt_factor, undercatch, rain_multiplier, rain_melt_factor,
    density_function, snow_conductivity, canopy_interception,
)
micro_config = MicroConfig(;
    convergence = convergence_choice, rainfall_schedule = rainfall_schedule_choice,
    soil_moisture_strategy = soil_moisture_strategy_choice,
)
micro_model = MicroModel(;
    hours = collect(0.0:1:23.0), depths, heights = canopy_heights,
    soil_properties_model, soil_hydraulic_model, snow_model, canopy_model,
    config = micro_config,
)
# NoCanopy baseline (defaults to NoCanopy(), profile_heights matches the
# ordinary bare-ground comparison.jl setup) -- run alongside the canopy
# version so the soil comparison can isolate the canopy's own effect from any
# other model difference vs. micropoint.
micro_model_nocanopy = MicroModel(;
    hours = collect(0.0:1:23.0), depths, heights = profile_heights,
    soil_properties_model, soil_hydraulic_model, snow_model,
    config = micro_config,
)

# ── Fetch real weather + build the problem for one site (reuses the pipeline's
# own prepare_site -- same fetch/cache/obs/init logic comparison.jl already
# validates, just handed this canopy-enabled micro_model instead) ───────────
meta_all = DataFrame(CSV.File(joinpath(SNOTEL_DIR, "Map metadata export.csv")))
meta_all[!, :ID] = strip.(string.(meta_all[!, :ID]))
weather_cache = Dict{Tuple{Int,Date,Date}, Any}()

prep = prepare_site(site_num, meta_all, weather_cache, micro_model, soil_profile;
    sim_start = run_start, sim_end, auto_date_range = false, max_sim_years = spinup_years + 1)
isnothing(prep) && error("Site $site_num failed to prepare -- see warnings above.")

println("\nSolving Microclimate.jl (MultilayerCanopy, $(length(canopy_heights)-1) layers, site $site_num)...")
julia_elapsed = @elapsed (out_canopy = Microclimate.solve(prep.problem))
@printf("Solved in %.2f s (%.4f s/day)\n", julia_elapsed, julia_elapsed / prep.ndays)

# NoCanopy baseline, same site/weather/soil -- isolates the canopy's own
# effect on soil temperature from any other Julia-vs-micropoint difference.
prep_nocanopy = prepare_site(site_num, meta_all, weather_cache, micro_model_nocanopy, soil_profile;
    sim_start = run_start, sim_end, auto_date_range = false, max_sim_years = spinup_years + 1)
isnothing(prep_nocanopy) && error("Site $site_num failed to prepare (NoCanopy) -- see warnings above.")

println("\nSolving Microclimate.jl (NoCanopy, site $site_num)...")
julia_elapsed_nocanopy = @elapsed (out_nocanopy = Microclimate.solve(prep_nocanopy.problem))
@printf("Solved in %.2f s (%.4f s/day)\n", julia_elapsed_nocanopy, julia_elapsed_nocanopy / prep_nocanopy.ndays)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  micropoint (bare ground) on the same above-canopy forcing --            ║
# ║  soil-temperature sanity check (adapted from comparison.jl)              ║
# ╚══════════════════════════════════════════════════════════════════════════╝

mp_outdir = mcf_outdir(site_num)
write_micropoint_inputs(out_canopy, prep, mp_outdir)

# Canopy spec for run_micropoint_vegetated.R -- only the fields that map onto
# createvegp()'s own parameters (see that script for the full mapping/caveats).
CSV.write(joinpath(mp_outdir, "canopy_params.csv"), DataFrame(
    canopy_height_m = [ustrip(u"m", canopy_model.canopy_height)],
    pai = [canopy_model.plant_area_index],
    n_layers = [n_layers],
    leaf_length_m = [ustrip(u"m", canopy_model.leaf_parameters.leaf_length)],
    leaf_width_m = [ustrip(u"m", canopy_model.leaf_parameters.leaf_width)],
    leaf_emissivity = [canopy_model.leaf_parameters.leaf_emissivity],
    leaf_reflectance = [canopy_model.shortwave_model.leaf_reflectance],
    leaf_transmittance = [canopy_model.shortwave_model.leaf_transmittance],
    leaf_angle_x = [canopy_model.leaf_parameters.canopy_projection_ratio],
))

mp_rscript = joinpath(@__DIR__, "run_micropoint.R")
println("\nRunning micropoint (bare ground)...")
mp_elapsed = @elapsed run(`Rscript $mp_rscript $mp_outdir`)
@printf("micropoint solved in %.2f s\n", mp_elapsed)

veg_rscript = joinpath(@__DIR__, "run_micropoint_vegetated.R")
println("\nRunning micropoint (vegetated, $n_layers layers, matching the Julia canopy spec)...")
veg_elapsed = @elapsed run(`Rscript $veg_rscript $mp_outdir`)
@printf("micropoint (vegetated) solved in %.2f s\n", veg_elapsed)

mcf_out = DataFrame(CSV.File(joinpath(mp_outdir, "micropoint_out.csv"); missingstring = "NA"))
mcf_out.obs_time = DateTime.(mcf_out.obs_time, dateformat"yyyy-mm-dd HH:MM:SS")
mcf_out.obs_time = mcf_out.obs_time .- Hour(round(Int, -prep.lon_dd / 15))
for col in names(mcf_out)
    col == "obs_time" && continue
    mcf_out[!, col] = Float64.(coalesce.(mcf_out[!, col], NaN))
end

depths_m = ustrip.(u"m", depths)
mc_idx(target_cm) = argmin(abs.(depths_m .* 100 .- target_cm))

t_model_full = [DateTime(d) + Hour(h) for d in prep.dates_vec for h in 0:23]
eval_idx_mc = findall(t -> DateTime(sim_start) <= t <= DateTime(sim_end, Time(23,0,0)), t_model_full)
eval_idx_r  = findall(t -> DateTime(sim_start) <= t <= DateTime(sim_end, Time(23,0,0)), mcf_out.obs_time)
n_common = min(length(eval_idx_mc), length(eval_idx_r))
eval_idx_mc, eval_idx_r = eval_idx_mc[1:n_common], eval_idx_r[1:n_common]

jl_soilT = ustrip.(u"°C", out_canopy.soil_temperature[eval_idx_mc, :])
jl_soilT_nocanopy = ustrip.(u"°C", out_nocanopy.soil_temperature[eval_idx_mc, :])
r_soilT  = mcf_out[eval_idx_r, :]

soil_depth_cols = [("5cm", mc_idx(5), :T_5cm), ("20cm", mc_idx(20), :T_20cm), ("50cm", mc_idx(50), :T_50cm)]

println("\n── Soil temperature: Microclimate.jl (NoCanopy) vs micropoint (bare ground), site $site_num ──")
println("  (cleanest apples-to-apples baseline -- both bare ground)")
for (label, jcol, rcol) in soil_depth_cols
    s = compute_stats(r_soilT[!, rcol], jl_soilT_nocanopy[:, jcol])
    println("  " * rpad(label, 8) * fmt_stat(s))
end

println("\n── Soil temperature: Microclimate.jl (canopy) vs micropoint (bare ground), site $site_num ──")
println("  (canopy vs. bare ground -- expect a systematic offset, not just noise)")
for (label, jcol, rcol) in soil_depth_cols
    s = compute_stats(r_soilT[!, rcol], jl_soilT[:, jcol])
    println("  " * rpad(label, 8) * fmt_stat(s))
end

println("\n── Soil temperature: Microclimate.jl (canopy) vs Microclimate.jl (NoCanopy), site $site_num ──")
println("  (isolates the canopy's own effect, no micropoint involved)")
for (label, jcol, rcol) in soil_depth_cols
    s = compute_stats(jl_soilT_nocanopy[:, jcol], jl_soilT[:, jcol])
    println("  " * rpad(label, 8) * fmt_stat(s))
end

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Canopy leaf/air temperature diagnostic (the actual point of this        ║
# ║  script) -- restricted to the reporting window (spin-up discarded)       ║
# ╚══════════════════════════════════════════════════════════════════════════╝

t_model = t_model_full[eval_idx_mc]

# Date window shared by every time-series plot below (profile plots are keyed
# to a single hour instead, unaffected) -- reuses plot_start/plot_end from
# micropoint/config.jl (defaults to the full reporting year); narrow those to
# zoom in, e.g. `plot_start = Date(2015, 6, 1); plot_end = Date(2015, 6, 7)`.
plot_start = Date(2015, 1, 1)
plot_end   = Date(2015, 3, 31)
plot_idx = findall(t -> DateTime(plot_start) <= t <= DateTime(plot_end, Time(23, 0, 0)), t_model)

leaf_top = ustrip.(u"°C", out_canopy.canopy.leaf_temperature[eval_idx_mc, 1])
air_top  = ustrip.(u"°C", out_canopy.canopy.air_temperature[eval_idx_mc, 1])
leaf_bot = ustrip.(u"°C", out_canopy.canopy.leaf_temperature[eval_idx_mc, end])
air_bot  = ustrip.(u"°C", out_canopy.canopy.air_temperature[eval_idx_mc, end])
gap_top  = leaf_top .- air_top
gap_bot  = leaf_bot .- air_bot

println("\n── Leaf - air gap (site $site_num, $sim_start to $sim_end) ──")
@printf("  Top layer:    mean=%.3f K  max=%.3f K  min=%.3f K\n", mean(gap_top), maximum(gap_top), minimum(gap_top))
@printf("  Bottom layer: mean=%.3f K  max=%.3f K  min=%.3f K\n", mean(gap_bot), maximum(gap_bot), minimum(gap_bot))

sw = ustrip.(u"W/m^2", out_canopy.global_radiation[eval_idx_mc])
night = sw .<= 0.0
@printf("  Top layer, night hours only:    mean=%.3f K  max=%.3f K  (n=%d)\n",
    mean(gap_top[night]), maximum(abs.(gap_top[night])), count(night))
@printf("  Bottom layer, night hours only: mean=%.3f K  max=%.3f K  (n=%d)\n",
    mean(gap_bot[night]), maximum(abs.(gap_bot[night])), count(night))

fig = plot(size = (1100, 400), dpi = 120, title = "Canopy leaf/air temperature, site $site_num, $plot_start to $plot_end",
    xlabel = "Date", ylabel = "°C")
plot!(fig, t_model[plot_idx], air_top[plot_idx], label = "air (top)", color = :steelblue)
plot!(fig, t_model[plot_idx], leaf_top[plot_idx], label = "leaf (top)", color = :steelblue, ls = :dash)
plot!(fig, t_model[plot_idx], air_bot[plot_idx], label = "air (bottom)", color = :tomato)
plot!(fig, t_model[plot_idx], leaf_bot[plot_idx], label = "leaf (bottom)", color = :tomato, ls = :dash)
savefig(fig, joinpath(outdir, "canopy_site_check_$(site_num).png"))
display(fig)
println("\nPlot saved to $(joinpath(outdir, "canopy_site_check_$(site_num).png"))")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  micropoint (vegetated) vs. Microclimate.jl -- same comparison set as     ║
# ║  micropoint/vignette/vignette_comparison.jl, but on real GRIDMET forcing  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

read_veg(name) = Matrix{Float64}(coalesce.(Matrix(DataFrame(
    CSV.File(joinpath(mp_outdir, "veg_" * name * ".csv"); missingstring = "NA"))), NaN))
r_tair_full, r_tleaf_full, r_relhum_full, r_windspeed_full =
    read_veg("tair"), read_veg("tleaf"), read_veg("relhum"), read_veg("windspeed")
r_Rdirdown_full, r_Rdifdown_full, r_Rswup_full, r_Rlwdown_full, r_Rlwup_full =
    read_veg("Rdirdown"), read_veg("Rdifdown"), read_veg("Rswup"), read_veg("Rlwdown"), read_veg("Rlwup")

# micropoint's canopy-layer columns run bottom(1)-to-top(n) -- flip once so
# every column index lines up with Microclimate.jl's top(1)-to-bottom(n).
flip(m) = m[:, end:-1:1]
r_tair_full, r_tleaf_full, r_relhum_full, r_windspeed_full =
    flip(r_tair_full), flip(r_tleaf_full), flip(r_relhum_full), flip(r_windspeed_full)
r_Rdirdown_full, r_Rdifdown_full, r_Rswup_full, r_Rlwdown_full, r_Rlwup_full =
    flip(r_Rdirdown_full), flip(r_Rdifdown_full), flip(r_Rswup_full), flip(r_Rlwdown_full), flip(r_Rlwup_full)

# veg_*.csv rows are climdata.csv's own row order (all of t_model_full, no
# separate obs_time to align) -- eval_idx_mc indexes directly.
r_tair, r_tleaf, r_relhum, r_windspeed =
    r_tair_full[eval_idx_mc, :], r_tleaf_full[eval_idx_mc, :], r_relhum_full[eval_idx_mc, :], r_windspeed_full[eval_idx_mc, :]
r_downSW = r_Rdirdown_full[eval_idx_mc, :] .+ r_Rdifdown_full[eval_idx_mc, :]
r_upSW   = r_Rswup_full[eval_idx_mc, :]
r_downLW = r_Rlwdown_full[eval_idx_mc, :]
r_upLW   = r_Rlwup_full[eval_idx_mc, :]

jl_air_h  = ustrip.(u"°C", out_canopy.canopy.air_temperature[eval_idx_mc, :])
jl_leaf_h = ustrip.(u"°C", out_canopy.canopy.leaf_temperature[eval_idx_mc, :])
jl_rh_h   = out_canopy.canopy.relative_humidity[eval_idx_mc, :] .* 100.0
jl_wind_h = ustrip.(u"m/s", out_canopy.canopy.wind_speed[eval_idx_mc, :])
jl_downSW = ustrip.(u"W/m^2", out_canopy.canopy.boundary_downward_shortwave[eval_idx_mc, 1:n_layers])
jl_upSW   = ustrip.(u"W/m^2", out_canopy.canopy.boundary_upward_shortwave[eval_idx_mc, 1:n_layers])
jl_downLW = ustrip.(u"W/m^2", out_canopy.canopy.boundary_downward_longwave[eval_idx_mc, 1:n_layers])
jl_upLW   = ustrip.(u"W/m^2", out_canopy.canopy.boundary_upward_longwave[eval_idx_mc, 1:n_layers])

println("\n── Canopy top vs. bottom (Microclimate.jl vs. micropoint vegetated, hourly, direct) ──")
for (label, col) in [("top", 1), ("bottom", n_layers)]
    println("  Layer: $label")
    for (name, jlm, rm) in [
        ("Air temp (°C)",  jl_air_h,  r_tair),
        ("Leaf temp (°C)", jl_leaf_h, r_tleaf),
        ("Rel. humidity (%)", jl_rh_h, r_relhum),
        ("Wind speed (m/s)", jl_wind_h, r_windspeed),
    ]
        s = compute_stats(rm[:, col], jlm[:, col])
        println("    " * rpad(name, 20) * fmt_stat(s))
    end
end

# ── Height-profile plot at a representative daytime hour (peak swdown) ───────
hour_idx = argmax(sw)

fig_profile = plot(layout = (2, 3), size = (1300, 700), dpi = 120,
    left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
    plot_title = "Profile at hour $(hour_idx) ($(t_model[hour_idx])) — site $site_num, " *
        "Microclimate.jl ($(round(julia_elapsed,digits=1))s) vs micropoint ($(round(veg_elapsed,digits=1))s)")

plot!(fig_profile[1], jl_air_h[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
    title = "Air temp", xlabel = "°C", ylabel = "Height (m)")
plot!(fig_profile[1], r_tair[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

plot!(fig_profile[2], jl_leaf_h[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
    title = "Leaf temp", xlabel = "°C")
plot!(fig_profile[2], r_tleaf[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

plot!(fig_profile[3], jl_rh_h[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
    title = "Relative humidity", xlabel = "%")
plot!(fig_profile[3], r_relhum[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

plot!(fig_profile[4], jl_downSW[hour_idx, :], heights_desc; label = "Julia down", color = :tomato, marker = :circle,
    title = "Shortwave", xlabel = "W/m²", ylabel = "Height (m)")
plot!(fig_profile[4], r_downSW[hour_idx, :], heights_desc; label = "micropoint down", color = :steelblue, marker = :circle)
plot!(fig_profile[4], jl_upSW[hour_idx, :], heights_desc; label = "Julia up", color = :tomato, ls = :dash, marker = :diamond)
plot!(fig_profile[4], r_upSW[hour_idx, :], heights_desc; label = "micropoint up", color = :steelblue, ls = :dash, marker = :diamond)

plot!(fig_profile[5], jl_downLW[hour_idx, :], heights_desc; label = "Julia down", color = :tomato, marker = :circle,
    title = "Longwave", xlabel = "W/m²")
plot!(fig_profile[5], r_downLW[hour_idx, :], heights_desc; label = "micropoint down", color = :steelblue, marker = :circle)
plot!(fig_profile[5], jl_upLW[hour_idx, :], heights_desc; label = "Julia up", color = :tomato, ls = :dash, marker = :diamond)
plot!(fig_profile[5], r_upLW[hour_idx, :], heights_desc; label = "micropoint up", color = :steelblue, ls = :dash, marker = :diamond)

plot!(fig_profile[6], jl_wind_h[hour_idx, :], heights_desc; label = "Julia", color = :tomato, marker = :circle,
    title = "Wind speed", xlabel = "m/s")
plot!(fig_profile[6], r_windspeed[hour_idx, :], heights_desc; label = "micropoint", color = :steelblue, marker = :circle)

savefig(fig_profile, joinpath(outdir, "canopy_site_check_profile_$(site_num).png"))
display(fig_profile)
println("Profile plot saved to $(joinpath(outdir, "canopy_site_check_profile_$(site_num).png"))")

# ── Time-series plot: canopy-top air/leaf temp over a representative week ────
week_range = clamp.(hour_idx-84:hour_idx+84, 1, length(t_model))

fig_ts = plot(size = (1100, 400), dpi = 120, title = "Canopy-top air/leaf temperature, centred on hour $hour_idx",
    xlabel = "Hour", ylabel = "°C")
plot!(fig_ts, week_range, jl_air_h[week_range, 1]; label = "Julia air (top)", color = :tomato)
plot!(fig_ts, week_range, r_tair[week_range, 1]; label = "micropoint air (top)", color = :steelblue)
plot!(fig_ts, week_range, jl_leaf_h[week_range, 1]; label = "Julia leaf (top)", color = :tomato, ls = :dash)
plot!(fig_ts, week_range, r_tleaf[week_range, 1]; label = "micropoint leaf (top)", color = :steelblue, ls = :dash)
savefig(fig_ts, joinpath(outdir, "canopy_site_check_timeseries_$(site_num).png"))
display(fig_ts)
println("Time-series plot saved to $(joinpath(outdir, "canopy_site_check_timeseries_$(site_num).png"))")

# ── Time series diagnostic: wind speed / humidity / radiation at top+bottom ──
fig_diag = plot(layout = (2, 2), size = (1300, 700), dpi = 120,
    left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
    plot_title = "Diagnostic time series, centred on hour $hour_idx")

plot!(fig_diag[1], week_range, jl_wind_h[week_range, 1]; label = "Julia (top)", color = :tomato,
    title = "Wind speed", ylabel = "m/s")
plot!(fig_diag[1], week_range, r_windspeed[week_range, 1]; label = "micropoint (top)", color = :steelblue)
plot!(fig_diag[1], week_range, jl_wind_h[week_range, n_layers]; label = "Julia (bottom)", color = :tomato, ls = :dash)
plot!(fig_diag[1], week_range, r_windspeed[week_range, n_layers]; label = "micropoint (bottom)", color = :steelblue, ls = :dash)

plot!(fig_diag[2], week_range, jl_rh_h[week_range, 1]; label = "Julia (top)", color = :tomato,
    title = "Relative humidity", ylabel = "%")
plot!(fig_diag[2], week_range, r_relhum[week_range, 1]; label = "micropoint (top)", color = :steelblue)
plot!(fig_diag[2], week_range, jl_rh_h[week_range, n_layers]; label = "Julia (bottom)", color = :tomato, ls = :dash)
plot!(fig_diag[2], week_range, r_relhum[week_range, n_layers]; label = "micropoint (bottom)", color = :steelblue, ls = :dash)

plot!(fig_diag[3], week_range, jl_downSW[week_range, 1]; label = "Julia", color = :tomato,
    title = "Downward shortwave (canopy top)", ylabel = "W/m²")
plot!(fig_diag[3], week_range, r_downSW[week_range, 1]; label = "micropoint", color = :steelblue)

plot!(fig_diag[4], week_range, jl_downLW[week_range, 1]; label = "Julia", color = :tomato,
    title = "Downward longwave (canopy top)", ylabel = "W/m²", xlabel = "Hour")
plot!(fig_diag[4], week_range, r_downLW[week_range, 1]; label = "micropoint", color = :steelblue)

savefig(fig_diag, joinpath(outdir, "canopy_site_check_diagnostic_$(site_num).png"))
display(fig_diag)
println("Diagnostic time-series plot saved to $(joinpath(outdir, "canopy_site_check_diagnostic_$(site_num).png"))")

# ── Time series: soil temperature + soil moisture vs micropoint (bare ground)
# -- one figure with canopy, one without, so each is a clean two-model
# comparison rather than three lines competing on one axis ──────────────────
jl_soilT5, jl_soilT50 = jl_soilT[:, mc_idx(5)], jl_soilT[:, mc_idx(50)]
jl_soilT5_nc, jl_soilT50_nc = jl_soilT_nocanopy[:, mc_idx(5)], jl_soilT_nocanopy[:, mc_idx(50)]
jl_soilM = out_canopy.soil_moisture[eval_idx_mc, :]
jl_soilM5, jl_soilM50 = jl_soilM[:, mc_idx(5)], jl_soilM[:, mc_idx(50)]
jl_soilM_nc = out_nocanopy.soil_moisture[eval_idx_mc, :]
jl_soilM5_nc, jl_soilM50_nc = jl_soilM_nc[:, mc_idx(5)], jl_soilM_nc[:, mc_idx(50)]
r_soilM = mcf_out[eval_idx_r, :]

# prep.obs (from prepare_site) is the same real SCAN/SNOTEL station regardless
# of canopy -- shared by both plots below. om restricts it to the plot window;
# obs_vals guards for depths/columns some sites don't have.
om = (prep.obs.DateTime .>= DateTime(plot_start)) .& (prep.obs.DateTime .<= DateTime(plot_end, Time(23, 0, 0)))
obs_vals(col) = Symbol(col) in propertynames(prep.obs) ?
    Float64.(coalesce.(prep.obs[om, Symbol(col)], NaN)) : Float64[]

function soil_comparison_plot(label, jsoilT5, jsoilT50, jsoilM5, jsoilM50, fname)
    fig = plot(layout = (2, 1), size = (1200, 700), dpi = 120,
        left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
        plot_title = "Soil temperature and moisture ($label), $plot_start to $plot_end — site $site_num")
    plot!(fig[1], t_model[plot_idx], jsoilT5[plot_idx]; label = "Julia (5cm)", color = :tomato,
        title = "Soil temperature", ylabel = "°C")
    plot!(fig[1], t_model[plot_idx], r_soilT[plot_idx, :T_5cm]; label = "micropoint (5cm)", color = :steelblue)
    plot!(fig[1], t_model[plot_idx], jsoilT50[plot_idx]; label = "Julia (50cm)", color = :tomato, ls = :dash)
    plot!(fig[1], t_model[plot_idx], r_soilT[plot_idx, :T_50cm]; label = "micropoint (50cm)", color = :steelblue, ls = :dash)
    !isempty(obs_vals("STO_5cm")) &&
        scatter!(fig[1], prep.obs.DateTime[om], obs_vals("STO_5cm"); color = :green, ms = 1, alpha = 0.7, label = "obs (5cm)")
    !isempty(obs_vals("STO_50cm")) &&
        scatter!(fig[1], prep.obs.DateTime[om], obs_vals("STO_50cm"); color = :darkgreen, ms = 1, alpha = 0.7, label = "obs (50cm)")

    plot!(fig[2], t_model[plot_idx], jsoilM5[plot_idx]; label = "Julia (5cm)", color = :tomato,
        title = "Soil moisture", ylabel = "m³/m³", xlabel = "Date")
    plot!(fig[2], t_model[plot_idx], r_soilM[plot_idx, :soilm_bulk]; label = "micropoint (bulk 0-10cm)", color = :steelblue)
    plot!(fig[2], t_model[plot_idx], jsoilM50[plot_idx]; label = "Julia (50cm)", color = :tomato, ls = :dash)
    !isempty(obs_vals("SMS_5cm")) &&
        scatter!(fig[2], prep.obs.DateTime[om], obs_vals("SMS_5cm") ./ 100.0; color = :green, ms = 1, alpha = 0.7, label = "obs (5cm)")
    !isempty(obs_vals("SMS_50cm")) &&
        scatter!(fig[2], prep.obs.DateTime[om], obs_vals("SMS_50cm") ./ 100.0; color = :darkgreen, ms = 1, alpha = 0.7, label = "obs (50cm)")

    savefig(fig, joinpath(outdir, fname))
    display(fig)
    println("Soil time-series plot ($label) saved to $(joinpath(outdir, fname))")
end

soil_comparison_plot("canopy", jl_soilT5, jl_soilT50, jl_soilM5, jl_soilM50,
    "canopy_site_check_soil_canopy_$(site_num).png")
soil_comparison_plot("NoCanopy", jl_soilT5_nc, jl_soilT50_nc, jl_soilM5_nc, jl_soilM50_nc,
    "canopy_site_check_soil_nocanopy_$(site_num).png")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Ablation test: cloud_cover=1, blackbody longwave (sigma*T_ref^4), zero   ║
# ║  shortwave/rain, RH=100%, real air temp/wind/pressure -- same recipe the ║
# ║  user ran by hand against vignette_comparison.jl (where it triggered a   ║
# ║  DomainError before the relaxation fix). Re-run here on the real GRIDMET ║
# ║  series (via prepare_site's own already-fetched environment_hourly, not  ║
# ║  a hand-built one) to check whether that instability was specific to the ║
# ║  vignette's own forcing translation or is a real Microclimate.jl fragility║
# ╚══════════════════════════════════════════════════════════════════════════╝

println("\n" * "="^72)
println("Ablation test: cloud_cover=1, blackbody LW, zero SW/rain, RH=100%")
println("="^72)

# GRIDMET-based problems drive off environment_minmax (daily min/max) with
# hourly temperature synthesized internally during solve -- inputs.environment_hourly
# is nothing here, not a literal per-hour series. Use the already-solved
# out_canopy's own per-hour reference_temperature/pressure/reference_wind_speed
# outputs instead (real hourly values from the non-ablated run), same pattern
# write_micropoint_inputs.jl already uses to build climdata.csv for R.
old_inputs = prep.problem.inputs
nhours_full = length(out_canopy.reference_temperature)

env_h_ablation = HourlyTimeseries(;
    pressure = out_canopy.pressure[1:nhours_full],
    reference_temperature = out_canopy.reference_temperature[1:nhours_full],
    reference_humidity = fill(1.0, nhours_full),
    reference_wind_speed = out_canopy.reference_wind_speed[1:nhours_full],
    global_radiation = zeros(typeof(0.0u"W/m^2"), nhours_full),
    longwave_radiation = SIGMA_SB_MP .* ustrip.(u"K", out_canopy.reference_temperature[1:nhours_full]) .^ 4 .* 1.0u"W/m^2",
    cloud_cover = fill(1.0, nhours_full),
    rainfall = zeros(typeof(0.0u"kg/m^2"), nhours_full),
    zenith_angle = nothing,
)
inputs_ablation = MicroInputs(;
    site = old_inputs.site, soil_profile = old_inputs.soil_profile,
    # nothing, not old_inputs.environment_minmax -- interpolate_minmax! prefers
    # environment_minmax over environment_hourly for everything except
    # pressure when both are set, which would silently ignore this ablation.
    environment_minmax = nothing, environment_hourly = env_h_ablation,
    environment_daily = old_inputs.environment_daily,
    initial_soil_temperature = old_inputs.initial_soil_temperature,
    initial_soil_moisture = old_inputs.initial_soil_moisture,
    initial_snow_depth = old_inputs.initial_snow_depth,
    initial_snow_temperature = old_inputs.initial_snow_temperature,
    initial_snow_density = old_inputs.initial_snow_density,
)
problem_ablation = MicroProblem(micro_model, inputs_ablation; days = prep.problem.days, time_mode = prep.problem.time_mode)

ablation_ok = true
try
    global out_ablation
    t_abl = @elapsed (out_ablation = Microclimate.solve(problem_ablation))
    @printf("Julia (canopy) ablation solve: completed in %.2f s\n", t_abl)
catch e
    global ablation_ok = false
    println("Julia (canopy) ablation solve: FAILED -- ", sprint(showerror, e))
end

if ablation_ok
    abl_leaf_top = ustrip.(u"°C", out_ablation.canopy.leaf_temperature[eval_idx_mc, 1])
    abl_air_top  = ustrip.(u"°C", out_ablation.canopy.air_temperature[eval_idx_mc, 1])
    abl_gap_top  = abl_leaf_top .- abl_air_top
    abl_soilT    = ustrip.(u"°C", out_ablation.soil_temperature[eval_idx_mc, :])
    @printf("  Leaf - air gap (top): mean=%.3f K  max=%.3f K  min=%.3f K\n",
        mean(abl_gap_top), maximum(abl_gap_top), minimum(abl_gap_top))
    @printf("  Soil temperature: min=%.2f degC  max=%.2f degC (all depths, all hours)\n",
        minimum(abl_soilT), maximum(abl_soilT))
    @printf("  All finite: leaf=%s air=%s soilT=%s\n",
        all(isfinite, abl_leaf_top), all(isfinite, abl_air_top), all(isfinite, abl_soilT))
end

# ── Same ablation on micropoint (bare ground + vegetated), same real
# air temp/wind/pressure series, for a symmetric check ──────────────────────
abl_outdir = joinpath(mp_outdir, "ablation")
mkpath(abl_outdir)

utc_offset_hours = round(Int, -prep.lon_dd / 15)
t_model_utc = t_model_full .+ Hour(utc_offset_hours)
obs_time_str = Dates.format.(t_model_utc, dateformat"yyyy-mm-dd HH:MM:SS")

abl_climdata = DataFrame(
    obs_time  = obs_time_str,
    temp      = ustrip.(u"°C", uconvert.(u"°C", out_canopy.reference_temperature[1:nhours_full])),
    relhum    = fill(100.0, nhours_full),
    pres      = ustrip.(u"kPa", uconvert.(u"kPa", out_canopy.pressure[1:nhours_full])),
    swdown    = zeros(nhours_full),
    difrad    = zeros(nhours_full),
    lwdown    = ustrip.(u"W/m^2", env_h_ablation.longwave_radiation),
    windspeed = ustrip.(u"m/s", out_canopy.reference_wind_speed[1:nhours_full]),
    winddir   = zeros(nhours_full),
    precip    = zeros(nhours_full),
)
CSV.write(joinpath(abl_outdir, "climdata.csv"), abl_climdata)
cp(joinpath(mp_outdir, "params.csv"), joinpath(abl_outdir, "params.csv"); force = true)
cp(joinpath(mp_outdir, "canopy_params.csv"), joinpath(abl_outdir, "canopy_params.csv"); force = true)

println("\nRunning micropoint (bare ground) ablation...")
mp_abl_ok = try
    run(`Rscript $mp_rscript $abl_outdir`)
    true
catch e
    println("  FAILED -- ", sprint(showerror, e))
    false
end
if mp_abl_ok
    mp_abl = DataFrame(CSV.File(joinpath(abl_outdir, "micropoint_out.csv"); missingstring = "NA"))
    n_na = sum(ismissing, Matrix(mp_abl[!, Not(:obs_time)]))
    @printf("  %d missing/NA cells of %d\n", n_na, (ncol(mp_abl)-1) * nrow(mp_abl))
end

println("\nRunning micropoint (vegetated) ablation...")
veg_abl_ok = try
    run(`Rscript $veg_rscript $abl_outdir`)
    true
catch e
    println("  FAILED -- ", sprint(showerror, e))
    false
end
if veg_abl_ok
    veg_abl_tair = Matrix{Float64}(coalesce.(Matrix(DataFrame(
        CSV.File(joinpath(abl_outdir, "veg_tair.csv"); missingstring = "NA"))), NaN))
    @printf("  tair: %d NaN of %d\n", count(isnan, veg_abl_tair), length(veg_abl_tair))
end

# ── Ablation plots ────────────────────────────────────────────────────────
if ablation_ok
    abl_leaf_bot = ustrip.(u"°C", out_ablation.canopy.leaf_temperature[eval_idx_mc, end])
    abl_air_bot  = ustrip.(u"°C", out_ablation.canopy.air_temperature[eval_idx_mc, end])

    fig_abl = plot(size = (1100, 400), dpi = 120,
        title = "Ablation (cloud_cover=1, blackbody LW, no SW/rain): leaf/air temperature, $plot_start to $plot_end",
        xlabel = "Date", ylabel = "°C")
    plot!(fig_abl, t_model[plot_idx], abl_air_top[plot_idx]; label = "air (top)", color = :steelblue)
    plot!(fig_abl, t_model[plot_idx], abl_leaf_top[plot_idx]; label = "leaf (top)", color = :steelblue, ls = :dash)
    plot!(fig_abl, t_model[plot_idx], abl_air_bot[plot_idx]; label = "air (bottom)", color = :tomato)
    plot!(fig_abl, t_model[plot_idx], abl_leaf_bot[plot_idx]; label = "leaf (bottom)", color = :tomato, ls = :dash)
    savefig(fig_abl, joinpath(outdir, "canopy_site_check_ablation_leafair_$(site_num).png"))
    display(fig_abl)
    println("\nAblation leaf/air plot saved to $(joinpath(outdir, "canopy_site_check_ablation_leafair_$(site_num).png"))")

    if mp_abl_ok
        mp_abl.obs_time = DateTime.(mp_abl.obs_time, dateformat"yyyy-mm-dd HH:MM:SS")
        mp_abl.obs_time = mp_abl.obs_time .- Hour(round(Int, -prep.lon_dd / 15))
        for col in names(mp_abl)
            col == "obs_time" && continue
            mp_abl[!, col] = Float64.(coalesce.(mp_abl[!, col], NaN))
        end
        eval_idx_mp_abl = findall(t -> DateTime(sim_start) <= t <= DateTime(sim_end, Time(23, 0, 0)), mp_abl.obs_time)
        n_common_abl = min(length(eval_idx_mc), length(eval_idx_mp_abl))
        idx_j_abl, idx_r_abl = eval_idx_mc[1:n_common_abl], eval_idx_mp_abl[1:n_common_abl]
        r_soilT_abl = mp_abl[idx_r_abl, :]

        fig_abl_soil = plot(size = (1100, 400), dpi = 120,
            title = "Ablation: soil temperature, Julia (canopy) vs micropoint (bare ground), $plot_start to $plot_end",
            xlabel = "Date", ylabel = "°C")
        t_abl_soil = t_model_full[idx_j_abl]
        plot!(fig_abl_soil, t_abl_soil, ustrip.(u"°C", out_ablation.soil_temperature[idx_j_abl, mc_idx(5)]);
            label = "Julia (5cm)", color = :tomato)
        plot!(fig_abl_soil, t_abl_soil, r_soilT_abl[!, :T_5cm]; label = "micropoint (5cm)", color = :steelblue)
        plot!(fig_abl_soil, t_abl_soil, ustrip.(u"°C", out_ablation.soil_temperature[idx_j_abl, mc_idx(50)]);
            label = "Julia (50cm)", color = :tomato, ls = :dash)
        plot!(fig_abl_soil, t_abl_soil, r_soilT_abl[!, :T_50cm]; label = "micropoint (50cm)", color = :steelblue, ls = :dash)
        savefig(fig_abl_soil, joinpath(outdir, "canopy_site_check_ablation_soil_$(site_num).png"))
        display(fig_abl_soil)
        println("Ablation soil plot saved to $(joinpath(outdir, "canopy_site_check_ablation_soil_$(site_num).png"))")
    end
end
