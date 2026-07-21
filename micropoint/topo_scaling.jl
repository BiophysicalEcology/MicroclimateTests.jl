# topo_scaling.jl — Julia side (STAGE 1: SRTM+CRUCL2 terrain-corrected run
# and its own snapshot plot, confirmed working) plus STAGE 2: exporting a
# coarse, pre-correction CRUCL2 forcing grid + elevation reference for
# microclimf's own topographic correction (run_microclimf_topo.R) to consume
# independently. Both sides then correct the SAME coarse climate onto their
# own fine DEM -- microclimf via its `altcorrect` elevation lapse-rate
# correction (fixed or humidity-dependent, per rpubs.com/ilyamaclean/
# microclimlearn: this is microclimf's own sanctioned mechanism, distinct
# from the separate `mesoclim` package's more elaborate downscaling, which
# MicroclimateMapper.jl has no equivalent counterpart for and would make
# this an unfair comparison) plus its own slope/aspect/horizon-angle solar
# and wind-shelter correction (`runmicro`'s grid machinery, already
# confirmed real from `.modelina`'s source).
#
# Unlike grid_scaling.jl (flat terrain, CRUCL2 for both DEM and weather --
# a pure pixel-count timing test with no real spatial structure),
# this uses SRTM (fine, real elevation) as the terrain template with CRUCL2
# (coarse, ~18km) as the weather source and `compute_terrain = true` --
# matching demos/monthly_climate.jl's own SRTM-DEM/CRUCL2-weather example.
# The spatial detail visible in the output comes from MicroclimateMapper.jl's
# own lapse-rate elevation correction and slope/aspect/horizon-angle solar
# correction downscaling CRUCL2's coarse climate onto SRTM's fine terrain --
# not from CRUCL2 itself, which is flat at this resolution.
#
# Timed as the whole pipeline (DEM fetch + weather fetch + terrain build +
# lapse-rate/terrain-corrected solve), not just the solve step, since the
# terrain-correction machinery is exactly what's being exercised here.

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Rasters: X, Y, Ti
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW

gr()

const SNOTEL_DIR = joinpath(@__DIR__, "..", "comparisons", "scan_snotel")
include(joinpath(SNOTEL_DIR, "utils.jl"))
include(joinpath(SNOTEL_DIR, "config.jl"))
include(joinpath(SNOTEL_DIR, "pipeline.jl"))
include(joinpath(@__DIR__, "config.jl"))

outputs_dir = joinpath(@__DIR__, "outputs")
mkpath(outputs_dir)

meta_all = DataFrame(CSV.File(joinpath(SNOTEL_DIR, "Map metadata export.csv")))
meta_all[!, :ID] = strip.(string.(meta_all[!, :ID]))
site_row = filter(r -> strip(string(r.ID)) == string(site_num), meta_all)[1, :]
lat_dd = parse(Float64, string(site_row.Latitude))
lon_dd = parse(Float64, string(site_row.Longitude))

# SRTM is 3 arcsecond = 1/1200 degree/pixel (uniform in degrees, same as
# CRUCL2's own angular sampling) -- a bare degree value here badly
# undersells how many pixels that becomes (0.15 degrees half-width = ~130000
# pixels: (2*0.15)/(1/1200) squared). Size from a target pixel count instead,
# same approach grid_scaling.jl uses for CRUCL2. compute_terrain=true is real
# per-pixel slope/aspect/horizon-angle work (pricier than the flat CRUCL2
# case that already needed care about grid size) -- start small.
target_n_topo = 30
half_deg_area = target_n_topo * (1 / 1200) / 2
area = Extent(X = (lon_dd - half_deg_area, lon_dd + half_deg_area),
              Y = (lat_dd - half_deg_area, lat_dd + half_deg_area))

# CRUCL2's Monthly calendar always solves all 12 representative days
# regardless of the requested range (see grid_scaling.jl's header for why),
# and a narrower range breaks its internal solar-geometry array sizing --
# match demos/monthly_climate.jl's own full-year CRUCL2 pattern.
dates_full_year = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)

# Representative day/hour to snapshot -- July (day 7 of 12), local noon.
# Ti index follows DEFAULT_DAYS' month order: (month-1)*24 + hour + 1.
snap_month, snap_hour = 7, 12
snap_ti = (snap_month - 1) * 24 + snap_hour + 1

micro_model, soil_profile, _ = build_micro_model()

# Diagnostic: swap dynamic soil moisture for prescribed (CPCSoil, real
# observed data, time-varying -- matches demos/monthly_climate.jl's own
# `soil_moisture_source = CPCSoil` usage), to check whether the scattered
# "hot pixel" outliers seen in the terrain-corrected output are an artifact
# of the dynamic hydraulic solver interacting with terrain-corrected forcing
# at specific pixels, rather than the radiation/lapse-rate correction itself.
# Everything else (soil hydraulic calibration, depths, snow model) stays
# exactly as build_micro_model() set it up -- a controlled, single-variable
# change, not a switch to a different soil setup.
micro_model = MicroModel(;
    hours                    = micro_model.hours,
    depths                   = micro_model.depths,
    heights                  = micro_model.heights,
    soil_properties_model    = micro_model.soil_properties_model,
    soil_hydraulic_model     = micro_model.soil_hydraulic_model,
    radiation                = micro_model.radiation,
    snow_model               = micro_model.snow_model,
    vapour_pressure_equation = micro_model.vapour_pressure_equation,
    boundary_layer_model     = micro_model.boundary_layer_model,
    evaporation_model        = micro_model.evaporation_model,
    soil_energy_model        = micro_model.soil_energy_model,
    config = MicroConfig(;
        convergence            = micro_model.config.convergence,
        rainfall_schedule      = micro_model.config.rainfall_schedule,
        soil_moisture_strategy = PrescribedSoilMoisture(),
        max_surface_pool       = micro_model.config.max_surface_pool,
    ),
)

mapper_model_topo = MicroMapModel(;
    micro_model,
    dem_source              = SRTM,
    weather_source          = CRUCL2,
    # No soil_moisture_source: CPCSoil is a ~0.5deg global product, far
    # coarser than this ~0.025deg test area -- cropping it down that far hits
    # a 0-length range in Rasters' resample/warp path. PrescribedSoilMoisture
    # (set on micro_model above) falls back to the constant init.soil_moisture
    # given to MicroRasterProblem below, which is enough for the actual
    # diagnostic question (does removing the *dynamic* hydraulic solver change
    # the hot-pixel pattern) without needing CPCSoil's real fetched values.
    surface_albedo_source   = albedo,
    roughness_height_source = roughness_height,
    compute_terrain          = true,
)

println("=== SRTM (fine terrain) + CRUCL2 (coarse climate) around ($lat_dd, $lon_dd) ===")
jl_time = @elapsed begin
    raster_problem = MicroRasterProblem(;
        model = mapper_model_topo, area, dates = dates_full_year, template = SRTM,
        soil_profile, init = (; soil_moisture = fill(0.2, length(depths))))
    jl_result = MicroclimateMapper.solve(raster_problem)
end
nx, ny = length(dims(jl_result, X)), length(dims(jl_result, Y))
@printf("  MicroclimateMapper.jl (SRTM+CRUCL2, terrain-corrected): %d x %d pixels, %.2f s\n", nx, ny, jl_time)

# ── Snapshot the 4 requested variables at the chosen hour ────────────────────
jl_soilT = ustrip.(u"°C", jl_result.soil_temperature[depth = 1, Ti(snap_ti)])
jl_airT  = ustrip.(u"°C", jl_result.air_temperature[height = 1, Ti(snap_ti)])
jl_soilm = collect(jl_result.soil_moisture[depth = 1, Ti(snap_ti)])
jl_snow  = ustrip.(u"cm", jl_result.snow_depth[Ti(snap_ti)])

# Elevation, purely as context for why the above four vary spatially --
# loaded separately since MicroclimateMapper.solve()'s return doesn't carry
# the DEM itself, only derived output layers (matches demos/monthly_climate.jl's
# own separate `load_template` call for plotting the DEM).
srtm_template = load_template(SRTM, area)
jl_elev = collect(srtm_template)

timing_row = DataFrame(
    site = [site_num], nx = [nx], ny = [ny], n_points = [nx * ny],
    microclimatemapper_s = [jl_time], run_at = [string(now())],
)
CSV.write(joinpath(outputs_dir, "topo_timing_julia.csv"), timing_row; append = isfile(joinpath(outputs_dir, "topo_timing_julia.csv")))

# ── Diagnostics: is the elevation-driven gradient actually there, and are
# the standout hot pixels physically sane (steep/low-elevation) or arbitrary
# (would point at a bug rather than a real microsite)? A handful of extreme
# pixels alone will stretch a plain min/max color scale enough to make a
# real but subtler lapse-rate gradient look flat everywhere else.
println("\n── Diagnostics ──")
@printf("  Elevation range:  %.0f to %.0f m\n", minimum(jl_elev), maximum(jl_elev))
@printf("  Soil T range:     %.1f to %.1f °C (median %.1f)\n", minimum(jl_soilT), maximum(jl_soilT), median(jl_soilT))
@printf("  Air T range:      %.1f to %.1f °C (median %.1f)\n", minimum(jl_airT), maximum(jl_airT), median(jl_airT))
@printf("  corr(elevation, soil T) = %.2f   corr(elevation, air T) = %.2f   (negative = colder at higher elevation, as lapse rate would predict)\n",
        cor(vec(jl_elev), vec(jl_soilT)), cor(vec(jl_elev), vec(jl_airT)))

# Simple local-slope proxy (gradient magnitude of elevation) to check whether
# the hottest pixels sit on steep ground -- a plausible, real reason for
# extreme local heating (strong direct-beam exposure, little wind cooling),
# vs. scattered/edge pixels, which would look more like an artifact.
slope_proxy = zeros(size(jl_elev))
for i in 2:(size(jl_elev, 1) - 1), j in 2:(size(jl_elev, 2) - 1)
    slope_proxy[i, j] = sqrt((jl_elev[i+1, j] - jl_elev[i-1, j])^2 + (jl_elev[i, j+1] - jl_elev[i, j-1])^2)
end
n_top = min(5, length(jl_soilT))
hot_idx = partialsortperm(vec(jl_soilT), 1:n_top; rev = true)
println("  Hottest $n_top soil-T pixels (row, col, T °C, elevation m, slope-proxy):")
for k in hot_idx
    i, j = Tuple(CartesianIndices(jl_soilT)[k])
    @printf("    (%3d,%3d)  T=%.1f°C  elev=%.0fm  slope~%.0f\n", i, j, jl_soilT[i, j], jl_elev[i, j], slope_proxy[i, j])
end

# Percentile-clipped color limits so a few extreme pixels don't wash out the
# gradient across the rest of the grid.
clip(v, lo = 0.02, hi = 0.98) = (quantile(vec(v), lo), quantile(vec(v), hi))

fig = plot(layout = (2, 3), size = (1500, 900), dpi = 120)
heatmap!(fig[1], jl_elev;  title = "Elevation (m)", color = :terrain, yflip = true)
heatmap!(fig[2], jl_soilT; title = "Soil surface T (°C), clipped 2-98%ile", color = :thermal, yflip = true, clims = clip(jl_soilT))
heatmap!(fig[3], jl_airT;  title = "Air T, near-surface (°C), clipped 2-98%ile", color = :thermal, yflip = true, clims = clip(jl_airT))
heatmap!(fig[4], jl_soilm; title = "Soil moisture (m³/m³)", color = :viridis, yflip = true)
heatmap!(fig[5], jl_snow;  title = "Snow depth (cm)", color = :ice, yflip = true)
plot!(fig, plot_title = "MicroclimateMapper.jl — SRTM+CRUCL2, terrain-corrected " *
    "($(nx)x$(ny), $(round(jl_time, digits=1))s) — site $site_num, day $snap_month hour $snap_hour")

fig_path = joinpath(outputs_dir, "topo_julia_$(site_num).png")
savefig(fig, fig_path)
display(fig)
println("\nPlot saved to $fig_path")
println("If this doesn't show real spatial variation (e.g. the site's relief is " *
        "too subtle at $(half_deg_area)° half-width), widen half_deg_area or pick a " *
        "more mountainous site_num before moving on to the microclimf side.")

# ══════════════════════════════════════════════════════════════════════════
# STAGE 2 — Export a coarse, pre-correction CRUCL2 forcing grid + elevation
# reference for microclimf to independently downscale.
#
# Our SRTM test area (~0.025 deg span) is smaller than one CRUCL2 pixel
# (1/6 deg), so the "coarse grid" fetched here is deliberately buffered a
# bit wider (by one CRUCL2 pixel each side) to get a handful of coarse
# points rather than a degenerate single pixel -- dtmc's job is just to give
# each fine SRTM pixel a real reference elevation to lapse-correct against
# (elevd = dtmc_resampled - dtm), so a few coarse points spanning a slightly
# wider area than the fine target is exactly what's wanted, not a mismatch.
#
# Point-solved with compute_terrain=false (flat) and CRUCL2 as both
# dem_source and weather_source, so what's exported is the coarse climate
# itself, not something Julia has already terrain-corrected for us.
# ══════════════════════════════════════════════════════════════════════════

coarse_buffer = 1 / 6   # one CRUCL2 pixel
coarse_area = Extent(X = (lon_dd - half_deg_area - coarse_buffer, lon_dd + half_deg_area + coarse_buffer),
                      Y = (lat_dd - half_deg_area - coarse_buffer, lat_dd + half_deg_area + coarse_buffer))
crucl2_template = load_template(CRUCL2, coarse_area)
coarse_xs = collect(lookup(crucl2_template, X))
coarse_ys = collect(lookup(crucl2_template, Y))
n_cx, n_cy = length(coarse_xs), length(coarse_ys)
println("\nCoarse CRUCL2 grid for export: $n_cx x $n_cy pixels")

# Extra output layers beyond the defaults -- pressure and diffuse_fraction
# are real MicroResult fields (confirmed from the original Microclimate.jl
# API exploration) but aren't in MicroclimateMapper's default output set.
coarse_output_layers = (
    LayerSpec(:air_temperature, :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:pressure, :scalar),
    LayerSpec(:diffuse_fraction, :scalar),
)

mapper_model_coarse = MicroMapModel(;
    micro_model,
    dem_source               = CRUCL2,
    weather_source           = CRUCL2,
    surface_albedo_source    = albedo,
    roughness_height_source  = roughness_height,
    compute_terrain           = false,
    output_layers             = coarse_output_layers,
)

coarse_points = [GIW.Point((x, y)) for y in coarse_ys for x in coarse_xs]
coarse_vec_problem = MicroVectorProblem(;
    model = mapper_model_coarse, points = coarse_points, dates = dates_full_year,
    soil_profile, init = (; soil_moisture = fill(0.2, length(depths))))
coarse_out = MicroclimateMapper.solve(coarse_vec_problem)

sigma_sb = 5.670374419e-8
# `end` doesn't reliably resolve to the *named* dimension's own size when
# mixed with a positional Ti(...) selector in the same indexing expression
# (confirmed: it fell back to the array's first-dimension size instead) --
# use an explicit index into `profile_heights` instead. Index `ref_height_idx`
# = the reference (tallest) height, matching what "end" was meant to select.
ref_height_idx = length(profile_heights)
temp_c   = ustrip.(u"°C", coarse_out.air_temperature[height = ref_height_idx, Ti(snap_ti)])
relhum_c = ustrip.(coarse_out.relative_humidity[height = ref_height_idx, Ti(snap_ti)]) .* 100.0
wind_c   = ustrip.(u"m/s", coarse_out.wind_speed[height = ref_height_idx, Ti(snap_ti)])
pres_c   = ustrip.(u"kPa", uconvert.(u"kPa", coarse_out.pressure[Ti(snap_ti)]))
sw_c     = ustrip.(u"W/m^2", coarse_out.global_radiation[Ti(snap_ti)])
difrad_c = sw_c .* ustrip.(coarse_out.diffuse_fraction[Ti(snap_ti)])
lw_c     = sigma_sb .* ustrip.(u"K", coarse_out.sky_temperature[Ti(snap_ti)]).^4

# Precip: not part of MicroResult's output (it's forcing, not a model
# result), and this snapshot's core purpose (temperature/radiation response
# to terrain correction) isn't sensitive to it -- left at zero rather than
# chasing a separate CRUCL2 :pre band read for a variable that won't
# materially change this comparison. If a rain-sensitive snapshot is wanted
# later, this is the place to add a real value.
precip_c  = zeros(n_cx * n_cy)
winddir_c = zeros(n_cx * n_cy)   # flat terrain -- direction has no local effect at this stage

# ── Write GeoTIFFs for R to read (terra::rast reads these directly) ─────────
topo_dir = joinpath(outputs_dir, "topo_$(site_num)")
mkpath(topo_dir)

function write_coarse(name, vec_vals)
    # coarse_points iterates x fastest (inner loop), y slowest -- matches
    # Julia's column-major reshape (first dim varies fastest), and
    # crucl2_template's own dims are (X, Y) in that same fast/slow order, so
    # no transpose should be needed. NOT independently verified against a
    # live run -- if the written GeoTIFFs come out transposed/flipped when
    # read back in R, this reshape (not the point-list order) is the first
    # place to check.
    r = Raster(reshape(vec_vals, n_cx, n_cy), dims(crucl2_template); crs = crs(crucl2_template))
    Rasters.write(joinpath(topo_dir, "coarse_$name.tif"), r; force = true)
end
write_coarse("temp", temp_c);   write_coarse("relhum", relhum_c); write_coarse("pres", pres_c)
write_coarse("swdown", sw_c);   write_coarse("difrad", difrad_c); write_coarse("lwdown", lw_c)
write_coarse("windspeed", wind_c); write_coarse("winddir", winddir_c); write_coarse("precip", precip_c)
Rasters.write(joinpath(topo_dir, "dtmc.tif"), crucl2_template; force = true)
Rasters.write(joinpath(topo_dir, "dtm.tif"), srtm_template; force = true)

# Same soil/veg translation already established in config.jl (MCF_* globals)
# -- reused here rather than re-specified in R, so the topographic-correction
# comparison uses the identical calibration as the point/grid-timing scripts.
topo_params = DataFrame(
    height1_m = [ustrip(u"m", profile_heights[1])],
    smax = [MCF_SMAX], b = [MCF_B], psi_e = [MCF_PSI_E], ksat = [MCF_KSAT], groundr = [MCF_GROUNDR],
    rho_Mgm3 = [ustrip(u"Mg/m^3", uconvert(u"Mg/m^3", bulk_density))],
    pai = [MCF_PAI], hgt = [MCF_HGT], x = [MCF_X], gsmax = [MCF_GSMAX],
    clump = [MCF_CLUMP], leafr = [MCF_LEAFR], leafd = [MCF_LEAFD], leaft = [MCF_LEAFT],
)
CSV.write(joinpath(topo_dir, "topo_params.csv"), topo_params)

println("Coarse forcing grid + dtmc + fine dtm + params written to $topo_dir")

# ══════════════════════════════════════════════════════════════════════════
# STAGE 3 — Run microclimf's own topographic correction on that coarse grid,
# time it, and compare against Julia's terrain-corrected snapshot.
# ══════════════════════════════════════════════════════════════════════════

rscript = joinpath(@__DIR__, "run_microclimf_topo.R")
println("\nRunning microclimf's own topographic correction (Rscript)...")
r_time = @elapsed run(`Rscript $rscript $topo_dir`)
r_timing = DataFrame(CSV.File(joinpath(topo_dir, "r_timing.csv")))
r_solve_s = r_timing.elapsed_s[1]
@printf("  microclimf (altcorrect + terrain grid): %.2f s (headline solve: %.2f s)\n", r_time, r_solve_s)

speedup = r_solve_s / jl_time
faster  = speedup >= 1 ? "MicroclimateMapper.jl" : "microclimf"
@printf("  %s is %.1fx faster on this comparison\n", faster, max(speedup, 1 / speedup))

timing_row2 = DataFrame(
    site = [site_num], microclimatemapper_s = [jl_time], microclimf_s = [r_solve_s],
    microclimf_wall_s = [r_time], run_at = [string(now())],
)
timing_file2 = joinpath(outputs_dir, "topo_timing_comparison.csv")
CSV.write(timing_file2, timing_row2; append = isfile(timing_file2))

# ── Read microclimf's snapshots back and plot side by side with Julia's ────
# Plain numeric GeoTIFFs from R -- already °C / m³/m³, no Unitful units
# attached, so just collect into plain arrays (no uconvert/ustrip needed).
r_airT  = collect(Raster(joinpath(topo_dir, "r_airT.tif")))
r_soilT = collect(Raster(joinpath(topo_dir, "r_soilT.tif")))
r_soilm = collect(Raster(joinpath(topo_dir, "r_soilm.tif")))

fig2 = plot(layout = (2, 3), size = (1500, 900), dpi = 120)
heatmap!(fig2[1], jl_soilT; title = "Julia soil surface T (°C)", color = :thermal, yflip = true, clims = clip(jl_soilT))
heatmap!(fig2[2], jl_airT;  title = "Julia air T, near-surface (°C)", color = :thermal, yflip = true, clims = clip(jl_airT))
heatmap!(fig2[3], jl_soilm; title = "Julia soil moisture (m³/m³)", color = :viridis, yflip = true)
heatmap!(fig2[4], r_soilT;  title = "microclimf soil surface T (°C)", color = :thermal, yflip = true, clims = clip(r_soilT))
heatmap!(fig2[5], r_airT;   title = "microclimf air T, near-surface (°C)", color = :thermal, yflip = true, clims = clip(r_airT))
heatmap!(fig2[6], r_soilm;  title = "microclimf soil moisture (m³/m³)", color = :viridis, yflip = true)
plot!(fig2, plot_title = "Topographic correction — MicroclimateMapper.jl ($(round(jl_time,digits=1))s) vs " *
    "microclimf ($(round(r_solve_s,digits=1))s) — site $site_num, day $snap_month hour $snap_hour")

fig2_path = joinpath(outputs_dir, "topo_comparison_$(site_num).png")
savefig(fig2, fig2_path)
display(fig2)
println("\nComparison plot saved to $fig2_path")
println("Timing saved to $timing_file2")
