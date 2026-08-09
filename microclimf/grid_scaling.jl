# grid_scaling.jl — times MicroclimateMapper.jl's raster mode against
# microclimf's native grid mode (runmicro) on one reasonably large grid.
# Solve-time comparison only — no accuracy comparison, no observations
# (none exist at grid scale). Secondary/exploratory deliverable relative to
# comparison.jl — see README.md's caveats.
#
# Requires comparison.jl to have been run at least once first (reuses its
# climdata.csv, sliced to the first 24 hours, as microclimf's forcing —
# timing is the point here, not matching the two models' data sources or
# spatial realism).
#
# Uses CRUCL2 (CRU CL 2.0 monthly climatology, 10-arcminute resolution,
# ~18 km) rather than GRIDMET/SRTM for the Julia raster: no per-day
# historical fetch (one pre-downloaded global NetCDF, cropped per area),
# matching demos/monthly_climate.jl's own CRUCL2 raster pattern
# (dem_source = weather_source = CRUCL2).
#
# CRUCL2's `Monthly` calendar (MicroclimateMapper.jl/src/climate/weather.jl)
# always solves all `days_per_year(Monthly()) = 12` representative days
# (one per month) for whatever calendar year(s) the requested `dates` span —
# `_years_from_dates` collapses any date range down to its touched years,
# and `_days_of_year(Monthly(), years) = repeat(DEFAULT_DAYS, length(years))`
# doesn't look inside a year at all. So narrowing `dates` to a single day
# does NOT reduce the per-pixel solve to 1 day — it's always 12 days
# (288 hours) per pixel regardless. Grid size (`target_n`) is therefore the
# only real lever for keeping total memory/solve cost down.
#
# MicroclimateMapper.jl's `MicroRasterProblem` sizes its raster off the DEM
# source's own resolution cropped to `area` — no direct "give me exactly
# NxN pixels" knob. CRUCL2's resolution is 10 arcminutes = 1/6 degree per
# pixel (uniform in degrees, like SRTM's angular sampling), so for a target
# of `n` pixels per side: `half_deg = n * (1/6) / 2`. The actual resulting
# size is still read back from the solved `RasterStack`'s own (X, Y) dims
# rather than assumed, and microclimf's synthetic grid is matched to that
# real count.

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Rasters: X, Y
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW

gr()

const SNOTEL_DIR = joinpath(@__DIR__, "..", "comparisons", "scan_snotel")
include(joinpath(SNOTEL_DIR, "utils.jl"))
include(joinpath(SNOTEL_DIR, "config.jl"))
include(joinpath(SNOTEL_DIR, "pipeline.jl"))
include(joinpath(@__DIR__, "config.jl"))

isfile(joinpath(MCF_OUTDIR, "climdata.csv")) ||
    error("Run comparison.jl first — grid_scaling.jl reuses its climdata.csv as forcing.")

outputs_dir = joinpath(@__DIR__, "outputs")
mkpath(outputs_dir)

# Reuse the same site-calibrated micro_model/soil_profile as the point
# comparison; only dem_source/weather_source change, to CRUCL2.
micro_model, soil_profile, _ = build_micro_model()
mapper_model_cru = MicroMapModel(;
    micro_model,
    dem_source              = CRUCL2,
    weather_source          = CRUCL2,
    surface_albedo_source   = albedo,
    roughness_height_source = roughness_height,
    compute_terrain          = false,
)

meta_all = DataFrame(CSV.File(joinpath(SNOTEL_DIR, "Map metadata export.csv")))
meta_all[!, :ID] = strip.(string.(meta_all[!, :ID]))
site_row = filter(r -> strip(string(r.ID)) == string(site_num), meta_all)[1, :]
lat_dd = parse(Float64, string(site_row.Latitude))
lon_dd = parse(Float64, string(site_row.Longitude))

# Every pixel solves 12 representative days (288 hours) regardless of grid
# size (see header) -- 100x100 crashed (likely OOM: 10000 pixels x 288
# hours x several depth-resolved variables). Start much smaller and scale
# up only after confirming a size actually completes.
target_n = 20
half_deg = target_n * (1 / 6) / 2

area = Extent(X = (lon_dd - half_deg, lon_dd + half_deg), Y = (lat_dd - half_deg, lat_dd + half_deg))
# CRUCL2 is a 1961-1990 climatology -- the year is ignored, and narrowing
# this below a full year doesn't reduce the actual solve (see header: always
# 12 representative days regardless) -- it just breaks solar-geometry array
# sizing internally (confirmed: a single-day range throws a BoundsError
# indexing day 2 of a 1-element solar array). Match demos/monthly_climate.jl,
# which never uses anything narrower than a full year for CRUCL2.
dates_range = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)

println("=== Grid comparison: targeting ~$(target_n)x$(target_n), ±$(round(half_deg, digits=3))° box around ($lat_dd, $lon_dd), CRUCL2 ===")
# CRUCL2 needs a pre-loaded template (unlike SRTM, which accepts the bare
# type and auto-loads/crops) -- its own load_template method reads the :elv
# band directly since the generic RasterDataSources extent-keyword path
# doesn't apply to its single-file format (see crucl2.jl's own comment).
crucl2_template = load_template(CRUCL2, area)
# CRUCL2 doesn't provide soil_moisture as a canonical variable (unlike e.g.
# TerraClimate), so it must be supplied explicitly -- matches
# demos/monthly_climate.jl's own CRUCL2 raster example.
raster_problem = MicroRasterProblem(;
    model = mapper_model_cru, area, dates = dates_range, template = crucl2_template, soil_profile,
    init = (; soil_moisture = fill(0.2, length(depths))))

jl_time = @elapsed jl_result = MicroclimateMapper.solve(raster_problem)
nx, ny = length(dims(jl_result, X)), length(dims(jl_result, Y))
n_points = nx * ny
@printf("  MicroclimateMapper.jl: %d x %d = %d points, %.2f s\n", nx, ny, n_points, jl_time)

# Match microclimf's synthetic grid to the same actual pixel count. Floored
# at 10: `.windsheltera()` (R/internal.R), used for wind-shelter effects on
# every grid run (not just snow), aggregates a same-size copy of the DTM by
# a factor hardcoded to 10 whenever pixel resolution is <=100m (our case),
# so a smaller grid means aggregating by more than the raster's own size.
grid_n = max(round(Int, sqrt(n_points)), 10)
rscript = joinpath(@__DIR__, "run_microclimf_grid.R")
r_time = @elapsed run(`Rscript $rscript $MCF_OUTDIR $grid_n $lat_dd $lon_dd 24`)
r_timing = DataFrame(CSV.File(joinpath(MCF_OUTDIR, "microclimf_grid_timing.csv")))
r_solve_s = r_timing.elapsed_s[1]
@printf("  microclimf:             %d x %d = %d points, %.2f s (solve only: %.2f s)\n",
        grid_n, grid_n, grid_n^2, r_time, r_solve_s)

row = (n_points_julia = n_points, nx_julia = nx, ny_julia = ny,
       grid_n_microclimf = grid_n, n_points_microclimf = grid_n^2,
       microclimatemapper_s = jl_time,
       microclimf_s = r_solve_s, microclimf_wall_s = r_time)

df = DataFrame([row])
CSV.write(joinpath(outputs_dir, "grid_timing.csv"), df)
println("\nGrid timing saved to $(joinpath(outputs_dir, "grid_timing.csv"))")

speedup = r_solve_s / jl_time
faster  = speedup >= 1 ? "MicroclimateMapper.jl" : "microclimf"
@printf("  %s is %.1fx faster on this grid\n", faster, max(speedup, 1 / speedup))

fig = bar(["MicroclimateMapper.jl\n($(nx)x$(ny))", "microclimf\n($(grid_n)x$(grid_n))"],
    [jl_time, r_solve_s];
    color = [:tomato, :steelblue], legend = false, ylabel = "solve time (s)",
    title = "Grid solve time — site $site_num (models vs each other only)")
savefig(fig, joinpath(outputs_dir, "grid_timing.png"))
display(fig)
println("Plot saved to $(joinpath(outputs_dir, "grid_timing.png"))")
