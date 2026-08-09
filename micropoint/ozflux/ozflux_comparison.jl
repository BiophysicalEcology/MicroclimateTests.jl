# ozflux_comparison.jl — entry point: runs Microclimate.jl (comparisons/
# ozflux's own run_site_gapfilled) and micropoint (R, RunModelFull) on the
# identical real forcing, then reports 3-way stats/plots (obs, Microclimate.jl,
# micropoint) via report_site_results3 (report.jl), plus a solve-time
# comparison. Self-contained (not added to comparisons/ozflux/report.jl) so
# the main ozflux comparison pipeline can't regress.
#
# Usage: julia ozflux_comparison.jl (or include() from the REPL)

using Microclimate, MicroclimateMapper, Unitful, FluidProperties
using CSV, DataFrames, Dates, Statistics, Printf, Plots, Measures
using Rasters, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

const OZFLUX_DIR = joinpath(@__DIR__, "..", "..", "comparisons", "ozflux")
include(joinpath(OZFLUX_DIR, "config.jl"))
include(joinpath(OZFLUX_DIR, "utils.jl"))
include(joinpath(OZFLUX_DIR, "pipeline.jl"))
include(joinpath(OZFLUX_DIR, "report.jl"))
include(joinpath(@__DIR__, "write_ozflux_micropoint_inputs.jl"))
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "report.jl"))

# site_name: any comparisons/ozflux/data site with a SITE_UTC_OFFSET_MINUTES
# entry (config.jl) -- "CapeTribulation", "Calperum", "Whroo", "Wallaby",
# "GWW", "Longreach", "TiTreeEast", "AliceSpringsMulga".
site_name = "Calperum"  # "Calperum"  # "Whroo"  # "Wallaby"  # "GWW"  # "Longreach"  # "TiTreeEast"  # "AliceSpringsMulgaa"
years = [2015]
outdir = joinpath(@__DIR__, "outputs", "$(site_name)_$(join(years, '-'))")
save_outputs_3way = true
display_plots_3way = false
plot_start = Date(years[1], 1, 1)
plot_end   = Date(years[1], 12, 30)
# Arbitrary snapshot window for the vertical-profile plots -- needs to land
# in a cleanly-solved chunk (see README); empty panels just mean pick a
# different date, not an error.
profile_times = DateTime(years[1], 11, 1, 4):Hour(1):DateTime(years[1], 11, 1, 15)

# ── Run both models on identical forcing ─────────────────────────────────
result = run_site_gapfilled(site_name, years; canopy_mode=:full)
write_ozflux_micropoint_inputs(result, outdir)

rscript = joinpath(@__DIR__, "run_micropoint_ozflux_vegetated.R")
println("\nRunning micropoint (vegetated) on the same forcing...")
run(`Rscript $rscript $outdir`)

timing = CSV.read(joinpath(outdir, "veg_timing.csv"), DataFrame)
println("\n== Timing (model solve only -- excludes SILO gap-fill donor solve, ")
println("   SLGA/data fetch, and Rscript subprocess startup, for both sides) ==")
@printf("  Julia (Microclimate.jl, Microclimate.solve only): %.2f s\n", result.solve_time)
@printf("  R (micropoint RunModelFull, internal system.time only): %.2f s\n", timing.elapsed_s[1])
CSV.write(joinpath(outdir, "julia_timing.csv"),
    DataFrame(; solve_time_s=result.solve_time, convergence_model=string(typeof(canopy_convergence_model_choice)),
        air_profile_model=string(typeof(canopy_air_profile_model_choice))))

# rows == climdata.csv's own order == result.t_model, 1:1 -- no spin-up
# period here, so no realignment needed.
mp = read_micropoint_output(outdir)

stats_df = report_site_results3(result, mp; outdir, plot_start, plot_end, profile_times, save_outputs_3way, display_plots_3way)
println("\nDone. Outputs in $outdir")
