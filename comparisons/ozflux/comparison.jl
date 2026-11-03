# comparison.jl
#
# Batch entry point: runs the canopy model on tower forcing at every site in
# data/, over each site's best-covered year or, where shorter, its whole record
# (sites.jl's run_years). Writes per-site stats (report.jl) to all_stats.csv and prints a cross-site summary. All model
# setup and per-site logic lives in pipeline.jl/sites.jl/report.jl (shared
# with single_site.jl) and config.jl. SILO forcing runs are single_site.jl.
#
# Usage:
#   julia comparison.jl              ← every site in data/
#   include("comparison.jl")         ← from a Julia REPL

using Microclimate, MicroclimateMapper, Unitful, FluidProperties, SolarRadiation, OrdinaryDiffEqSDIRK
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

include("ozflux.jl")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  RUN CONFIGURATION                                                       ║
# ╚══════════════════════════════════════════════════════════════════════════╝

sites = filter(!excluded, sort(collect(keys(discover_site_years(OZFLUX_DATA_DIR)))))
make_plots = false   # per-site plots for every run (slow over every site)
plot_forcing = false
plot_variable = nothing

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN LOOP                                                               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

all_stats = DataFrame()
failures = Tuple{String,String}[]
fail!(site_name, err) = (push!(failures, (site_name, sprint(showerror, err))); @warn "Site $site_name failed" exception=err)

# Prepare serially (NetCDF and the caches aren't thread-safe), solve on every thread, report serially
prepared = []
for site_name in sites
    try
        push!(prepared, announce_and_prepare_site(site_name, run_years(site_name)))
    catch err
        fail!(site_name, err)
    end
end
println("\nSolving $(length(prepared)) sites on $(Threads.nthreads()) threads...")
results = Vector{Any}(undef, length(prepared))
Threads.@threads :dynamic for i in eachindex(prepared)
    results[i] = try
        solve_site(prepared[i])
    catch err
        err
    end
end
for (prep, result) in zip(prepared, results)
    result isa Exception && (fail!(prep.site_name, result); continue)
    save_site_result(result, "tower")
    try
        stats = report_site_results(result, result.output; make_plots, display_plots=false, plot_forcing, plot_variable)
        stats.solve_time .= result.solve_time
        append!(all_stats, stats; cols=:union)
    catch err
        fail!(prep.site_name, err)
    end
end
CSV.write(joinpath(outputs_dir, "all_stats.csv"), all_stats)

# ── Cross-site summary: mean bias and RMSE of each variable over the sites measuring it ──
if nrow(all_stats) > 0
    println("\n" * "="^72)
    println("SUMMARY — $(length(unique(all_stats.site))) site(s) run, $(length(failures)) failed")
    println("="^72)
    sub = filter(r -> isfinite(r.rmse), all_stats)
    println(rpad("Variable", 34) * "  sites   bias / RMSE    units")
    for v in unique(sub.variable)
        rows = filter(r -> r.variable == v, sub)
        @printf("  %-32s  %5d   %+8.3f / %7.3f  %s\n", v, nrow(rows), mean(rows.bias), mean(rows.rmse), only(unique(rows.units)))
    end
    isempty(failures) || (println("\nFailed:"); foreach(((s, e),) -> println("  $s: ", first(split(e, '\n'))), failures))
end
