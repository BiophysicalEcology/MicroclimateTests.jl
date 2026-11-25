# Entry point: runs every site on tower forcing over its run years, writes all_stats.csv and prints a cross-site summary

using Microclimate, MicroclimateMapper, Unitful, FluidProperties, SolarRadiation, OrdinaryDiffEqSDIRK
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Serialization

include("ozflux.jl")

sites = filter(!excluded, sort(collect(keys(discover_site_years(OZFLUX_DATA_DIR)))))
# Per-site plots, slow over every site
make_plots = false
plot_forcing = false
plot_variable = nothing
# Report on the saved results (save_site_result) instead of solving
rescore = get(ENV, "MICROCLIMATETESTS_RESCORE", "false") == "true"

all_stats = DataFrame()
failures = Tuple{String,String}[]
function fail!(site_name, err)
    push!(failures, (site_name, sprint(showerror, err)))
    @warn "Site $site_name failed" exception=err
end

# Prepare serially (NetCDF and the caches aren't thread-safe), solve on every thread, report serially
prepared = []
for site_name in sites
    try
        push!(prepared, rescore ? prepare_saved_site(site_name, run_years(site_name), "tower") :
            announce_and_prepare_site(site_name, run_years(site_name)))
    catch err
        fail!(site_name, err)
    end
end
results = Vector{Any}(undef, length(prepared))
if rescore
    results .= prepared
else
    println("\nSolving $(length(prepared)) sites on $(Threads.nthreads()) threads...")
    Threads.@threads :dynamic for i in eachindex(prepared)
        results[i] = try
            solve_site(prepared[i])
        catch err
            err
        end
    end
end
for (prep, result) in zip(prepared, results)
    if result isa Exception
        fail!(prep.site_name, result)
        continue
    end
    rescore || save_site_result(result, "tower")
    try
        stats = report_site_results(result, result.output; make_plots, display_plots=false, plot_forcing, plot_variable)
        stats.solve_time .= result.solve_time
        append!(all_stats, stats; cols=:union)
    catch err
        fail!(prep.site_name, err)
    end
end
mkpath(outputs_dir)
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
    isempty(failures) || println("\nFailed:")
    for (s, e) in failures
        println("  $s: ", first(split(e, '\n')))
    end
end
