# comparison.jl
#
# Batch entry point: runs the canopy model on tower forcing at every site in
# data/, over each site's best-covered year (sites.jl's best_year), with and
# without litter where the site has it (woody canopies). Writes per-site stats
# (report.jl) to all_stats.csv and prints a cross-site summary. All model
# setup and per-site logic lives in pipeline.jl/sites.jl/report.jl (shared
# with single_site.jl) and config.jl. SILO forcing runs need mapper.jl, which
# needs MicroclimateMapper ported to the current Microclimate first.
#
# Usage:
#   julia comparison.jl              ← every site in data/
#   include("comparison.jl")         ← from a Julia REPL

using Microclimate, Unitful, FluidProperties, OrdinaryDiffEqSDIRK
using CSV, DataFrames, Dates, Statistics, Printf, Plots, NCDatasets

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

# Each site's litter, and none where that differs, so litter's effect shows at every site with it
litter_runs(site_name) = unique([site_litter(site_name), NoLitter()])
litter_label(litter_model) = litter_model isa NoLitter ? "no_litter" : "litter"

all_stats = DataFrame()
failures = Tuple{String,String}[]
for site_name in sites
    try
        years = [best_year(site_name)]
        for litter_model in litter_runs(site_name)
            result = run_site(site_name, years; litter_model)
            stats = report_site_results(result, result.output; make_plots, display_plots=false,
                source=litter_label(litter_model), plot_forcing, plot_variable)
            stats.solve_time .= result.solve_time
            append!(all_stats, stats; cols=:union)
        end
    catch err
        push!(failures, (site_name, sprint(showerror, err)))
        @warn "Site $site_name failed" exception=err
    end
end
CSV.write(joinpath(outputs_dir, "all_stats.csv"), all_stats)

# ── Cross-site summary: litter vs none at the sites run both ways ───────────
if nrow(all_stats) > 0
    both = intersect(unique(all_stats.site[all_stats.source .== "litter"]), unique(all_stats.site[all_stats.source .== "no_litter"]))
    println("\n" * "="^72)
    println("SUMMARY — $(length(unique(all_stats.site))) site(s) run, $(length(failures)) failed, $(length(both)) with and without litter")
    println("="^72)
    sub = filter(r -> r.site in both && isfinite(r.rmse), all_stats)
    println(rpad("Variable", 34) * "  sites   no-litter bias / RMSE    litter bias / RMSE")
    for v in unique(sub.variable)
        rows(source) = filter(r -> r.variable == v && r.source == source, sub)
        none, with = rows("no_litter"), rows("litter")
        common = intersect(none.site, with.site)
        isempty(common) && continue
        pick(df) = df[indexin(common, df.site), :]
        none, with = pick(none), pick(with)
        @printf("  %-32s  %5d   %+8.3f / %7.3f      %+8.3f / %7.3f\n", v, length(common),
            mean(none.bias), mean(none.rmse), mean(with.bias), mean(with.rmse))
    end
    isempty(failures) || (println("\nFailed:"); foreach(((s, e),) -> println("  $s: ", first(split(e, '\n'))), failures))
end
