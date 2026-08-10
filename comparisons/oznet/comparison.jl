# comparison.jl
#
# Julia-first microclimate comparison against OzNet soil moisture/temperature
# field observations (Smith et al. 2012, http://www.oznet.org.au/), the same
# dataset used in the original NicheMapR-based study (Kearney & Maino 2018)
# whose R prototype lives alongside this folder (Australia_Soil_Test.Rmd,
# compare.R, plot_moist.R, plot_temp.R, soilprops.txt — kept for reference/
# provenance; this pipeline doesn't read any of them). Forcing is fetched via
# MicroclimateMapper's weather-source support (see `weather_source_choice` in
# config.jl — SILO/BARRA/NCEP{SurfaceFlux}); soil properties come from a live
# SLGA fetch + Cosby pedotransfer per site. NicheMapR hourly outputs are
# loaded as a secondary comparison if present in nmr_out_dir.
#
# All model setup and per-site logic lives in pipeline.jl (shared with
# single_site.jl) and config.jl (shared model/soil-source/weather-source
# parameters — change `weather_source_choice`, `pedotransfer_model_choice`,
# etc. there to compare variants).
#
# Requires the OzNet dataset (oznetsiteinfo2.csv + per-site *_sm.txt files)
# at OZNET_DATA_DIR (config.jl) — see README.md; not committed to the repo.
#
# Usage:
#   julia comparison.jl              ← runs sites list defined below
#   include("comparison.jl")         ← from a Julia REPL

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization
using Random

include("config.jl")
include("utils.jl")
include("pipeline.jl")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  RUN CONFIGURATION — what varies per run, as opposed to the shared model ║
# ║  parameters in config.jl                                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── Sites ─────────────────────────────────────────────────────────────────────
#   :all           — every site in oznetsiteinfo2.csv
#   ["a1","k6"]     — only these site names
#   5               — randomly sample this many sites (reproducible via
#                      site_subset_seed) — handy for a quick multi-site test
site_subset      = :all #["a1"]
site_subset_seed = 1234

# ── Simulation period ─────────────────────────────────────────────────────────
# Used when auto_date_range = false.
sim_start = Date(2007, 1, 1)
sim_end   = Date(2010, 12, 31)

# auto_date_range = true: detect sim_start/sim_end from each site's own obs
# file date span (intersected with the chosen weather source's availability
# — WEATHER_SOURCE_START in config.jl — through yesterday), capped to
# max_sim_years (config.jl).
auto_date_range = true

# Plot window — set either to nothing to show the full simulation period.
plot_start = Date(2007, 1, 1)
plot_end   = Date(2007, 7, 31)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN LOOP                                                               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

all_stats = DataFrame()
siteinfo = DataFrame(CSV.File(OZNET_SITEINFO))

rows = collect(eachrow(siteinfo))
if site_subset isa Integer
    Random.seed!(site_subset_seed)
    rows = rows[sort(randperm(length(rows))[1:min(site_subset, length(rows))])]
    println("site_subset: randomly sampled $(length(rows)) of $(nrow(siteinfo)) site(s) (seed=$site_subset_seed).")
elseif site_subset != :all
    rows = filter(r -> r.name in site_subset, rows)
    println("site_subset: restricted to $(length(rows)) of $(nrow(siteinfo)) requested site(s).")
end
println("Running $(length(rows)) OzNet site(s).")

@isdefined(weather_cache) || (global weather_cache = Dict{Tuple{Symbol,String,Date,Date}, Any}())

make_plots_flag = save_outputs || length(rows) == 1
display_plots_flag = length(rows) == 1

# ── Phase 1: prepare every site (weather fetch, SLGA soil profile, obs,
#    initial conditions, build the MicroProblem). Not batched across sites
#    like scan_snotel — OzNet is 38 sites, not 600+, and each site's own
#    depth-node profile differs, so there's no shared soil_profile to
#    amortise a batched weather fetch against. ───────────────────────────────
preps = []
for row in rows
    try
        prep = prepare_site(row, weather_cache; sim_start, sim_end, auto_date_range, max_sim_years)
        isnothing(prep) || push!(preps, prep)
    catch _err
        @warn "Site $(row.name) failed to prepare: $(_err)"
    end
end
println("\n$(length(preps)) of $(length(rows)) site(s) prepared and ready to run.")

# ── Phase 2: run all Julia models (parallel across Threads.nthreads()) and all
#    NicheMapR runs (concurrent Rscript processes, capped) ────────────────────
julia_results = run_julia_batch!(preps)
run_nmr_batch!(preps)

# ── Phase 3: stats + plots per site ───────────────────────────────────────────
for prep in preps
    try
        micro_out, julia_solve_time = julia_results[prep.site_name]
        site_rows = report_site_results(prep, micro_out, julia_solve_time;
            plot_start, plot_end, make_plots = make_plots_flag, display_plots = display_plots_flag)
        isnothing(site_rows) || append!(all_stats, DataFrame(site_rows))
    catch _err
        @warn "Site $(prep.site_name) failed to report: $(_err)"
    end
end

# ── Cross-site summary (multi-site runs) ──────────────────────────────────────
if length(rows) > 1 && nrow(all_stats) > 0
    println("\n" * "="^72)
    println("SUMMARY — $(length(rows)) sites — weather source: $(nameof(weather_source_choice))")
    println("="^72)
    for kind in ("temperature", "moisture")
        sub_kind = filter(r -> r.kind == kind, all_stats)
        nrow(sub_kind) == 0 && continue
        println("\n$(uppercasefirst(kind)):")
        println(rpad("Depth (cm)", 12) * "  " *
                "Julia: mean r  mean RMSE  mean bias     " *
                "NMR:  mean r  mean RMSE  mean bias")
        println("-"^90)
        for d in sort(unique(sub_kind.depth_cm))
            sub = filter(r -> r.depth_cm == d, sub_kind)
            fmtv(v) = all(isnan, v) ? "  —  " : @sprintf("%.4g", mean(filter(isfinite, v)))
            jr  = all(isnan, sub.j_r)     ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.j_r)))
            jb  = all(isnan, sub.j_bias)  ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.j_bias)))
            nr_ = all(isnan, sub.nmr_r)   ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.nmr_r)))
            nb  = all(isnan, sub.nmr_bias) ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.nmr_bias)))
            println("  $(rpad(d,10))  $jr    $(fmtv(sub.j_rmse))    $jb      $nr_    $(fmtv(sub.nmr_rmse))    $nb")
        end
    end
end
