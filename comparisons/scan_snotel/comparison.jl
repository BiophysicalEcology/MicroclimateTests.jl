# comparison.jl
#
# Julia-first microclimate comparison against SCAN/SNOTEL field observations.
# Forcing is fetched via MicroclimateMapper's weather-source support (see
# `weather_source_choice` in config.jl — GRIDMET/ERA5/NCEP{SurfaceFlux}).
# NicheMapR hourly outputs are loaded as a secondary comparison if present in
# nmr_out_dir.
#
# All model setup and per-site logic lives in pipeline.jl (shared with
# debug_site.jl) and config.jl (shared model/soil/snow parameters — change
# `infiltration_algorithm_choice`, `weather_source_choice`, etc. there to
# compare variants).
#
# Requires `Map metadata export.csv` (site lat/lon/elevation) and per-site
# observation CSVs under `observations/` — see README.md for how to obtain
# these; neither is committed to the repo.
#
# Usage:
#   julia comparison.jl              ← runs sites list defined below
#   include("comparison.jl")         ← from a Julia REPL

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
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
# auto_sites = true : scan observations/ for all sites that have an obs CSV,
#                      then restrict to site_subset (see below).
# auto_sites = false: use the sites list below verbatim (site_subset ignored).
auto_sites = true
sites      = [310]   # ignored when auto_sites = true

# When auto_sites = true, site_subset narrows the discovered sites:
#   :all           — use every discovered site
#   [304, 306,310] — use only these IDs (intersected with discovered sites)
#   20             — randomly sample this many discovered sites (reproducible
#                     via site_subset_seed) — handy for a quick multi-site
#                     test without waiting on the full set
site_subset      = 16 #:all
site_subset_seed = 1234

# ── Simulation period ─────────────────────────────────────────────────────────
# Used when auto_date_range = false.
sim_start = Date(2013, 1, 1)
sim_end   = Date(2013, 12, 31)

# auto_date_range = true: detect sim_start/sim_end from the obs file date span
# (intersected with the chosen weather source's availability — see
# WEATHER_SOURCE_START in config.jl — through yesterday).
auto_date_range = true
max_sim_years   = 5    # cap simulation to this many most-recent years of obs data

# Plot window — set either to nothing to show the full simulation period.
plot_start = nothing
plot_end   = nothing

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN LOOP                                                               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

all_stats = DataFrame()
meta_all = DataFrame(CSV.File(joinpath(@__DIR__, "Map metadata export.csv")))
meta_all[!, :ID] = strip.(string.(meta_all[!, :ID]))

# ── Auto-populate sites if requested ─────────────────────────────────────────
if auto_sites
    obs_dir = joinpath(@__DIR__, "observations")
    sites = Int[]

    if isdir(obs_dir)
        for f in readdir(obs_dir)
            endswith(f, ".csv") || continue

            sid = tryparse(Int, first(splitext(f)))
            isnothing(sid) && continue

            push!(sites, sid)
        end
        sort!(sites)
    end

    n_found = length(sites)
    println("auto_sites: found $n_found sites with obs data.")

    if site_subset isa Integer
        Random.seed!(site_subset_seed)
        sites = sort(shuffle(sites)[1:min(site_subset, n_found)])
        println("site_subset: randomly sampled $(length(sites)) of $n_found (seed=$site_subset_seed).")
    elseif site_subset != :all
        sites = sort(intersect(sites, site_subset))
        println("site_subset: restricted to $(length(sites)) of $n_found requested site(s).")
    end
end

@isdefined(weather_cache) || (global weather_cache = Dict{Tuple{Int,Date,Date}, Any}())

micro_model, soil_profile, mapper_model = build_micro_model()

make_plots_flag = save_outputs || length(sites) == 1
display_plots_flag = length(sites) == 1

# ── Phase 1: batch-fetch weather forcing for all sites at once ────────────────
prefetch_weather_batch!(weather_cache, meta_all, sites, mapper_model, soil_profile;
    sim_start, sim_end, auto_date_range, max_sim_years)

# ── Phase 2: prepare every site (metadata, weather from cache, obs, initial
#    conditions, build the MicroProblem) — sequential, but cheap; the two
#    expensive phases (Julia solve, NicheMapR) are batched below.
preps = []
for site_num in sites
    try
        prep = prepare_site(site_num, meta_all, weather_cache, micro_model, soil_profile;
            sim_start, sim_end, auto_date_range, max_sim_years)
        isnothing(prep) || push!(preps, prep)
    catch _err
        @warn "Site $site_num failed to prepare: $(_err)"
    end
end
println("\n$(length(preps)) of $(length(sites)) site(s) prepared and ready to run.")

# ── Phase 3: run all Julia models (parallel across Threads.nthreads()) and all
#    NicheMapR runs (concurrent Rscript processes, capped) ────────────────────
julia_results = run_julia_batch!(preps)
run_nmr_batch!(preps)

# ── Phase 4: stats + plots per site ───────────────────────────────────────────
for prep in preps
    try
        micro_out, julia_solve_time = julia_results[prep.site_num]
        site_rows = report_site_results(prep, micro_out, julia_solve_time;
            plot_start, plot_end, make_plots = make_plots_flag, display_plots = display_plots_flag)
        isnothing(site_rows) || append!(all_stats, DataFrame(site_rows))
    catch _err
        @warn "Site $(prep.site_num) failed to report: $(_err)"
    end
end

# ── Cross-site summary (multi-site runs) ──────────────────────────────────────
if length(sites) > 1 && nrow(all_stats) > 0
    println("\n" * "="^72)
    println("SUMMARY — $(length(sites)) sites — weather source: $(nameof(weather_source_choice))")
    println("="^72)
    println(rpad("Variable", 24) * "  " *
            "Julia: mean r  mean RMSE  mean bias     " *
            "NMR:  mean r  mean RMSE  mean bias")
    println("-"^100)
    for var in unique(all_stats.variable)
        sub = filter(r -> r.variable == var, all_stats)
        fmtv(v) = all(isnan, v) ? "  —  " : @sprintf("%.4g", mean(filter(isfinite, v)))
        jr  = all(isnan, sub.j_r)     ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.j_r)))
        jr2 = fmtv(sub.j_rmse)
        jb  = all(isnan, sub.j_bias)  ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.j_bias)))
        nr_ = all(isnan, sub.nmr_r)   ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.nmr_r)))
        nr2 = fmtv(sub.nmr_rmse)
        nb  = all(isnan, sub.nmr_bias) ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.nmr_bias)))
        println("  $(rpad(var,24))  $jr    $jr2    $jb      $nr_    $nr2    $nb")
    end
end
