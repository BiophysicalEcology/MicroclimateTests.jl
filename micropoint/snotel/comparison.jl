# comparison.jl — Microclimate.jl (via MicroclimateMapper.jl/GRIDMET) vs
# micropoint (R, bare-ground point mode) vs real SCAN/SNOTEL observations.
# Runs all sites, one site, or a random subset -- same selection pattern as
# comparisons/scan_snotel/comparison.jl (see the "Sites" block below).
#
# Run once install_micropoint.R has installed micropoint, with:
#   cd c:/git/BiophysicalEcologyEnv
#   julia --project=. c:/git/MicroclimateTests.jl/micropoint/comparison.jl
#
# Covers soil temperature (5/20/50cm + surface), soil moisture (bulk 0-10cm
# vs Microclimate's 5cm node -- a model-capability mismatch, not a bug, see
# README.md), downward shortwave/longwave radiation, and air temperature at
# the lower profile height. No snow: micropoint has no snow submodel (unlike
# microclimf), so snow/SWE are left out of this comparison.

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots, Serialization
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Random

gr()

# ── Reuse comparisons/scan_snotel's machinery, same as ClimaLand/comparison.jl
const SNOTEL_DIR = joinpath(@__DIR__, "..", "comparisons", "scan_snotel")
include(joinpath(SNOTEL_DIR, "utils.jl"))
include(joinpath(SNOTEL_DIR, "config.jl"))
include(joinpath(SNOTEL_DIR, "pipeline.jl"))

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "write_micropoint_inputs.jl"))

outputs_dir = joinpath(@__DIR__, "outputs")
mkpath(outputs_dir)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SITES — same selection pattern as scan_snotel/comparison.jl             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# auto_sites = true : scan scan_snotel/observations/ for all sites with an obs
#                      CSV, then restrict to site_subset (see below).
# auto_sites = false: use the sites list below verbatim (site_subset ignored).
auto_sites = true
sites      = [site_num]   # site_num comes from config.jl; ignored when auto_sites = true

# When auto_sites = true, site_subset narrows the discovered sites:
#   :all           — use every discovered site
#   [304, 306,310] — use only these IDs (intersected with discovered sites)
#   20             — randomly sample this many discovered sites (reproducible
#                     via site_subset_seed)
site_subset      = 30 #:all
site_subset_seed = 1234

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PER-SITE PIPELINE FUNCTIONS                                             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── Write micropoint inputs for many prepared sites, run R concurrently ──────
# Mirrors scan_snotel/pipeline.jl's run_nmr_batch! -- same capped-concurrent
# Rscript pattern, since micropoint is also a per-site external R process.
function run_micropoint_batch!(preps, julia_results; max_concurrent = Sys.CPU_THREADS)
    rscript = joinpath(@__DIR__, "run_micropoint.R")
    to_run = Tuple{Int,String}[]
    for prep in preps
        micro_out, _ = julia_results[prep.site_num]
        outdir = mcf_outdir(prep.site_num)
        write_micropoint_inputs(micro_out, prep, outdir)
        push!(to_run, (prep.site_num, outdir))
    end

    println("\nLaunching micropoint runs for $(length(to_run)) site(s), up to $max_concurrent concurrent...")
    queue = copy(to_run)
    running = Dict{Int, Base.Process}()
    print_lock = ReentrantLock()

    batch_time = @elapsed begin
        while !isempty(queue) || !isempty(running)
            while !isempty(queue) && length(running) < max_concurrent
                site_num_i, outdir = popfirst!(queue)
                lock(print_lock) do
                    println("  Launching micropoint for site $site_num_i...")
                end
                running[site_num_i] = run(`Rscript $rscript $outdir`; wait = false)
            end
            sleep(0.5)
            for (site_num_i, proc) in collect(running)
                if process_exited(proc)
                    lock(print_lock) do
                        status = success(proc) ? "ok" : "FAILED (exit $(proc.exitcode))"
                        println("  micropoint finished for site $site_num_i ($status)")
                    end
                    delete!(running, site_num_i)
                end
            end
        end
    end
    @printf("micropoint batch total: %.2f s (%d site(s))\n", batch_time, length(to_run))
    return nothing
end

# ── Stats + timing + plot for one prepared site ───────────────────────────────
function report_micropoint_results(prep, micro_out, julia_solve_time;
                                    make_plots = true, display_plots = true)
    s_num = prep.site_num
    outdir = mcf_outdir(s_num)

    mcf_out = DataFrame(CSV.File(joinpath(outdir, "micropoint_out.csv"); missingstring = "NA"))
    mcf_out.obs_time = DateTime.(mcf_out.obs_time, dateformat"yyyy-mm-dd HH:MM:SS")
    # write_micropoint_inputs.jl shifted obs_time from Microclimate.jl's local
    # solar time to true UTC (micropoint's required convention) before writing
    # climdata.csv -- micropoint just echoes those labels straight back, so
    # shift back here to realign with t_model/obs (both local solar time).
    mcf_out.obs_time = mcf_out.obs_time .- Hour(round(Int, -prep.lon_dd / 15))
    for col in names(mcf_out)
        col == "obs_time" && continue
        mcf_out[!, col] = Float64.(coalesce.(mcf_out[!, col], NaN))
    end
    mcf_timing = DataFrame(CSV.File(joinpath(outdir, "micropoint_timing.csv")))
    mcf_solve_time = mcf_timing.elapsed_s[1]

    depths_mc_m = ustrip.(u"m", depths)
    mc_idx(target_cm) = argmin(abs.(depths_mc_m .* 100 .- target_cm))

    t_model_full = [DateTime(d) + Hour(h) for d in prep.dates_vec for h in 0:23]

    # Discard the spin-up period: both models ran run_start..sim_end, but
    # stats/plots below cover only sim_start..sim_end.
    eval_start_dt = DateTime(sim_start)
    eval_end_dt   = DateTime(sim_end, Time(23, 0, 0))
    eval_idx_mc = findall(t -> eval_start_dt <= t <= eval_end_dt, t_model_full)
    eval_idx_r  = findall(t -> eval_start_dt <= t <= eval_end_dt, mcf_out.obs_time)
    t_model = t_model_full[eval_idx_mc]
    mcf_out = mcf_out[eval_idx_r, :]
    obs     = prep.obs[(prep.obs.DateTime .>= eval_start_dt) .& (prep.obs.DateTime .<= eval_end_dt), :]

    n_common = min(length(t_model), nrow(mcf_out))

    mc_T(cm)  = ustrip.(u"°C", micro_out.soil_temperature[eval_idx_mc, mc_idx(cm)])
    mc_WC(cm) = micro_out.soil_moisture[eval_idx_mc, mc_idx(cm)] .* 100.0

    mc_cols = Dict{Symbol, Vector{Float64}}(
        :D5cm    => mc_T(5),  :D20cm => mc_T(20), :D50cm => mc_T(50),
        :WC5cm   => mc_WC(5),
        :Tsurf   => mc_T(0),
        :Tz_h1   => ustrip.(u"°C", micro_out.profile.air_temperature[eval_idx_mc, 1]),
        :SWdown  => ustrip.(u"W/m^2", micro_out.global_radiation[eval_idx_mc]),
        :LWdown  => 5.670374419e-8 .* ustrip.(u"K", micro_out.sky_temperature[eval_idx_mc]).^4,
    )

    mcf_cols = Dict{Symbol, Vector{Float64}}(
        :D5cm   => mcf_out.T_5cm,  :D20cm => mcf_out.T_20cm, :D50cm => mcf_out.T_50cm,
        :WC5cm  => mcf_out.soilm_bulk .* 100.0,
        :Tsurf  => mcf_out.Tsurf,
        :Tz_h1  => mcf_out.Tz_h1,
        :SWdown => mcf_out.Rdirdown .+ mcf_out.Rdifdown,
        :LWdown => mcf_out.Rlwdown,
    )

    jal = align_model_to_obs(obs.DateTime, t_model, mc_cols)
    ral = align_model_to_obs(obs.DateTime, mcf_out.obs_time, mcf_cols)

    obs_col_map = [
        ("Soil temp 5 cm (°C)",   :D5cm,   :STO_5cm),
        ("Soil temp 20 cm (°C)",  :D20cm,  :STO_20cm),
        ("Soil temp 50 cm (°C)",  :D50cm,  :STO_50cm),
        ("Soil moist 0-10cm/5cm (%)", :WC5cm, :SMS_5cm),
    ]
    direct_only_map = [
        ("Soil surface temp (°C)", :Tsurf),
        ("Air temp height 1 (°C)", :Tz_h1),
        ("SW down (W/m²)",         :SWdown),
        ("LW down (W/m²)",         :LWdown),
    ]

    println("\n── Model vs observations (site $s_num, $sim_start to $sim_end) ──")
    hdr = rpad("Variable", 28) * "  " * rpad("Microclimate.jl", 54) * "  " * rpad("micropoint", 54)
    println(hdr); println("-"^length(hdr))

    site_rows = []
    for (label, sym, obscol) in obs_col_map
        ov = obscol in propertynames(obs) ? obs[!, obscol] : fill(missing, nrow(obs))
        s_j = compute_stats(ov, jal[sym])
        s_r = compute_stats(ov, ral[sym])
        println(rpad(label, 28) * "  " * rpad(fmt_stat(s_j), 54) * "  " * rpad(fmt_stat(s_r), 54))
        push!(site_rows, (site = s_num, variable = label,
            j_r = s_j.r, j_rmse = s_j.rmse, j_bias = s_j.bias, j_n = s_j.n,
            r_r = s_r.r, r_rmse = s_r.rmse, r_bias = s_r.bias, r_n = s_r.n))
    end

    println("\n── Microclimate.jl vs micropoint (hourly, direct — no observations for these) ──")
    for (label, sym) in direct_only_map
        s = compute_stats(mc_cols[sym][1:n_common], mcf_cols[sym][1:n_common])
        println("  " * rpad(label, 26) * fmt_stat(s))
        push!(site_rows, (site = s_num, variable = label * " [direct]",
            j_r = s.r, j_rmse = s.rmse, j_bias = s.bias, j_n = s.n,
            r_r = NaN, r_rmse = NaN, r_bias = NaN, r_n = 0))
    end

    CSV.write(joinpath(outputs_dir, "stats_$(s_num)_$(sim_start)_$(sim_end).csv"), DataFrame(site_rows))

    speedup = mcf_solve_time / julia_solve_time
    faster  = speedup >= 1 ? "Microclimate.jl" : "micropoint"
    println("\n── Solver runtime (site $s_num, $(prep.ndays) days incl. $spinup_years yr spin-up) ──")
    @printf("  Microclimate.jl : %8.3f s  (%.4f s/day)\n", julia_solve_time, julia_solve_time / prep.ndays)
    @printf("  micropoint      : %8.3f s  (%.4f s/day)\n", mcf_solve_time, mcf_solve_time / prep.ndays)
    @printf("  %s is %.1fx faster on this run\n", faster, max(speedup, 1 / speedup))

    timing_row = DataFrame(
        site = [s_num], ndays = [prep.ndays],
        microclimate_s = [julia_solve_time], micropoint_s = [mcf_solve_time],
        microclimate_s_per_day = [julia_solve_time / prep.ndays],
        micropoint_s_per_day   = [mcf_solve_time / prep.ndays],
        speedup_micropoint_over_microclimate = [speedup],
        run_at = [string(now())],
    )
    timing_file = joinpath(outputs_dir, "timings.csv")
    CSV.write(timing_file, timing_row; append = isfile(timing_file))

    if make_plots
        ps = isnothing(plot_start) ? DateTime(sim_start) : DateTime(plot_start)
        pe = isnothing(plot_end)   ? DateTime(sim_end)   : DateTime(plot_end)
        hm  = findall(t -> ps <= t <= pe, t_model)
        hmr = findall(t -> ps <= t <= pe, mcf_out.obs_time)
        om  = (obs.DateTime .>= ps) .& (obs.DateTime .<= pe)

        obs_vals(col) = Symbol(col) in propertynames(obs) ?
            Float64.(coalesce.(obs[om, Symbol(col)], NaN)) : Float64[]

        fig = plot(layout = (2, 4), size = (1700, 800), dpi = 120,
                   left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)

        titles = ["Soil T 5 cm", "Soil T 20 cm", "Soil T 50 cm", "Soil moisture (0-10cm / 5cm)",
                  "Air T height 1", "Soil surface T", "SW down", "LW down"]
        for i in 1:8
            plot!(fig[i], title = titles[i])
        end

        for (i, (sym, obscol)) in enumerate([(:D5cm, "STO_5cm"), (:D20cm, "STO_20cm"), (:D50cm, "STO_50cm")])
            plot!(fig[i], t_model[hm], mc_cols[sym][hm]; color = :tomato, lw = 1.2,
                  label = i == 1 ? "Microclimate.jl" : "", ylabel = i == 1 ? "T (°C)" : "")
            plot!(fig[i], mcf_out.obs_time[hmr], mcf_cols[sym][hmr]; color = :steelblue, lw = 1.2,
                  label = i == 1 ? "micropoint" : "")
            !isempty(obs_vals(obscol)) &&
                scatter!(fig[i], obs.DateTime[om], obs_vals(obscol);
                         color = :green, ms = 1, alpha = 0.7, label = i == 1 ? "obs" : "")
        end

        plot!(fig[4], t_model[hm], mc_cols[:WC5cm][hm]; color = :tomato, lw = 1.2, label = "", ylabel = "%")
        plot!(fig[4], mcf_out.obs_time[hmr], mcf_cols[:WC5cm][hmr]; color = :steelblue, lw = 1.2, label = "")
        !isempty(obs_vals("SMS_5cm")) &&
            scatter!(fig[4], obs.DateTime[om], obs_vals("SMS_5cm"); color = :green, ms = 1, alpha = 0.7, label = "")

        for (i, sym) in [(5, :Tz_h1), (6, :Tsurf)]
            plot!(fig[i], t_model[hm], mc_cols[sym][hm]; color = :tomato, lw = 1.2, label = "", ylabel = i == 5 ? "T (°C)" : "")
            plot!(fig[i], mcf_out.obs_time[hmr], mcf_cols[sym][hmr]; color = :steelblue, lw = 1.2, label = "")
        end

        for (i, sym) in [(7, :SWdown), (8, :LWdown)]
            plot!(fig[i], t_model[hm], mc_cols[sym][hm]; color = :tomato, lw = 1.2, label = "", ylabel = i == 7 ? "W/m²" : "")
            plot!(fig[i], mcf_out.obs_time[hmr], mcf_cols[sym][hmr]; color = :steelblue, lw = 1.2, label = "")
        end

        plot!(fig, plot_title = "Microclimate.jl ($(round(julia_solve_time, digits=1))s) vs " *
            "micropoint ($(round(mcf_solve_time, digits=1))s) — SNOTEL $s_num, $sim_start to $sim_end")

        fig_path = joinpath(outputs_dir, "comparison_$(s_num)_$(sim_start)_$(sim_end).png")
        savefig(fig, fig_path)
        display_plots && display(fig)
        println("\nPlot saved to $fig_path")
    end

    return site_rows
end

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN                                                                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝

all_stats = DataFrame()
meta_all = DataFrame(CSV.File(joinpath(SNOTEL_DIR, "Map metadata export.csv")))
meta_all[!, :ID] = strip.(string.(meta_all[!, :ID]))

if auto_sites
    obs_dir = joinpath(SNOTEL_DIR, "observations")
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
    sim_start = run_start, sim_end, auto_date_range = false, max_sim_years = spinup_years + 1)

# ── Phase 2: prepare every site (metadata, weather from cache, obs, initial
#    conditions, per-site SoilGrids profile, build the MicroProblem) ─────────
preps = []
for s in sites
    try
        prep = prepare_site(s, meta_all, weather_cache, micro_model, soil_profile;
            sim_start = run_start, sim_end, auto_date_range = false, max_sim_years = spinup_years + 1)
        isnothing(prep) || push!(preps, prep)
    catch _err
        @warn "Site $s failed to prepare: $(_err)"
    end
end
println("\n$(length(preps)) of $(length(sites)) site(s) prepared and ready to run.")

# ── Phase 3: run Microclimate.jl (parallel across Threads.nthreads()) and
#    micropoint (concurrent Rscript processes, capped) ───────────────────────
julia_results = run_julia_batch!(preps)
run_micropoint_batch!(preps, julia_results)

# ── Phase 4: stats + plots per site ───────────────────────────────────────────
for prep in preps
    try
        micro_out, julia_solve_time = julia_results[prep.site_num]
        site_rows = report_micropoint_results(prep, micro_out, julia_solve_time;
            make_plots = make_plots_flag, display_plots = display_plots_flag)
        isnothing(site_rows) || append!(all_stats, DataFrame(site_rows))
    catch _err
        @warn "Site $(prep.site_num) failed to report: $(_err)"
    end
end

println("\nStats saved per-site to $(joinpath(outputs_dir, "stats_<site>_$(sim_start)_$(sim_end).csv"))")

# ── Cross-site summary (multi-site runs) ──────────────────────────────────────
if length(sites) > 1 && nrow(all_stats) > 0
    println("\n" * "="^72)
    println("SUMMARY — $(length(sites)) sites")
    println("="^72)
    println(rpad("Variable", 28) * "  " *
            "Microclimate.jl: mean r  mean RMSE  mean bias     " *
            "micropoint: mean r  mean RMSE  mean bias")
    println("-"^100)
    for var in unique(all_stats.variable)
        sub = filter(r -> r.variable == var, all_stats)
        fmtv(v) = all(isnan, v) ? "  —  " : @sprintf("%.4g", mean(filter(isfinite, v)))
        jr  = all(isnan, sub.j_r)    ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.j_r)))
        jr2 = fmtv(sub.j_rmse)
        jb  = all(isnan, sub.j_bias) ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.j_bias)))
        rr  = all(isnan, sub.r_r)    ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.r_r)))
        rr2 = fmtv(sub.r_rmse)
        rb  = all(isnan, sub.r_bias) ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.r_bias)))
        println("  $(rpad(var,28))  $jr    $jr2    $jb      $rr    $rr2    $rb")
    end
end
