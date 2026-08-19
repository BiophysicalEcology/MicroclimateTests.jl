# comparison.jl — ClimaLand Soil.EnergyHydrology vs Microclimate.jl (driven
# via MicroclimateMapper.jl/GRIDMET) vs real SCAN/SNOTEL observations.
# Runs all sites, one site, or a random subset -- same selection pattern as
# comparisons/scan_snotel/comparison.jl (see the "Sites" block below).
#
# Run once `setup.jl` has installed ClimaLand globally (see that file), with:
#   cd c:/git/BiophysicalEcologyEnv
#   julia --project=. c:/git/MicroclimateTests.jl/ClimaLand/comparison.jl
#
# ClimaLand's own solve step runs sequentially across sites (not
# Threads.@threads).

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots, Serialization
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Random
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaLand.Soil
using ClimaLand.Domains: Column
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaDiagnostics
using Interpolations
import FluidProperties

gr()

# ── Reuse comparisons/scan_snotel's machinery instead of duplicating it.
const SNOTEL_DIR = joinpath(@__DIR__, "..", "comparisons", "scan_snotel")
include(joinpath(SNOTEL_DIR, "utils.jl"))
include(joinpath(SNOTEL_DIR, "config.jl"))
include(joinpath(SNOTEL_DIR, "pipeline.jl"))

# ── ClimaLand-specific configuration + code, local to this folder.
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "forcing_climaland.jl"))
include(joinpath(@__DIR__, "run_climaland.jl"))

outputs_dir = joinpath(@__DIR__, "outputs")
mkpath(outputs_dir)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SITES — same selection pattern as scan_snotel/comparison.jl             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

auto_sites = true
sites      = [site_num]   # site_num comes from config.jl; ignored when auto_sites = true

site_subset      = 30 #:all
site_subset_seed = 1234

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PER-SITE REPORTING                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

function report_climaland_results(prep, micro_out, julia_solve_time, cl_result;
                                   make_plots = true, display_plots = true)
    s_num = prep.site_num

    depths_mc_m = ustrip.(u"m", depths)                # 19-node Microclimate depths
    depths_cl_m = cl_result.depths_cl                   # 40-node ClimaLand depths (shallow-first)
    mc_idx(target_cm) = argmin(abs.(depths_mc_m .* 100 .- target_cm))
    cl_idx(target_cm) = argmin(abs.(depths_cl_m .* 100 .- target_cm))

    depth_labels = [("5cm", 5), ("20cm", 20), ("50cm", 50), ("100cm", 100)]

    t_model_full = [DateTime(d) + Hour(h) for d in prep.dates_vec for h in 0:23]
    t_cl_full    = cl_result.start_dt .+ Second.(round.(Int, cl_result.t))

    eval_start_dt = DateTime(sim_start)
    eval_end_dt   = DateTime(sim_end, Time(23, 0, 0))
    eval_idx_mc = findall(t -> eval_start_dt <= t <= eval_end_dt, t_model_full)
    eval_idx_cl = findall(t -> eval_start_dt <= t <= eval_end_dt, t_cl_full)
    t_model = t_model_full[eval_idx_mc]
    t_cl    = t_cl_full[eval_idx_cl]
    obs     = prep.obs[(prep.obs.DateTime .>= eval_start_dt) .& (prep.obs.DateTime .<= eval_end_dt), :]
    theta_i_cl = cl_result.theta_i[eval_idx_cl, :]

    julia_cols = Dict{Symbol, Vector{Float64}}()
    cl_cols    = Dict{Symbol, Vector{Float64}}()
    for (label, cm) in depth_labels
        mi, ci = mc_idx(cm), cl_idx(cm)
        julia_cols[Symbol("D$label")]  = ustrip.(u"°C", micro_out.soil_temperature[eval_idx_mc, mi])
        julia_cols[Symbol("WC$label")] = micro_out.soil_moisture[eval_idx_mc, mi] .* 100.0
        cl_cols[Symbol("D$label")]  = cl_result.T[eval_idx_cl, ci] .- 273.15
        cl_cols[Symbol("WC$label")] = cl_result.theta_l[eval_idx_cl, ci] .* 100.0
    end
    julia_cols[:snow] = ustrip.(u"cm", micro_out.snow_depth[eval_idx_mc])
    julia_cols[:swe]  = julia_cols[:snow] .* ustrip.(u"g/cm^3", micro_out.snow_density[eval_idx_mc])

    jal = align_model_to_obs(obs.DateTime, t_model, julia_cols)
    cal = align_model_to_obs(obs.DateTime, t_cl,    cl_cols)

    obs_col_map = [
        ("Soil temp 5 cm (°C)",    :D5cm,   :STO_5cm),
        ("Soil temp 20 cm (°C)",   :D20cm,  :STO_20cm),
        ("Soil temp 50 cm (°C)",   :D50cm,  :STO_50cm),
        ("Soil temp 100 cm (°C)",  :D100cm, :STO_100cm),
        ("Soil moist 5 cm (%)",    :WC5cm,  :SMS_5cm),
        ("Soil moist 20 cm (%)",   :WC20cm, :SMS_20cm),
        ("Soil moist 50 cm (%)",   :WC50cm, :SMS_50cm),
        ("Soil moist 100 cm (%)",  :WC100cm,:SMS_100cm),
        ("Snow depth (cm)",        :snow,   :SNWD_cm),
        ("SWE (cm)",               :swe,    :WTEQ_cm),
    ]

    println("\n── Model vs observations (site $s_num, $sim_start to $sim_end) ──")
    hdr = rpad("Variable", 24) * "  " * rpad("ClimaLand", 54) * "  " * rpad("Microclimate.jl", 54)
    println(hdr); println("-"^length(hdr))

    site_rows = []
    for (label, sym, obscol) in obs_col_map
        ov = obscol in propertynames(obs) ? obs[!, obscol] : fill(missing, nrow(obs))
        s_cl = haskey(cal, sym) ? compute_stats(ov, cal[sym]) : ModelStats(NaN, NaN, NaN, 0)
        s_j  = compute_stats(ov, jal[sym])
        println(rpad(label, 24) * "  " * rpad(fmt_stat(s_cl), 54) * "  " * rpad(fmt_stat(s_j), 54))
        push!(site_rows, (site = s_num, variable = label,
            cl_r = s_cl.r, cl_rmse = s_cl.rmse, cl_bias = s_cl.bias, cl_n = s_cl.n,
            j_r  = s_j.r,  j_rmse  = s_j.rmse,  j_bias  = s_j.bias,  j_n  = s_j.n))
    end

    println("\n── ClimaLand vs Microclimate.jl (hourly, direct) ──")
    n_common = min(length(t_model), length(t_cl))
    for (label, cm) in depth_labels
        st = compute_stats(julia_cols[Symbol("D$label")][1:n_common], cl_cols[Symbol("D$label")][1:n_common])
        sm = compute_stats(julia_cols[Symbol("WC$label")][1:n_common], cl_cols[Symbol("WC$label")][1:n_common])
        println("  Soil temp $label   " * fmt_stat(st))
        println("  Soil moist $label  " * fmt_stat(sm))
    end

    CSV.write(joinpath(outputs_dir, "stats_$(s_num)_$(sim_start)_$(sim_end).csv"), DataFrame(site_rows))

    speedup = cl_result.cl_solve_time / julia_solve_time
    faster  = speedup >= 1 ? "Microclimate.jl" : "ClimaLand"
    println("\n── Solver runtime (site $s_num, $(prep.ndays) days incl. $spinup_years yr spin-up) ──")
    @printf("  Microclimate.jl : %8.3f s  (%.4f s/day)\n", julia_solve_time, julia_solve_time / prep.ndays)
    @printf("  ClimaLand       : %8.3f s  (%.4f s/day)\n", cl_result.cl_solve_time, cl_result.cl_solve_time / prep.ndays)
    @printf("  %s is %.1fx faster on this run\n", faster, max(speedup, 1 / speedup))

    timing_row = DataFrame(
        site              = [s_num],
        ndays             = [prep.ndays],
        microclimate_s    = [julia_solve_time],
        climaland_s       = [cl_result.cl_solve_time],
        microclimate_s_per_day = [julia_solve_time / prep.ndays],
        climaland_s_per_day    = [cl_result.cl_solve_time / prep.ndays],
        speedup_climaland_over_microclimate = [speedup],
        run_at            = [string(now())],
    )
    timing_file = joinpath(outputs_dir, "timings.csv")
    CSV.write(timing_file, timing_row; append = isfile(timing_file))

    if make_plots
        ps = isnothing(plot_start) ? DateTime(sim_start) : DateTime(plot_start)
        pe = isnothing(plot_end)   ? DateTime(sim_end)   : DateTime(plot_end)
        hm  = findall(t -> ps <= t <= pe, t_model)
        hmc = findall(t -> ps <= t <= pe, t_cl)
        om  = (obs.DateTime .>= ps) .& (obs.DateTime .<= pe)

        obs_vals(col) = Symbol(col) in propertynames(obs) ?
            Float64.(coalesce.(obs[om, Symbol(col)], NaN)) : Float64[]

        fig = plot(layout = (3, 3), size = (1500, 1100), dpi = 120,
                   left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)

        titles = ["Soil T 5 cm", "Soil T 20 cm", "Soil T 50 cm",
                  "Soil moisture 5 cm", "Soil moisture 20 cm", "Soil moisture 50 cm",
                  "Snow depth", "ClimaLand ice content 5 cm", "ClimaLand vs Micro T 5 cm"]
        for i in 1:9
            plot!(fig[i], title = titles[i])
        end

        for (i, (label, cm, obscol)) in enumerate([("5cm", 5, "STO_5cm"), ("20cm", 20, "STO_20cm"), ("50cm", 50, "STO_50cm")])
            plot!(fig[i], t_cl[hmc], cl_cols[Symbol("D$label")][hmc]; color = :steelblue, lw = 1.2,
                  label = i == 1 ? "ClimaLand" : "", ylabel = i == 1 ? "T (°C)" : "")
            plot!(fig[i], t_model[hm], julia_cols[Symbol("D$label")][hm]; color = :tomato, lw = 1.2,
                  label = i == 1 ? "Microclimate.jl" : "")
            !isempty(obs_vals(obscol)) &&
                scatter!(fig[i], obs.DateTime[om], obs_vals(obscol);
                         color = :green, ms = 1, alpha = 0.7, label = i == 1 ? "obs" : "")
        end

        for (j, (label, cm, obscol)) in enumerate([("5cm", 5, "SMS_5cm"), ("20cm", 20, "SMS_20cm"), ("50cm", 50, "SMS_50cm")])
            i = 3 + j
            plot!(fig[i], t_cl[hmc], cl_cols[Symbol("WC$label")][hmc]; color = :steelblue, lw = 1.2,
                  label = "", ylabel = j == 1 ? "WC (%)" : "")
            plot!(fig[i], t_model[hm], julia_cols[Symbol("WC$label")][hm]; color = :tomato, lw = 1.2, label = "")
            !isempty(obs_vals(obscol)) &&
                scatter!(fig[i], obs.DateTime[om], obs_vals(obscol);
                         color = :green, ms = 1, alpha = 0.7, label = "")
        end

        plot!(fig[7], t_model[hm], julia_cols[:snow][hm]; color = :tomato, lw = 1.2, label = "Microclimate.jl", ylabel = "cm")
        !isempty(obs_vals("SNWD_cm")) &&
            scatter!(fig[7], obs.DateTime[om], obs_vals("SNWD_cm"); color = :green, ms = 1, alpha = 0.7, label = "obs")

        ci_5cm = cl_idx(5)
        plot!(fig[8], t_cl[hmc], theta_i_cl[hmc, ci_5cm]; color = :steelblue, lw = 1.2, label = "", ylabel = "θi (m³/m³)")

        n_sc = min(length(julia_cols[:D5cm]), length(cl_cols[:D5cm]))
        lims = extrema(vcat(julia_cols[:D5cm][1:n_sc], cl_cols[:D5cm][1:n_sc]))
        scatter!(fig[9], julia_cols[:D5cm][1:n_sc], cl_cols[:D5cm][1:n_sc];
                 color = :steelblue, ms = 1.5, alpha = 0.4,
                 xlabel = "Microclimate T (°C)", ylabel = "ClimaLand T (°C)", label = "")
        plot!(fig[9], collect(lims), collect(lims); color = :black, lw = 1, ls = :dash, label = "1:1")

        plot!(fig, plot_title = "ClimaLand ($(round(cl_result.cl_solve_time, digits=1))s) vs " *
            "Microclimate.jl ($(round(julia_solve_time, digits=1))s) — SNOTEL $s_num, $sim_start to $sim_end")

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

# ── Phase 3: run Microclimate.jl (parallel across Threads.nthreads()) ────────
julia_results = run_julia_batch!(preps)

# ── Phase 4: run ClimaLand per site (sequential -- see header) + report ──────
for prep in preps
    try
        micro_out, julia_solve_time = julia_results[prep.site_num]
        cl_forcing = build_climaland_forcing(micro_out, prep)
        cl_result  = run_climaland(cl_forcing, prep)
        site_rows = report_climaland_results(prep, micro_out, julia_solve_time, cl_result;
            make_plots = make_plots_flag, display_plots = display_plots_flag)
        isnothing(site_rows) || append!(all_stats, DataFrame(site_rows))
    catch _err
        @warn "Site $(prep.site_num) failed to run/report: $(_err)"
    end
end

println("\nStats saved per-site to $(joinpath(outputs_dir, "stats_<site>_$(sim_start)_$(sim_end).csv"))")

# ── Cross-site summary (multi-site runs) ──────────────────────────────────────
if length(sites) > 1 && nrow(all_stats) > 0
    println("\n" * "="^72)
    println("SUMMARY — $(length(sites)) sites")
    println("="^72)
    println(rpad("Variable", 24) * "  " *
            "ClimaLand: mean r  mean RMSE  mean bias     " *
            "Microclimate.jl: mean r  mean RMSE  mean bias")
    println("-"^100)
    for var in unique(all_stats.variable)
        sub = filter(r -> r.variable == var, all_stats)
        fmtv(v) = all(isnan, v) ? "  —  " : @sprintf("%.4g", mean(filter(isfinite, v)))
        clr  = all(isnan, sub.cl_r)    ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.cl_r)))
        clr2 = fmtv(sub.cl_rmse)
        clb  = all(isnan, sub.cl_bias) ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.cl_bias)))
        jr   = all(isnan, sub.j_r)     ? "  —  " : @sprintf("%+.3f", mean(filter(isfinite, sub.j_r)))
        jr2  = fmtv(sub.j_rmse)
        jb   = all(isnan, sub.j_bias)  ? "  —  " : @sprintf("%+.4f", mean(filter(isfinite, sub.j_bias)))
        println("  $(rpad(var,24))  $clr    $clr2    $clb      $jr    $jr2    $jb")
    end
end
