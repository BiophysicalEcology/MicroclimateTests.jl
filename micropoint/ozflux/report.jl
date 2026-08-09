# report.jl — stats + plots for the micropoint (R) vs Microclimate.jl vs
# real OzFlux obs 3-way comparison, given a solved comparisons/ozflux
# `result` and a `read_micropoint_output` bundle. Requires comparisons/
# ozflux/{config,utils,pipeline,report}.jl and this directory's utils.jl
# already included (for _col, _height_field, _profile_series,
# _vertical_profile, _with_absolute_humidity, _ah_from_rh, compute_stats,
# fmt_stat, _hourly_depth_series, _interp_col).

# Single-height/flux timeseries panel, 3-way -- obs=red, Julia=black,
# micropoint=steelblue, same visual conventions as comparisons/ozflux/
# report.jl's _timeseries_panel.
function _timeseries_panel3(t_model, obs_vec, jl_vec, mp_vec; title, units, plot_start, plot_end)
    ps = isnothing(plot_start) ? t_model[1] : plot_start
    pe = isnothing(plot_end) ? t_model[end] : plot_end
    m = findall(t -> ps <= t <= pe, t_model)
    p = plot(t_model[m], Float64.(jl_vec[m]); label="Julia", color=:black, lw=1.2, title, ylabel=units, legend=:topleft)
    plot!(p, t_model[m], Float64.(mp_vec[m]); label="micropoint", color=:steelblue, lw=1)
    if obs_vec !== nothing
        obs_f = Float64[ismissing(v) ? NaN : Float64(v) for v in obs_vec]
        plot!(p, t_model[m], obs_f[m]; label="obs", color=:red, lw=1, alpha=0.7)
    end
    return p
end

function _scatter_panel3(obs_vec, model_vec; units, label, color)
    obs_vec === nothing && return nothing, plot(; title="$label\n(no obs)", titlefontsize=9, framestyle=:none)
    pairs = [(Float64(o), Float64(mv)) for (o, mv) in zip(obs_vec, model_vec) if !ismissing(o) && isfinite(Float64(o)) && isfinite(Float64(mv))]
    isempty(pairs) && return nothing, plot(; title="$label\n(no valid pairs)", titlefontsize=9, framestyle=:none)
    o = first.(pairs); mv = last.(pairs)
    lo, hi = extrema(vcat(o, mv))
    s = compute_stats(o, mv)
    p = scatter(o, mv; label=nothing, xlabel="obs ($units)", ylabel="$label ($units)",
        title=@sprintf("r=%.3f  RMSE=%.3g\nbias=%+.3g  n=%d", s.r, s.rmse, s.bias, s.n),
        titlefontsize=9, ms=2, alpha=0.35, color, xlims=(lo, hi), ylims=(lo, hi))
    plot!(p, [lo, hi], [lo, hi]; label="1:1", color=:gray, lw=1)
    return s, p
end

# One vertical-profile snapshot panel (Ta/Ws/AH vs height, at a single hour):
# Julia's own height grid (canopy-resolved below canopy_height, MOST above --
# see comparisons/ozflux/report.jl's _profile_series) as a line, micropoint's
# own evenly-spaced grid as a second line, tower obs at each configured
# sensor height as points.
function _profile_panel3(result, output, dt, field, jl_convert, mp_mat, mp_z_m, mp_convert; units)
    (; t_model, heights, resolved, hourly, height_series, canopy_mode) = result
    canopy_height = resolved.canopy_height
    title = Dates.format(dt, "yyyy-mm-dd HH:MM")
    i = findfirst(==(dt), t_model)
    if i === nothing
        return plot(; title="$title\noutside solved range", titlefontsize=7, framestyle=:none)
    end
    jl_vals = _vertical_profile(output, heights, canopy_height, i, field, canopy_mode, jl_convert)
    heights_m = ustrip.(u"m", heights)
    order = sortperm(heights_m)
    keep = [k for k in order if !ismissing(jl_vals[k])]
    p = plot(jl_vals[keep], heights_m[keep]; label="Julia", color=:black, lw=1.5,
        xlabel=units, ylabel="height (m)", title, titlefontsize=8, legend=false)
    mp_vals = mp_convert.(mp_mat[i, :])
    plot!(p, mp_vals, mp_z_m; label="micropoint", color=:steelblue, lw=1.5)
    hline!(p, [ustrip(u"m", canopy_height)]; label=nothing, color=:green, ls=:dash, lw=1)
    row = findfirst(==(dt), hourly.DateTime)
    if row !== nothing
        for (h, name) in height_series
            (_height_field(name) == field && name in names(hourly)) || continue
            v = hourly[row, name]
            ismissing(v) || scatter!(p, [Float64(v)], [h]; label=nothing, color=:red, ms=6,
                markerstrokewidth=1.5, markerstrokecolor=:white)
        end
    end
    return p
end

# Grid of _profile_panel3, one panel per profile_times entry.
function plot_profiles3(result, output, profile_times, field, jl_convert, mp_mat, mp_z_m, mp_convert;
                         units, save_dir=nothing, tag="", display_plots=false)
    panels = [_profile_panel3(result, output, dt, field, jl_convert, mp_mat, mp_z_m, mp_convert; units) for dt in profile_times]
    ncols = min(4, length(panels))
    nrows = cld(length(panels), ncols)
    p = plot(panels...; layout=(nrows, ncols), size=(320 * ncols, 260 * nrows),
        plot_title="$field profile — $(result.site_name) (Julia vs micropoint)")
    display_plots && display(p)
    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p, joinpath(save_dir, "$(tag)_$(field)_profiles.png"))
    end
    return nothing
end

# Top-level entry: stats/plots for every comparison target (radiation,
# height-series Ta/Ws, top-of-canopy leaf/air temperature, soil Ts/Sws depth
# panels, vertical profile snapshots), given a solved comparisons/ozflux
# `result` (run_site_gapfilled) and a `read_micropoint_output(outdir)` bundle.
function report_site_results3(result, mp; outdir, plot_start=nothing, plot_end=nothing, profile_times,
                               save_outputs_3way=true, display_plots_3way=false)
    (; output, t_model, resolved, hourly, depths, heights, height_series, canopy_mode, site_name, years) = result
    canopy_height = resolved.canopy_height
    canopy_height_m = ustrip(u"m", canopy_height)
    tag = "$(site_name)_$(join(years, '-'))"
    site_dir(sub) = joinpath(outdir, sub)
    mp_at_height(mat, h_m) = _interp_col(mat, mp.z_m, h_m)
    mp_at_depth(mat, d_m) = _interp_col(mat, mp.depth_m, d_m)

    all_stats = NamedTuple[]
    record!(label, kind, source, s) = push!(all_stats,
        (; site=site_name, years=join(years, '-'), variable=label, kind, source, r=s.r, rmse=s.rmse, bias=s.bias, n=s.n))

    println("\n" * "="^72)
    println("Site $site_name  years=$years — obs vs Microclimate.jl vs micropoint")
    println("="^72)

    function report3!(label, kind, key, obs_vec, jl_vec, mp_vec; units)
        sj, p_sc_j = _scatter_panel3(obs_vec, jl_vec; units, label="Julia", color=:black)
        sm, p_sc_m = _scatter_panel3(obs_vec, mp_vec; units, label="micropoint", color=:steelblue)
        sj === nothing || (println(rpad("$label (obs vs Julia)", 38) * fmt_stat(sj)); record!(label, kind, "julia", sj))
        sm === nothing || (println(rpad("$label (obs vs micropoint)", 38) * fmt_stat(sm)); record!(label, kind, "micropoint", sm))
        p_ts = _timeseries_panel3(t_model, obs_vec, jl_vec, mp_vec; title=label, units, plot_start, plot_end)
        p = plot(p_ts, p_sc_j, p_sc_m; layout=@layout([a{0.5w} b c]), size=(1400, 350))
        display_plots_3way && display(p)
        if save_outputs_3way
            d = site_dir(kind); mkpath(d)
            savefig(p, joinpath(d, "$(tag)_$(key).png"))
        end
        return nothing
    end

    # ── Fsu, Flu, Fn: canopy-top boundary ─────────────────────────────────────
    jl_fsu = ustrip.(u"W/m^2", output.canopy.boundary_upward_shortwave[:, 1])
    jl_flu = ustrip.(u"W/m^2", output.canopy.boundary_upward_longwave[:, 1])
    jl_fn = ustrip.(u"W/m^2",
        output.canopy.boundary_downward_shortwave[:, 1] .- output.canopy.boundary_upward_shortwave[:, 1] .+
        output.canopy.boundary_downward_longwave[:, 1] .- output.canopy.boundary_upward_longwave[:, 1])
    mp_fsu = mp.Rswup[:, end]
    mp_flu = mp.Rlwup[:, end]
    mp_fn = (mp.Rdirdown[:, end] .+ mp.Rdifdown[:, end]) .- mp.Rswup[:, end] .+ mp.Rlwdown[:, end] .- mp.Rlwup[:, end]
    for (label, key, jl_v, mp_v, obscol) in (
        ("Upward shortwave (Fsu)", "Fsu", jl_fsu, mp_fsu, "Fsu"),
        ("Upward longwave (Flu)", "Flu", jl_flu, mp_flu, "Flu"),
        ("Net radiation (Fn)", "Fn", jl_fn, mp_fn, "Fn"),
    )
        obs_vec = _col(hourly, t_model, obscol)
        report3!(label, "radiation", key, obs_vec, jl_v, mp_v; units="W/m^2")
    end

    # ── Ta/Ws at real sub-canopy sensor heights (below canopy_height only --
    # micropoint's output tops out at canopy_height; heights above it stay a
    # 2-way comparison, already covered by comparisons/ozflux/report.jl). ────
    for (h_m, name) in height_series
        h_m <= canopy_height_m || continue
        field = _height_field(name)
        field === nothing && continue
        obs_vec = _col(hourly, t_model, name)
        obs_vec === nothing && continue
        jl_vec = _profile_series(output, heights, canopy_height, h_m, field, canopy_mode)
        jl_vec === nothing && continue
        mp_mat = field == :air_temperature ? mp.tair : mp.wind
        mp_vec = mp_at_height(mp_mat, h_m)
        jl_conv, units = field == :air_temperature ? (ustrip.(u"°C", u"°C".(jl_vec)), "°C") : (ustrip.(u"m/s", jl_vec), "m/s")
        report3!("$name @ $(h_m) m", "profile_$(field)", name, obs_vec, jl_conv, mp_vec; units)
    end

    # ── Leaf and air temperature, top of canopy: Julia vs micropoint only, no
    # tower obs. ───────────────────────────────────────────────────────────
    jl_leaf = ustrip.(u"°C", u"°C".(output.canopy.leaf_temperature[:, 1]))
    jl_air = ustrip.(u"°C", u"°C".(output.canopy.air_temperature[:, 1]))
    mp_leaf = mp.tleaf[:, end]
    mp_air = mp.tair[:, end]
    s_leaf = compute_stats(jl_leaf, mp_leaf)
    s_air = compute_stats(jl_air, mp_air)
    println(rpad("Leaf temperature, top of canopy (Julia vs micropoint)", 55) * fmt_stat(s_leaf))
    println(rpad("Air temperature, top of canopy (Julia vs micropoint)", 55) * fmt_stat(s_air))
    record!("Leaf temperature (top)", "leaf_temperature", "micropoint", s_leaf)
    record!("Air temperature (top)", "air_temperature_top", "micropoint", s_air)

    ps = isnothing(plot_start) ? t_model[1] : plot_start
    pe = isnothing(plot_end) ? t_model[end] : plot_end
    m = findall(t -> ps <= t <= pe, t_model)
    p_leaf = plot(t_model[m], jl_leaf[m]; label="Julia leaf", color=:black, lw=1.2,
        title="Leaf and air temperature, top of canopy", ylabel="°C", legend=:topleft)
    plot!(p_leaf, t_model[m], jl_air[m]; label="Julia air", color=:gray, lw=1)
    plot!(p_leaf, t_model[m], mp_leaf[m]; label="micropoint leaf", color=:steelblue, lw=1.2)
    plot!(p_leaf, t_model[m], mp_air[m]; label="micropoint air", color=:orange, lw=1)
    if save_outputs_3way
        d = site_dir("leaf_temperature"); mkpath(d)
        savefig(p_leaf, joinpath(d, "$(tag)_leaf_temperature.png"))
    end
    display_plots_3way && display(p_leaf)

    # ══════════════════════════════════════════════════════════════════════
    # Soil temperature / moisture at every observed depth. Grid of timeseries
    # panels (one per depth, 3 lines each), mirroring comparisons/ozflux/
    # report.jl's _depth_panels! grid style.
    # ══════════════════════════════════════════════════════════════════════
    function depth_panels3!(prefix, kind, jl_mat, mp_mat, units, convert)
        depth_rows = _hourly_depth_series(hourly, prefix)
        isempty(depth_rows) && return nothing
        ts_panels, sc_j_panels, sc_m_panels = [], [], []
        for r in depth_rows
            depth_m = r.depth_m
            obs_vec = _col(hourly, t_model, "$(prefix)_$(depth_m)m")
            i = argmin(abs.(ustrip.(u"m", depths) .- depth_m))
            jl_vec = convert(jl_mat[:, i])
            mp_vec = mp_at_depth(mp_mat, depth_m)
            label = "$kind $(depth_m) m"
            sj = compute_stats(obs_vec, jl_vec)
            sm = compute_stats(obs_vec, mp_vec)
            println(rpad("$label (obs vs Julia)", 38) * fmt_stat(sj))
            println(rpad("$label (obs vs micropoint)", 38) * fmt_stat(sm))
            record!(label, kind, "julia", sj)
            record!(label, kind, "micropoint", sm)
            push!(ts_panels, _timeseries_panel3(t_model, obs_vec, jl_vec, mp_vec; title="$(depth_m) m", units, plot_start, plot_end))
            _, p_sc_j = _scatter_panel3(obs_vec, jl_vec; units, label="$(depth_m) m", color=:black)
            _, p_sc_m = _scatter_panel3(obs_vec, mp_vec; units, label="$(depth_m) m", color=:steelblue)
            push!(sc_j_panels, p_sc_j)
            push!(sc_m_panels, p_sc_m)
        end
        isempty(ts_panels) && return nothing
        ncols = min(4, length(ts_panels))
        nrows = cld(length(ts_panels), ncols)
        # Minimum sizes prevent a single-panel figure becoming only 320 × 240 px.
        ts_size = (max(900, 420 * ncols), max(500, 300 * nrows))
        scatter_size = (max(700, 360 * ncols), max(600, 360 * nrows))
        common_kwargs = (;
            layout=(nrows, ncols), dpi=150, margin=4mm, left_margin=7mm, right_margin=4mm,
            top_margin=8mm, bottom_margin=10mm, plot_titlefontsize=14, tickfontsize=8,
            guidefontsize=10, titlefontsize=10,
        )
        p_ts = plot(ts_panels...; common_kwargs..., size=ts_size, plot_title="$kind timeseries", xrotation=30)
        p_sc_j = plot(sc_j_panels...; common_kwargs..., size=scatter_size, plot_title="$kind scatter — obs vs Julia")
        p_sc_m = plot(sc_m_panels...; common_kwargs..., size=scatter_size, plot_title="$kind scatter — obs vs micropoint")
        if display_plots_3way
            display(p_ts)
            display(p_sc_j)
            display(p_sc_m)
        end
        if save_outputs_3way
            d = site_dir(kind)
            mkpath(d)
            savefig(p_ts, joinpath(d, "$(tag)_timeseries_panels.png"))
            savefig(p_sc_j, joinpath(d, "$(tag)_scatter_julia_panels.png"))
            savefig(p_sc_m, joinpath(d, "$(tag)_scatter_micropoint_panels.png"))
        end
        return nothing
    end
    depth_panels3!("Ts", "soil_temperature", output.soil_temperature, mp.tsoil, "°C", d -> ustrip.(u"°C", d))
    depth_panels3!("Sws", "soil_moisture", output.soil_moisture, mp.theta, "frac", identity)

    # ══════════════════════════════════════════════════════════════════════
    # Vertical profile snapshots (Ta, Ws, AH) at profile_times.
    # ══════════════════════════════════════════════════════════════════════
    plot_profiles3(result, output, profile_times, :air_temperature, t -> ustrip(u"°C", t), mp.tair, mp.z_m, identity;
        units="°C", save_dir=site_dir("profiles_vertical"), tag, display_plots=display_plots_3way)
    plot_profiles3(result, output, profile_times, :wind_speed, w -> ustrip(u"m/s", w), mp.wind, mp.z_m, identity;
        units="m/s", save_dir=site_dir("profiles_vertical"), tag, display_plots=display_plots_3way)
    output_ah = _with_absolute_humidity(output)
    mp_ah = _ah_from_rh.(mp.relhum ./ 100.0, mp.tair .* u"°C")
    plot_profiles3(result, output_ah, profile_times, :absolute_humidity, identity, mp_ah, mp.z_m, identity;
        units="g/m^3", save_dir=site_dir("profiles_vertical"), tag, display_plots=display_plots_3way)

    stats_df = DataFrame(all_stats)
    save_outputs_3way && CSV.write(joinpath(outdir, "$(tag)_stats.csv"), stats_df)
    return stats_df
end
