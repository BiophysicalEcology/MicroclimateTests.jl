# report.jl — stats + plots for every comparison target, given a solved
# `run_site` result. Requires config.jl/utils.jl/pipeline.jl already included.

const _FIG_SUBDIR = Dict(
    "Fsu" => "shortwave", "Flu" => "longwave", "Fn" => "net_radiation",
    "Fh" => "sensible_heat", "Fe" => "latent_heat", "Fg" => "ground_heat",
)

# Obs, reindexed onto `t_model` by DateTime (`missing` where t_model has no
# matching obs row) -- a no-op reindex when hourly.DateTime == t_model (the
# tower-forced run), but required for the SILO-forced run, whose model span
# (full calendar year) needn't match the tower file's own coverage.
function _col(hourly, t_model, name)
    name in names(hourly) || return nothing
    col = hourly[!, name]
    hourly.DateTime === t_model && return col
    idx = Dict(t => i for (i, t) in enumerate(hourly.DateTime))
    return [haskey(idx, t) ? col[idx[t]] : missing for t in t_model]
end

function _height_index(heights, h_m; atol=0.05)
    findfirst(x -> abs(ustrip(u"m", x) - h_m) <= atol, heights)
end

# RH (0-1) + air temperature -> absolute humidity (g/m^3) -- inverse of
# pipeline.jl's _rh_from_ah. Absolute humidity isolates moisture transport
# from the temperature profile's own shape (unlike RH = e/e_sat(T), which
# conflates the two), so it's the more diagnostic quantity for judging
# canopy air-profile physics rather than RH directly.
function _ah_from_rh(rh_frac, T, vapour_pressure_equation=GoffGratch())
    T_K = u"K"(T)
    p_v = rh_frac * FluidProperties.vapour_pressure(vapour_pressure_equation, T_K)
    return ustrip(u"g/m^3", p_v * _MW_WATER / (Unitful.R * T_K))
end

# output.canopy/output.profile as plain NamedTuples with an added
# absolute_humidity field (derived elementwise from relative_humidity +
# air_temperature) -- lets plot_canopy_profiles's existing field=Symbol
# machinery plot it without any changes there.
_as_namedtuple(x) = NamedTuple{fieldnames(typeof(x))}(getfield.(Ref(x), fieldnames(typeof(x))))
function _with_absolute_humidity(output)
    add_ah(sub) = merge(_as_namedtuple(sub),
        (; absolute_humidity = _ah_from_rh.(sub.relative_humidity, sub.air_temperature)))
    return merge(_as_namedtuple(output), (; canopy = add_ah(output.canopy), profile = add_ah(output.profile)))
end

# Above-canopy heights use the canopy-blind MOST profile (output.profile.*)
# and match exactly (that grid does include the exact tower/sensor heights,
# see _build_heights). Below-canopy heights use the canopy-resolved output
# (output.canopy.*) but the physics grid there is deliberately kept uniform
# (not sensor-height-injected, see _build_heights), so this interpolates
# linearly between the two nearest canopy layers bracketing h_m instead of
# requiring an exact match. NoCanopy (canopy_mode=:legacy) has no
# below-canopy resolution at all (output.canopy.* is zero-column) -- that's
# the whole point of comparing it.
#
# output.canopy's columns are top-to-bottom (canopy_layer_heights sorts
# descending, matching leaf_temperature's canopy-top-is-index-1 convention),
# while `heights`/below_idx are bottom-to-top (ascending) -- _canopy_col
# maps a position in the ascending below-canopy block to the matching
# (reversed) output.canopy column.
_canopy_col(j, n_canopy) = n_canopy - j + 1

function _profile_series(output, heights, canopy_height, h_m, field, canopy_mode)
    if h_m * u"m" > canopy_height
        k = _height_index(heights, h_m)
        k === nothing && return nothing
        return getproperty(output.profile, field)[:, k]
    end
    canopy_mode == :legacy && return nothing
    below_idx = findall(h -> h <= canopy_height, heights)
    isempty(below_idx) && return nothing
    n_canopy = length(below_idx)
    below_heights_m = ustrip.(u"m", heights[below_idx])
    mat = getproperty(output.canopy, field)
    j = searchsortedlast(below_heights_m, h_m)
    if j <= 0
        return mat[:, _canopy_col(1, n_canopy)]
    elseif j >= n_canopy
        return mat[:, _canopy_col(n_canopy, n_canopy)]
    else
        lo, hi = below_heights_m[j], below_heights_m[j + 1]
        t = (h_m - lo) / (hi - lo)
        return mat[:, _canopy_col(j, n_canopy)] .* (1 - t) .+ mat[:, _canopy_col(j + 1, n_canopy)] .* t
    end
end

_height_field(name) = startswith(name, "Ta") ? :air_temperature :
                       startswith(name, "Ws") ? :wind_speed : nothing  # AH needs conversion to RH -- not implemented, see README

# Model value at every height, for one hour (`i` into t_model) -- the
# vertical counterpart of _profile_series (which slices one height across
# every hour). `missing` for below-canopy heights under canopy_mode=:legacy
# (no canopy-resolved output there, same guard as _profile_series). See
# _canopy_col for why below-canopy heights need a reversed column index.
function _vertical_profile(output, heights, canopy_height, i, field, canopy_mode, convert)
    n_canopy = count(h -> h <= canopy_height, heights)
    vals = Vector{Union{Float64,Missing}}(missing, length(heights))
    for k in eachindex(heights)
        below_canopy = heights[k] <= canopy_height
        below_canopy && canopy_mode == :legacy && continue
        if below_canopy
            vals[k] = convert(getproperty(output.canopy, field)[i, _canopy_col(k, n_canopy)])
        else
            vals[k] = convert(getproperty(output.profile, field)[i, k])
        end
    end
    return vals
end

# One vertical profile (value vs height) at a single hour: the model's full
# height-grid profile (canopy-resolved below canopy_height, free-atmosphere
# MOST above -- see _profile_series) as a line, tower obs at each configured
# sensor height (SITE_HEIGHT_SERIES) as points. The sharpest available test
# of whether the *shape* of the profile is right, not just an aggregate flux.
function _profile_panel(prep, output, dt; field, units, convert, xlims=:auto)
    (; heights, resolved, hourly, height_series, canopy_mode, t_model) = prep
    canopy_height = resolved.canopy_height
    title = Dates.format(dt, "yyyy-mm-dd HH:MM")
    i = findfirst(==(dt), t_model)
    if i === nothing
        range_str = "$(first(t_model)) to $(last(t_model))"
        return plot(; title="$title\noutside solved range\n($range_str)", titlefontsize=7, framestyle=:none)
    end

    model_vals = _vertical_profile(output, heights, canopy_height, i, field, canopy_mode, convert)
    heights_m = ustrip.(u"m", heights)
    order = sortperm(heights_m)
    keep = [k for k in order if !ismissing(model_vals[k])]
    p = plot(model_vals[keep], heights_m[keep]; label="model", color=:black, lw=1.5,
        xlabel=units, ylabel="height (m)", title, titlefontsize=8, legend=false, xlims)
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

# Common x-axis range across every panel in one plot_canopy_profiles grid
# (model line + any obs scatter) -- panels default to independent auto-ranges
# otherwise, making the shape across hours hard to compare by eye.
function _shared_profile_xlims(prep, output, profile_times; field, convert)
    (; heights, resolved, hourly, height_series, canopy_mode, t_model) = prep
    canopy_height = resolved.canopy_height
    vals = Float64[]
    for dt in profile_times
        i = findfirst(==(dt), t_model)
        i === nothing && continue
        append!(vals, skipmissing(_vertical_profile(output, heights, canopy_height, i, field, canopy_mode, convert)))
        row = findfirst(==(dt), hourly.DateTime)
        row === nothing && continue
        for (h, name) in height_series
            (_height_field(name) == field && name in names(hourly)) || continue
            v = hourly[row, name]
            ismissing(v) || push!(vals, Float64(v))
        end
    end
    isempty(vals) && return :auto
    lo, hi = extrema(vals)
    pad = 0.05 * max(hi - lo, 1.0e-6)
    return (lo - pad, hi + pad)
end

# Grid of _profile_panel (as square a grid as fits, not one long column --
# stacking every hour into a single column makes an unreadably tall, squashed
# image), one panel per requested DateTime.
function plot_canopy_profiles(prep, output, profile_times; field, units, convert,
                               save_dir=nothing, tag="", display_plots=false)
    xlims = _shared_profile_xlims(prep, output, profile_times; field, convert)
    panels = [_profile_panel(prep, output, dt; field, units, convert, xlims) for dt in profile_times]
    ncols = min(4, length(panels))
    nrows = cld(length(panels), ncols)
    p = plot(panels...; layout=(nrows, ncols), size=(320 * ncols, 260 * nrows),
        plot_title="$field profile — $(prep.site_name)")
    display_plots && display(p)
    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p, joinpath(save_dir, "$(tag)_$(field)_profiles.$figure_format"))
    end
    return p
end

# Ta + Ws + RH profile grids for every hour in `profile_times` -- call this
# directly from debug_site.jl. Most meaningful at sites with real sub-canopy
# sensors (Whroo, Wallaby); harmless (just model-only panels, no obs points)
# at sites without SITE_HEIGHT_SERIES entries -- true for every site for RH
# specifically, since no site has a height-resolved RH/AH series wired into
# _height_field (AH->RH conversion not implemented, see README), so the RH
# panels are always model-only.
function plot_canopy_profiles_all(prep, output, profile_times; save_to_disk=save_outputs, display_plots=true)
    tag = "$(prep.site_name)_$(join(prep.years, '-'))"
    save_dir = save_to_disk ? joinpath(outputs_dir, "profiles_vertical") : nothing
    plot_canopy_profiles(prep, output, profile_times; field=:air_temperature, units="°C",
        convert=t -> ustrip(u"°C", t), save_dir, tag, display_plots)
    plot_canopy_profiles(prep, output, profile_times; field=:wind_speed, units="m/s",
        convert=w -> ustrip(u"m/s", w), save_dir, tag, display_plots)
    output_ah = _with_absolute_humidity(output)
    plot_canopy_profiles(prep, output_ah, profile_times; field=:absolute_humidity, units="g/m^3",
        convert=identity, save_dir, tag, display_plots)
    return nothing
end

# Top-of-canopy leaf temperature vs the local (same-layer) air temperature --
# direct diagnostic for how far leaves run from air, independent of any tower
# obs. Useful when tuning canopy_relaxation_choice/canopy_convergence_choice:
# a converged equilibrium that's still far from air temperature is a real
# physics/parameter issue, not a convergence one. NoCanopy (canopy_mode=
# :legacy) has no leaf temperature at all -- skipped there.
function _plot_leaf_vs_air(heights, canopy_height, canopy_mode, t_model, output;
                            plot_start=nothing, plot_end=nothing, save_dir=nothing, tag="", display_plots=false)
    canopy_mode == :legacy && return nothing
    count(h -> h <= canopy_height, heights) == 0 && return nothing
    # output.canopy's columns are top-to-bottom -- column 1 is canopy top.
    leaf_c = ustrip.(u"°C", u"°C".(output.canopy.leaf_temperature[:, 1]))
    air_c  = ustrip.(u"°C", u"°C".(output.canopy.air_temperature[:, 1]))

    ps = isnothing(plot_start) ? t_model[1] : plot_start
    pe = isnothing(plot_end) ? t_model[end] : plot_end
    m = findall(t -> ps <= t <= pe, t_model)
    p = plot(t_model[m], leaf_c[m]; label="leaf (top of canopy)", color=:darkorange, lw=1.2,
        title="Leaf vs air temperature — top of canopy", ylabel="°C")
    plot!(p, t_model[m], air_c[m]; label="air (same layer)", color=:black, lw=1)
    display_plots && display(p)
    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p, joinpath(save_dir, "$(tag).$figure_format"))
    end
    return nothing
end

# One combined figure per variable: time-varying obs-vs-model trace on the
# left, scatter/1:1/correlation panel on the right (not two separate files).
const _PAIR_LAYOUT = @layout [a{0.68w} b]

# Obs-vs-model timeseries, one panel. `legend`/`title` exposed so grid
# callers (_depth_panels!) can drop the per-panel legend (redundant across a
# grid) and use a short per-panel label instead of a full one.
function _timeseries_panel(t_model, model_vec, obs_vec; title, units, plot_start, plot_end, ylims=nothing, legend=:topleft)
    ps = isnothing(plot_start) ? t_model[1] : plot_start
    pe = isnothing(plot_end) ? t_model[end] : plot_end
    m = findall(t -> ps <= t <= pe, t_model)
    obs_f = Float64[ismissing(v) ? NaN : Float64(v) for v in obs_vec]
    p = plot(t_model[m], Float64.(model_vec[m]); label="model", color=:black, lw=1.2, title, ylabel=units, legend)
    plot!(p, t_model[m], obs_f[m]; label="obs", color=:red, lw=1, alpha=0.8)
    ylims === nothing || ylims!(p, ylims)
    return p
end

# Obs-vs-model 1:1 scatter with stats in the title, one panel.
# `title_prefix` lets grid callers prepend a per-panel label (depth, height).
function _scatter_panel(model_vec, obs_vec; units, ylims=nothing, title_prefix="")
    pairs = [(Float64(o), Float64(mv)) for (o, mv) in zip(obs_vec, model_vec)
             if !ismissing(o) && isfinite(Float64(o)) && isfinite(Float64(mv))]
    isempty(pairs) && return plot(; title="$(title_prefix)(no valid pairs)", titlefontsize=9, framestyle=:none)
    o = first.(pairs); mv = last.(pairs)
    lo, hi = ylims === nothing ? extrema(vcat(o, mv)) : ylims
    s = compute_stats(o, mv)
    p = scatter(o, mv; label=nothing, xlabel="obs ($units)", ylabel="model ($units)",
        title=title_prefix * @sprintf("r=%.3f  RMSE=%.3g\nbias=%+.3g  n=%d", s.r, s.rmse, s.bias, s.n),
        titlefontsize=9, ms=2, alpha=0.35, xlims=(lo, hi), ylims=(lo, hi))
    plot!(p, [lo, hi], [lo, hi]; label="1:1", color=:red, lw=1)
    return p
end

function _plot_pair(t_model, model_vec, obs_vec; label, units, plot_start, plot_end,
                     save_dir=nothing, tag="", display_plots=false, ylims=nothing)
    p_ts = _timeseries_panel(t_model, model_vec, obs_vec; title=label, units, plot_start, plot_end, ylims)
    p_sc = _scatter_panel(model_vec, obs_vec; units, ylims)
    p = plot(p_ts, p_sc; layout=_PAIR_LAYOUT, size=(1000, 350))
    display_plots && display(p)
    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p, joinpath(save_dir, "$(tag).$figure_format"))
    end
    return nothing
end

function report_site_results(prep, output; plot_start=nothing, plot_end=nothing,
                              make_plots=save_outputs, display_plots=false, source="tower", plot_forcing=true,
                              plot_variable=nothing)
    (; site_name, years, hourly, t_model, heights, depths, height_series, resolved, canopy_mode) = prep
    canopy_height = resolved.canopy_height
    tag = "$(site_name)_$(join(years, '-'))_$(source)"
    site_dir(sub) = joinpath(outputs_dir, sub)

    println("\n" * "="^72)
    println("Site $site_name  years=$years  forcing=$source — model vs observations")
    println("="^72)

    site_rows = NamedTuple[]
    record!(label, kind, s) = push!(site_rows,
        (site=site_name, years=join(years, '-'), source, variable=label, kind, r=s.r, rmse=s.rmse, bias=s.bias, n=s.n))
    # `plot_variable`: nothing plots every variable (default); a String
    # restricts plotting to the one variable matching that short key (Fsu,
    # Flu, Fn, Fh, Fe, Fg, Ts, Sws, or a height-series name like
    # "Ta_HMP_2m") -- stats are still computed/printed/recorded for all of
    # them regardless, only plot generation is restricted.
    report1!(label, kind, model_vec, obs_vec, units, save_dir, ptag, key; ylims=nothing) = begin
        s = compute_stats(obs_vec, model_vec)
        println(rpad(label, 32) * fmt_stat(s))
        record!(label, kind, s)
        want_plot = make_plots && (plot_variable === nothing || plot_variable == key)
        want_plot && _plot_pair(t_model, model_vec, obs_vec; label, units, plot_start, plot_end,
            save_dir=(save_outputs ? save_dir : nothing), tag=ptag, display_plots, ylims)
    end

    # Depth loop (Ts/Sws) as two panel grids -- one timeseries grid, one
    # scatter grid -- instead of one combined figure per depth, mirroring
    # plot_canopy_profiles's one-grid-per-variable style. Stats are still
    # computed/printed/recorded for every depth regardless of make_plots/
    # plot_variable, matching report1!'s behavior; only plot generation is
    # gated.
    function _depth_panels!(prefix, kind, model_mat, units, convert, ylims)
        depth_rows = _hourly_depth_series(hourly, prefix)
        isempty(depth_rows) && return nothing
        ts_panels = []
        sc_panels = []
        for r in depth_rows
            depth_m = r.depth_m
            obs_vec = _col(hourly, t_model, "$(prefix)_$(depth_m)m")
            i = argmin(abs.(ustrip.(u"m", depths) .- depth_m))
            model_vec = convert(model_mat[:, i])
            label = "$kind $(depth_m) m"
            s = compute_stats(obs_vec, model_vec)
            println(rpad(label, 32) * fmt_stat(s))
            record!(label, kind, s)
            make_plots && (plot_variable === nothing || plot_variable == prefix) || continue
            push!(ts_panels, _timeseries_panel(t_model, model_vec, obs_vec; title="$(depth_m) m", units,
                plot_start, plot_end, ylims, legend=false))
            push!(sc_panels, _scatter_panel(model_vec, obs_vec; units, ylims, title_prefix="$(depth_m) m\n"))
        end
        isempty(ts_panels) && return nothing
        ncols = min(4, length(ts_panels))
        nrows = cld(length(ts_panels), ncols)
        p_ts = plot(ts_panels...; layout=(nrows, ncols), size=(320 * ncols, 240 * nrows), plot_title="$kind timeseries")
        p_sc = plot(sc_panels...; layout=(nrows, ncols), size=(320 * ncols, 260 * nrows), plot_title="$kind scatter")
        display_plots && (display(p_ts); display(p_sc))
        if save_outputs
            mkpath(site_dir(kind))
            savefig(p_ts, joinpath(site_dir(kind), "$(tag)_timeseries_panels.$figure_format"))
            savefig(p_sc, joinpath(site_dir(kind), "$(tag)_scatter_panels.$figure_format"))
        end
        return nothing
    end

    # ── Forcing sanity: what actually drove the model, real-obs hours in black,
    # gap-filled hours (donor- or interpolation-filled) in orange -- confirms
    # both the read/aggregation pipeline and how much of the record is filled.
    # forcing_used/filled_mask are absent from a SILO-forced prep (tower
    # forcing isn't used at all there); falls back to raw obs-only in that case.
    forcing_used = get(prep, :forcing_used, nothing)
    filled_mask = get(prep, :filled_mask, nothing)
    for (name, units) in (("Ta","°C"), ("RH","%"), ("Ws","m/s"), ("Fsd","W/m^2"), ("Fld","W/m^2"), ("Precip","mm"))
        make_plots && plot_forcing && (plot_variable === nothing || plot_variable == name) || continue
        p = nothing
        if forcing_used !== nothing && haskey(forcing_used, Symbol(name))
            used = getproperty(forcing_used, Symbol(name))
            filled = getproperty(filled_mask, Symbol(name))
            obs_only = [f ? NaN : v for (v, f) in zip(used, filled)]
            filled_only = [f ? v : NaN for (v, f) in zip(used, filled)]
            p = plot(t_model, obs_only; label="obs", color=:black, lw=1, title="$name ($site_name)", ylabel=units)
            plot!(p, t_model, filled_only; label="filled", color=:orange, lw=1)
        else
            vals = _col(hourly, t_model, name)
            vals === nothing && continue
            p = plot(t_model, Float64.(coalesce.(vals, NaN)); label=nothing, title="$name ($site_name)", ylabel=units, color=:black, lw=1)
        end
        display_plots && display(p)
        if save_outputs
            mkpath(site_dir("forcing"))
            savefig(p, joinpath(site_dir("forcing"), "$(tag)_$(name).$figure_format"))
        end
    end

    # ── Fsu, Flu, Fn: canopy-top boundary (column 1 = canopy top -> ground).
    # NoCanopy (canopy_mode=:legacy) has zero boundary columns -- no two-stream
    # radiative transfer to compare, skip. ───────────────────────────────────
    if canopy_mode != :legacy
        fsu = ustrip.(u"W/m^2", output.canopy.boundary_upward_shortwave[:, 1])
        flu = ustrip.(u"W/m^2", output.canopy.boundary_upward_longwave[:, 1])
        fn = ustrip.(u"W/m^2",
            output.canopy.boundary_downward_shortwave[:, 1] .- output.canopy.boundary_upward_shortwave[:, 1] .+
            output.canopy.boundary_downward_longwave[:, 1] .- output.canopy.boundary_upward_longwave[:, 1])
        for (label, model_vec, obscol) in (("Upward shortwave (Fsu)", fsu, "Fsu"), ("Upward longwave (Flu)", flu, "Flu"), ("Net radiation (Fn)", fn, "Fn"))
            obs_vec = _col(hourly, t_model, obscol)
            obs_vec === nothing && continue
            report1!(label, "radiation", model_vec, obs_vec, "W/m^2", site_dir(_FIG_SUBDIR[obscol]), tag, obscol)
        end
    end

    make_plots && _plot_leaf_vs_air(heights, canopy_height, canopy_mode, t_model, output; plot_start, plot_end,
        save_dir=(save_outputs ? site_dir("leaf_temperature") : nothing), tag, display_plots)

    # ── Fh, Fe, Fg. Sign convention (leaf sensible/latent surface->atmosphere
    # positive, Fg downward positive) confirmed against the tower convention
    # directly -- no negation needed. canopy_sensible_heat_flux only sums the
    # leaf layers -- canopy_energy_balance! uses ground_temperature purely for
    # longwave exchange with the canopy, never computing the bare ground's own
    # convective loss to the atmosphere. That's computed separately, via the
    # same free-atmosphere MOST routine :legacy relies on entirely
    # (profile.convective_heat_flux, atmosphere->surface positive, hence the
    # sign flip) -- added back in here too, since a tower's Fh footprint
    # integrates canopy AND exposed ground together, and for a sparse canopy
    # (low PAI) the ground term dominates, not the leaves. There's no
    # equivalent ground-level Fe output to add for latent heat (bare-soil
    # evaporation isn't exposed as a separate flux here). ────────────────────
    fh = canopy_mode == :legacy ?
        -ustrip.(u"W/m^2", output.profile.convective_heat_flux) :
        ustrip.(u"W/m^2", output.canopy.canopy_sensible_heat_flux) .- ustrip.(u"W/m^2", output.profile.convective_heat_flux)
    fe = ustrip.(u"W/m^2", output.canopy.canopy_latent_heat_flux)
    fg = ustrip.(u"W/m^2", output.ground_heat_flux)
    flux_vars = canopy_mode == :legacy ?
        (("Sensible heat (Fh)", fh, "Fh"), ("Ground heat flux (Fg)", fg, "Fg")) :
        (("Sensible heat (Fh)", fh, "Fh"), ("Latent heat (Fe)", fe, "Fe"), ("Ground heat flux (Fg)", fg, "Fg"))
    for (label, model_vec, obscol) in flux_vars
        obs_vec = _col(hourly, t_model, obscol)
        obs_vec === nothing && continue
        report1!(label, "flux", model_vec, obs_vec, "W/m^2", site_dir(_FIG_SUBDIR[obscol]), tag, obscol)
    end

    # ── Soil temperature / moisture at every observed depth. Soil moisture
    # y-axis fixed (0-0.6 frac) across every depth/site so plots are directly
    # comparable, rather than each auto-scaling to its own range. Soil
    # temperature range varies too much by site/climate for a single
    # constant, so it's fixed to the shallowest depth's own range instead
    # (widest diurnal swing) -- deeper depths then visibly damp toward the
    # middle of that range rather than each auto-scaling to its own
    # (narrower) extent. ──────────────────────────────────────────────────
    ts_series = _hourly_depth_series(hourly, "Ts")
    ts_ylims = if isempty(ts_series)
        nothing
    else
        shallow_depth_m = minimum(r.depth_m for r in ts_series)
        i = argmin(abs.(ustrip.(u"m", depths) .- shallow_depth_m))
        vals = collect(skipmissing(ustrip.(u"°C", output.soil_temperature[:, i])))
        obs_vec = _col(hourly, t_model, "Ts_$(shallow_depth_m)m")
        obs_vec === nothing || append!(vals, collect(skipmissing(obs_vec)))
        if isempty(vals)
            nothing
        else
            lo, hi = extrema(vals)
            pad = 0.05 * max(hi - lo, 1.0e-6)
            (lo - pad, hi + pad)
        end
    end
    for (prefix, kind, model_mat, units, convert, ylims) in (
        ("Ts", "soil_temperature", output.soil_temperature, "°C", d -> ustrip.(u"°C", d), ts_ylims),
        ("Sws", "soil_moisture", output.soil_moisture, "frac", identity, (0.0, 0.6)),
    )
        _depth_panels!(prefix, kind, model_mat, units, convert, ylims)
    end

    # ── Multi-height Ta/Ws profiles -- routed to canopy- or profile-output
    # per height, to check the vertical shape (not just tower-height flux). ──
    for (height_m, name) in height_series
        field = _height_field(name)
        field === nothing && continue
        obs_vec = _col(hourly, t_model, name)
        obs_vec === nothing && continue
        raw = _profile_series(output, heights, canopy_height, height_m, field, canopy_mode)
        raw === nothing && continue
        model_vec, units = field == :air_temperature ? (ustrip.(u"°C", u"°C".(raw)), "°C") : (ustrip.(u"m/s", raw), "m/s")
        report1!("$name @ $(height_m) m", "profile_$(field)", model_vec, obs_vec, units, site_dir("profiles"), "$(tag)_$(name)", name)
    end

    stats_df = DataFrame(site_rows)
    if save_outputs
        mkpath(site_dir("stats"))
        CSV.write(joinpath(site_dir("stats"), "$(tag).csv"), stats_df)
    end
    return stats_df
end
