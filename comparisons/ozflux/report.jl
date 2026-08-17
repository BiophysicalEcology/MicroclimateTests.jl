# report.jl — stats + plots for every comparison target, given a solved
# `run_site` result. Requires config.jl/utils.jl/pipeline.jl already included.

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

# Rescales an obs column by a constant factor, `missing` untouched -- used to
# test whether a site's Sws sensor is reporting degree-of-saturation
# (fraction of pore space) rather than the volumetric water content its own
# metadata claims (see Longreach: obs peaks ~0.77, above plausible porosity).
_scale_obs(obs_vec::Nothing, scale) = nothing
_scale_obs(obs_vec, scale) = scale == 1.0 ? obs_vec : [ismissing(v) ? v : v * scale for v in obs_vec]

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
    (; heights, resolved, hourly, height_series, canopy_mode, t_model, depths) = prep
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
        marker=:circle, ms=3, markerstrokewidth=0,
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

    # Predicted soil surface temperature (shallowest depth node) as a single
    # point at height 0 -- context for the air profile's ground boundary,
    # without the clutter of the full depth profile (see plot_soil_profiles).
    if field == :air_temperature
        surface_i = argmin(ustrip.(u"m", depths))
        surface_val = convert(output.soil_temperature[i, surface_i])
        scatter!(p, [surface_val], [0.0]; label=nothing, color=:black, ms=4, markerstrokewidth=0)
    end
    return p
end

# Common x-axis range across every panel in one plot_canopy_profiles grid
# (model line + any obs scatter) -- panels default to independent auto-ranges
# otherwise, making the shape across hours hard to compare by eye.
function _shared_profile_xlims(prep, output, profile_times; field, convert)
    (; heights, resolved, hourly, height_series, canopy_mode, t_model, depths) = prep
    canopy_height = resolved.canopy_height
    vals = Float64[]
    for dt in profile_times
        i = findfirst(==(dt), t_model)
        i === nothing && continue
        append!(vals, Iterators.filter(isfinite, skipmissing(_vertical_profile(output, heights, canopy_height, i, field, canopy_mode, convert))))
        if field == :air_temperature
            v = convert(output.soil_temperature[i, argmin(ustrip.(u"m", depths))])
            isfinite(v) && push!(vals, v)
        end
        row = findfirst(==(dt), hourly.DateTime)
        row === nothing && continue
        for (h, name) in height_series
            (_height_field(name) == field && name in names(hourly)) || continue
            v = hourly[row, name]
            (!ismissing(v) && isfinite(Float64(v))) && push!(vals, Float64(v))
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

# Soil temperature vs depth at a single hour, plus a near-ground air
# temperature continuation (0 to 0.5 m, model + obs) -- same conventions as
# _profile_panel (model line+points, obs points), depth plotted as negative
# height so 0 is the surface. Kept as separate line segments either side of
# 0 (not one connected polyline) so any real soil/air discontinuity stays
# visible rather than being smoothed over. y-axis fixed to (-2, 0.5) m
# regardless of the model's own deepest node, so panels stay comparable
# across sites even when extra depth nodes are added for a deep observed
# sensor.
function _soil_profile_panel(prep, output, dt; xlims=:auto)
    (; hourly, t_model, depths, heights, resolved, canopy_mode, height_series) = prep
    canopy_height = resolved.canopy_height
    title = Dates.format(dt, "yyyy-mm-dd HH:MM")
    i = findfirst(==(dt), t_model)
    if i === nothing
        range_str = "$(first(t_model)) to $(last(t_model))"
        return plot(; title="$title\noutside solved range\n($range_str)", titlefontsize=7, framestyle=:none)
    end
    depths_m = -ustrip.(u"m", depths)
    soil_vals = ustrip.(u"°C", output.soil_temperature[i, :])
    order = sortperm(depths_m)
    p = plot(soil_vals[order], depths_m[order]; label=nothing, color=:black, lw=1.5,
        marker=:circle, ms=3, markerstrokewidth=0,
        xlabel="°C", ylabel="height (m)", title, titlefontsize=8, legend=false, xlims, ylims=(-2.0, 0.5))

    row = findfirst(==(dt), hourly.DateTime)
    if row !== nothing
        for r in _hourly_depth_series(hourly, "Ts")
            v = r.values[row]
            ismissing(v) || scatter!(p, [Float64(v)], [-r.depth_m]; label=nothing, color=:red, ms=6,
                markerstrokewidth=1.5, markerstrokecolor=:white)
        end
    end

    _plot_near_ground_air!(p, output, heights, canopy_height, i, canopy_mode, :black)
    row !== nothing && _scatter_near_ground_obs!(p, hourly, row, height_series)
    return p
end

# Model air temperature at heights <= NEAR_GROUND_TOP_M, appended to an
# existing soil-profile panel `p`. Shared by the single- and 3-way variants.
const NEAR_GROUND_TOP_M = 0.5

# Indices of z_m at or below NEAR_GROUND_TOP_M, falling back to just the
# single closest-to-ground index when none qualify -- a coarse independent
# grid (e.g. micropoint's own evenly-spaced canopy_height/50, whose lowest
# point can sit just above NEAR_GROUND_TOP_M for a tall canopy) would
# otherwise show nothing at all near the ground.
function _near_ground_indices(z_m, top_m)
    keep = findall(<=(top_m), z_m)
    return isempty(keep) ? [argmin(z_m)] : keep
end
function _plot_near_ground_air!(p, output, heights, canopy_height, i, canopy_mode, color)
    air_vals = _vertical_profile(output, heights, canopy_height, i, :air_temperature, canopy_mode, t -> ustrip(u"°C", t))
    heights_m = ustrip.(u"m", heights)
    valid = [k for k in eachindex(heights_m) if !ismissing(air_vals[k])]
    isempty(valid) && return nothing
    near = filter(k -> heights_m[k] <= NEAR_GROUND_TOP_M, valid)
    keep = isempty(near) ? [valid[argmin(heights_m[valid])]] : near
    ka = keep[sortperm(heights_m[keep])]
    plot!(p, air_vals[ka], heights_m[ka]; label=nothing, color, lw=1.5,
        marker=:circle, ms=3, markerstrokewidth=0)
    return nothing
end

# Observed air temperature at height_series sensors <= NEAR_GROUND_TOP_M.
function _scatter_near_ground_obs!(p, hourly, row, height_series)
    for (h, name) in height_series
        (h <= NEAR_GROUND_TOP_M && _height_field(name) == :air_temperature && name in names(hourly)) || continue
        v = hourly[row, name]
        ismissing(v) || scatter!(p, [Float64(v)], [h]; label=nothing, color=:red, ms=6,
            markerstrokewidth=1.5, markerstrokecolor=:white)
    end
    return nothing
end

# Shared x-axis range for a plot_soil_profiles grid, mirroring
# _shared_profile_xlims.
function _shared_soil_profile_xlims(prep, output, profile_times)
    (; hourly, t_model, heights, resolved, canopy_mode, height_series) = prep
    canopy_height = resolved.canopy_height
    vals = Float64[]
    for dt in profile_times
        i = findfirst(==(dt), t_model)
        i === nothing && continue
        append!(vals, Iterators.filter(isfinite, ustrip.(u"°C", output.soil_temperature[i, :])))
        air_vals = _vertical_profile(output, heights, canopy_height, i, :air_temperature, canopy_mode, t -> ustrip(u"°C", t))
        heights_m = ustrip.(u"m", heights)
        for k in eachindex(heights_m)
            (heights_m[k] <= NEAR_GROUND_TOP_M && !ismissing(air_vals[k]) && isfinite(air_vals[k])) && push!(vals, air_vals[k])
        end
        row = findfirst(==(dt), hourly.DateTime)
        row === nothing && continue
        for r in _hourly_depth_series(hourly, "Ts")
            v = r.values[row]
            (!ismissing(v) && isfinite(Float64(v))) && push!(vals, Float64(v))
        end
        for (h, name) in height_series
            (h <= NEAR_GROUND_TOP_M && _height_field(name) == :air_temperature && name in names(hourly)) || continue
            v = hourly[row, name]
            (!ismissing(v) && isfinite(Float64(v))) && push!(vals, Float64(v))
        end
    end
    isempty(vals) && return :auto
    lo, hi = extrema(vals)
    pad = 0.05 * max(hi - lo, 1.0e-6)
    return (lo - pad, hi + pad)
end

# Grid of _soil_profile_panel, one panel per profile_times entry.
function plot_soil_profiles(prep, output, profile_times; save_dir=nothing, tag="", display_plots=false)
    xlims = _shared_soil_profile_xlims(prep, output, profile_times)
    panels = [_soil_profile_panel(prep, output, dt; xlims) for dt in profile_times]
    ncols = min(4, length(panels))
    nrows = cld(length(panels), ncols)
    p = plot(panels...; layout=(nrows, ncols), size=(320 * ncols, 260 * nrows),
        plot_title="soil_temperature profile — $(prep.site_name)")
    display_plots && display(p)
    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p, joinpath(save_dir, "$(tag)_soil_temperature_profiles.$figure_format"))
    end
    return p
end

# Per-layer plant area index vs height -- the canopy structure the solve
# actually used, not time-varying (unlike the other profile plots), so a
# single panel rather than a profile_times grid. :legacy (NoCanopy) has no
# per-layer PAI -- skipped.
function plot_plant_area_index(prep; save_dir=nothing, tag="", display_plots=false)
    canopy_model = prep.problem.model.canopy_model
    canopy_model isa NoCanopy && return nothing
    pai = canopy_model.plant_area_index
    layer_heights = Microclimate.canopy_layer_heights(prep.heights, canopy_model.canopy_height, length(pai)).layer_heights
    heights_m = ustrip.(u"m", layer_heights)
    order = sortperm(heights_m)
    p = plot(pai[order], heights_m[order]; label=nothing, color=:black, lw=1.5,
        marker=:circle, ms=4, markerstrokewidth=0,
        xlabel="plant area index (m²/m²)", ylabel="height (m)", title="Plant area index — $(prep.site_name)")
    display_plots && display(p)
    if save_dir !== nothing
        mkpath(save_dir)
        savefig(p, joinpath(save_dir, "$(tag)_plant_area_index.$figure_format"))
    end
    return p
end

# Ta + Ws + RH profile grids for every hour in `profile_times` -- call this
# directly from single_site.jl. Most meaningful at sites with real sub-canopy
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
    plot_soil_profiles(prep, output, profile_times; save_dir, tag, display_plots)
    plot_plant_area_index(prep; save_dir, tag, display_plots)
    plot_canopy_profiles(prep, output, profile_times; field=:wind_speed, units="m/s",
        convert=w -> ustrip(u"m/s", w), save_dir, tag, display_plots)
    output_ah = _with_absolute_humidity(output)
    plot_canopy_profiles(prep, output_ah, profile_times; field=:absolute_humidity, units="g/m^3",
        convert=identity, save_dir, tag, display_plots)
    return nothing
end

# Top-of-canopy leaf temperature vs the local (same-layer) air temperature --
# direct diagnostic for how far leaves run from air, independent of any tower
# obs. Useful when tuning canopy_convergence_model_choice:
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

# Obs-vs-model timeseries, one panel. `legend`/`title` exposed so grid
# callers (_panel_grids!) can drop the per-panel legend (redundant across a
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

    # A group of same-`kind` variables as two panel grids -- one timeseries
    # grid, one scatter grid -- instead of one combined figure per variable
    # (mirrors plot_canopy_profiles's one-grid-per-variable style). Size
    # floors keep a single-panel grid from becoming a tiny image. Stats are
    # always computed/printed/recorded; plot_variable only gates which
    # panels enter the grid.
    function _panel_grids!(kind, vars; ylims=nothing)
        ts_panels = []
        sc_panels = []
        for (label, model_vec, obs_vec, units, key, panel_title) in vars
            obs_vec === nothing && continue
            s = compute_stats(obs_vec, model_vec)
            println(rpad(label, 32) * fmt_stat(s))
            record!(label, kind, s)
            make_plots && (plot_variable === nothing || plot_variable == key) || continue
            push!(ts_panels, _timeseries_panel(t_model, model_vec, obs_vec; title=panel_title, units,
                plot_start, plot_end, ylims, legend=false))
            push!(sc_panels, _scatter_panel(model_vec, obs_vec; units, ylims, title_prefix="$panel_title\n"))
        end
        isempty(ts_panels) && return nothing
        ncols = min(4, length(ts_panels))
        nrows = cld(length(ts_panels), ncols)
        # Minimum sizes prevent a single-panel figure becoming only 320 × 240 px.
        ts_size = (max(900, 420 * ncols), max(500, 300 * nrows),)
        scatter_size = (max(700, 360 * ncols), max(600, 360 * nrows),)
        common_kwargs = (; layout=(nrows, ncols), dpi=150, margin=4mm, left_margin=7mm,
                            right_margin=4mm, top_margin=8mm, bottom_margin=10mm,
                            plot_titlefontsize=14,tickfontsize=8, guidefontsize=10, titlefontsize=10,
                         )
        p_ts = plot(ts_panels...; common_kwargs..., size=ts_size, plot_title="$kind timeseries", xrotation=30,)
        p_sc = plot(sc_panels...; common_kwargs..., size=scatter_size, plot_title="$kind scatter")
        display_plots && (display(p_ts); display(p_sc))
        if save_outputs
            mkpath(site_dir(kind))
            savefig(p_ts, joinpath(site_dir(kind), "$(tag)_timeseries_panels.$figure_format"))
            savefig(p_sc, joinpath(site_dir(kind), "$(tag)_scatter_panels.$figure_format"))
        end
        return nothing
    end
    _depth_panels!(prefix, kind, model_mat, units, convert, ylims; obs_scale=1.0) = _panel_grids!(kind,
        [("$kind $(r.depth_m) m", convert(model_mat[:, argmin(abs.(ustrip.(u"m", depths) .- r.depth_m))]),
          _scale_obs(_col(hourly, t_model, "$(prefix)_$(r.depth_m)m"), obs_scale), units, prefix, "$(r.depth_m) m")
         for r in _hourly_depth_series(hourly, prefix)];
        ylims)

    # ── Forcing sanity: what actually drove the model, real-obs hours in black,
    # gap-filled hours (donor- or interpolation-filled) in orange -- confirms
    # both the read/aggregation pipeline and how much of the record is filled.
    # forcing_used/filled_mask are absent from a SILO-forced prep (tower
    # forcing isn't used there at all) -- falls back to what SILO's own model
    # actually consumed (output.reference_*/global_radiation, Fld backed out
    # of sky_temperature same as silo_gapfill_donor, environment_daily's
    # daily rainfall total) instead of the tower's raw columns, which a SILO
    # run never reads.
    forcing_used = get(prep, :forcing_used, nothing)
    filled_mask = get(prep, :filled_mask, nothing)
    environment_daily = get(prep, :environment_daily, nothing)
    sim_start = get(prep, :sim_start, nothing)
    for (name, units) in (("Ta","°C"), ("RH","%"), ("Ws","m/s"), ("Fsd","W/m^2"), ("Fld","W/m^2"), ("Precip","mm"))
        make_plots && plot_forcing && (plot_variable === nothing || plot_variable == name) || continue
        p = nothing
        ps = isnothing(plot_start) ? t_model[1] : plot_start
        pe = isnothing(plot_end) ? t_model[end] : plot_end
        m = findall(t -> ps <= t <= pe, t_model)
        if forcing_used !== nothing && haskey(forcing_used, Symbol(name))
            used = getproperty(forcing_used, Symbol(name))
            filled = getproperty(filled_mask, Symbol(name))
            obs_only = [f ? NaN : v for (v, f) in zip(used, filled)]
            filled_only = [f ? v : NaN for (v, f) in zip(used, filled)]
            p = plot(t_model[m], obs_only[m]; label="obs", color=:black, lw=1, title="$name ($site_name)", ylabel=units)
            plot!(p, t_model[m], filled_only[m]; label="filled", color=:orange, lw=1)
        elseif environment_daily !== nothing && name == "Precip"
            # SILO has no hourly rain, only a daily total -- plotted at daily
            # resolution (a step per day), not resampled onto t_model.
            daily_dates = collect(sim_start:Day(1):(sim_start + Day(length(environment_daily.rainfall) - 1)))
            vals = ustrip.(u"kg/m^2", environment_daily.rainfall)  # mm == kg/m^2 for water
            dm = findall(d -> Date(ps) <= d <= Date(pe), daily_dates)
            isempty(dm) && continue
            p = plot(daily_dates[dm], vals[dm]; label=nothing, title="$name ($site_name, SILO daily total)",
                ylabel=units, color=:black, lw=1, seriestype=:steppost)
        elseif environment_daily !== nothing
            vals = name == "Ta" ? ustrip.(u"°C", output.reference_temperature) :
                   name == "RH" ? output.reference_humidity .* 100.0 :
                   name == "Ws" ? ustrip.(u"m/s", output.reference_wind_speed) :
                   name == "Fsd" ? ustrip.(u"W/m^2", output.global_radiation) :
                   ustrip.(u"W/m^2", FluidProperties.σ .* output.sky_temperature .^ 4)  # Fld
            p = plot(t_model[m], vals[m]; label=nothing, title="$name ($site_name, SILO-derived)", ylabel=units, color=:black, lw=1)
        else
            vals = _col(hourly, t_model, name)
            vals === nothing && continue
            p = plot(t_model[m], Float64.(coalesce.(vals[m], NaN)); label=nothing, title="$name ($site_name)", ylabel=units, color=:black, lw=1)
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
        _panel_grids!("radiation", [
            (label, model_vec, _col(hourly, t_model, obscol), "W/m^2", obscol, obscol)
            for (label, model_vec, obscol) in
                (("Upward shortwave (Fsu)", fsu, "Fsu"), ("Upward longwave (Flu)", flu, "Flu"), ("Net radiation (Fn)", fn, "Fn"))
        ])
    end

    make_plots && _plot_leaf_vs_air(heights, canopy_height, canopy_mode, t_model, output; plot_start, plot_end,
        save_dir=(save_outputs ? site_dir("leaf_temperature") : nothing), tag, display_plots)

    # ── Fh, Fe, Fg. Sign convention (leaf sensible/latent surface->atmosphere
    # positive, Fg downward positive) confirmed against the tower convention
    # directly -- no negation needed. :legacy has one bare surface, so
    # profile.convective_heat_flux (atmosphere->surface positive, hence the
    # sign flip) is its Fh. MultilayerCanopy doesn't model an
    # exposed/gap ground fraction -- the whole ground surface exchanges with
    # the atmosphere *through* the canopy air column (ground_heat_conductance,
    # already reflected in canopy_sensible_heat_flux and the soil solve), not
    # via a separate direct pathway. profile.convective_heat_flux for this
    # canopy_mode is profile_surface_temperature's top-leaf-temperature MOST
    # profile (shapes the above-canopy air/wind profile continuously down to
    # canopy top) -- a second, independent, leaf-driven flux estimate, not a
    # ground term; subtracting it from canopy_sensible_heat_flux cancels two
    # correlated leaf-driven quantities instead of adding a real ground
    # contribution, so it's dropped here. ──────────────────────────────────
    fh = canopy_mode == :legacy ?
        -ustrip.(u"W/m^2", output.profile.convective_heat_flux) :
        ustrip.(u"W/m^2", output.canopy.canopy_sensible_heat_flux)
    fe = ustrip.(u"W/m^2", output.canopy.canopy_latent_heat_flux)
    fg = ustrip.(u"W/m^2", output.ground_heat_flux)
    flux_vars = canopy_mode == :legacy ?
        (("Sensible heat (Fh)", fh, "Fh"), ("Ground heat flux (Fg)", fg, "Fg")) :
        (("Sensible heat (Fh)", fh, "Fh"), ("Latent heat (Fe)", fe, "Fe"), ("Ground heat flux (Fg)", fg, "Fg"))
    _panel_grids!("flux", [
        (label, model_vec, _col(hourly, t_model, obscol), "W/m^2", obscol, obscol)
        for (label, model_vec, obscol) in flux_vars
    ])

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
    # sws_obs_scale (utils.jl): Longreach's Sws sensor appears to report
    # degree-of-saturation rather than true volumetric water content --
    # rescaled by the site's own void fraction, a no-op elsewhere. Also
    # applied to the initial-condition read (pipeline.jl), same raw column.
    for (prefix, kind, model_mat, units, convert, ylims, obs_scale) in (
        ("Ts", "soil_temperature", output.soil_temperature, "°C", d -> ustrip.(u"°C", d), ts_ylims, 1.0),
        ("Sws", "soil_moisture", output.soil_moisture, "frac", identity, (0.0, 0.6), sws_obs_scale(site_name, prep.problem.inputs.soil_profile)),
    )
        _depth_panels!(prefix, kind, model_mat, units, convert, ylims; obs_scale)
    end

    # ── Multi-height Ta/Ws profiles -- routed to canopy- or profile-output
    # per height, to check the vertical shape (not just tower-height flux).
    # Ta and Ws grouped into their own panel grid (not one figure per
    # height). ─────────────────────────────────────────────────────────────
    ta_vars = []
    ws_vars = []
    for (height_m, name) in height_series
        field = _height_field(name)
        field === nothing && continue
        obs_vec = _col(hourly, t_model, name)
        obs_vec === nothing && continue
        raw = _profile_series(output, heights, canopy_height, height_m, field, canopy_mode)
        raw === nothing && continue
        model_vec, units = field == :air_temperature ? (ustrip.(u"°C", u"°C".(raw)), "°C") : (ustrip.(u"m/s", raw), "m/s")
        label = "$name @ $(height_m) m"
        push!(field == :air_temperature ? ta_vars : ws_vars, (label, model_vec, obs_vec, units, name, label))
    end
    _panel_grids!("profile_air_temperature", ta_vars)
    _panel_grids!("profile_wind_speed", ws_vars)

    stats_df = DataFrame(site_rows)
    if save_outputs
        mkpath(site_dir("stats"))
        CSV.write(joinpath(site_dir("stats"), "$(tag).csv"), stats_df)
    end
    return stats_df
end
