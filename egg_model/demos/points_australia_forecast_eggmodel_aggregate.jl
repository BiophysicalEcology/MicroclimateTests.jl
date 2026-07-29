# Stage 8, HPC split egg-model stage (aggregate): loads the historical
# egg-model cache and every member's cached egg-model outcomes (both written
# by points_australia_forecast_eggmodel_historical.jl /
# points_australia_forecast_eggmodel_members.jl) and produces the combined
# stats CSV + maps. No egg-model solving here, so no egg_model/pars/etc.
# construction either -- ThermalPhysiology & friends are still `using`'d
# because the cached results deserialize into their types (EggState etc.).
#
# Safe to run standalone once every member's cache file exists.

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", "c:/Spatial_Data/")

using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Statistics
using NaturalEarth, GeoInterface
using Serialization

include(joinpath(@__DIR__, "points_australia_forecast_setup.jl"))

historical_egg_cache_path = joinpath(egg_dir, "$(historical_label())_$(diapause ? "dia" : "nodia")_egg_n$(n).jls")
isfile(historical_egg_cache_path) || error(
    "Historical egg-model cache not found at $historical_egg_cache_path -- " *
    "run points_australia_forecast_eggmodel_historical.jl first.")
println("Loading cached historical egg-model results...")
historical_egg_results = deserialize(historical_egg_cache_path)

n_hatched_historically = count(r -> r.hatched, historical_egg_results)
n_died_historically = count(r -> r.died, historical_egg_results)
println("$n_hatched_historically/$n hatched and $n_died_historically/$n died before $issue_date; " *
        "$(n - n_hatched_historically - n_died_historically)/$n continue into the forecast.")
forecast_point_indices = [i for i in 1:n
                           if !historical_egg_results[i].hatched && !historical_egg_results[i].died]

# historical_egg_results is concretely-typed (Vector{ResultType}, written by
# points_australia_forecast_eggmodel_historical.jl), so its eltype gives the
# same concrete result type each member's outcomes file was written with --
# recovered here without needing the egg_model/pars construction at all.
ResultType = eltype(historical_egg_results)
forecast_outcomes = Matrix{Union{Nothing,ResultType}}(nothing, n, length(members))

if !isempty(forecast_point_indices)
    for (member_index, member) in pairs(members)
        outcome_path = joinpath(egg_dir, "$(forecast_label(member))_$(diapause ? "dia" : "nodia")_eggout_n$(n).jls")
        isfile(outcome_path) || error(
            "Missing egg-model outcomes for member $member at $outcome_path -- " *
            "did every array task in points_australia_forecast_eggmodel_members.jl finish?")
        member_outcomes = deserialize(outcome_path)
        forecast_outcomes[forecast_point_indices, member_index] .= member_outcomes[forecast_point_indices]
    end
end

# ── per-point summary: median/p25/p75 hatching member, death-cause
# counts, saved to one combined CSV ──

function summarize_point(historical_egg_result, point_outcomes, members)
    if historical_egg_result.hatched
        return (; status=:hatched_historically, counts=Dict(:hatched => 1),
            p25=nothing, median=nothing, p75=nothing,
            cv_hatch_days=nothing, cv_hatch_mass_mg=nothing)
    elseif historical_egg_result.died
        return (; status=Symbol("died_historically_$(historical_egg_result.death_cause)"),
            counts=Dict(historical_egg_result.death_cause => 1),
            p25=nothing, median=nothing, p75=nothing,
            cv_hatch_days=nothing, cv_hatch_mass_mg=nothing)
    end
    counts = Dict{Symbol,Int}()
    for r in point_outcomes
        code = r.hatched ? :hatched : r.died ? r.death_cause : :timeout
        counts[code] = get(counts, code, 0) + 1
    end
    hatched = [(member, r) for (member, r) in zip(members, point_outcomes) if r.hatched]
    if isempty(hatched)
        return (; status=:no_hatch, counts, p25=nothing, median=nothing, p75=nothing,
            cv_hatch_days=nothing, cv_hatch_mass_mg=nothing)
    end
    sorted = sort(hatched; by=x -> x[2].hatch_time)
    n_hatched = length(sorted)
    # nearest-rank percentile: the member whose hatch_time sits at that rank,
    # so the reported mass stays paired with the same member as its date --
    # not an independently-interpolated order statistic (which could pick a
    # different member's mass than the one at that date rank).
    p25_idx = clamp(ceil(Int, 0.25 * n_hatched), 1, n_hatched)
    p75_idx = clamp(ceil(Int, 0.75 * n_hatched), 1, n_hatched)
    # CV (%) needs >=2 points to be meaningful -- nothing (blank in the CSV)
    # otherwise. Scale-free unlike raw std -- "8% CV" reads the same whether
    # hatching takes 40 days or 90.
    cv_hatch_days = if n_hatched >= 2
        day_vals = [ustrip(u"d", r.hatch_time) for (_, r) in hatched]
        100 * std(day_vals) / mean(day_vals)
    else
        nothing
    end
    cv_hatch_mass_mg = if n_hatched >= 2
        mass_vals = [ustrip(u"mg", r.final_state.egg_mass) for (_, r) in hatched]
        100 * std(mass_vals) / mean(mass_vals)
    else
        nothing
    end
    (; status=:forecast, counts,
        p25=sorted[p25_idx], median=sorted[cld(n_hatched, 2)], p75=sorted[p75_idx],
        cv_hatch_days, cv_hatch_mass_mg)
end

hatch_date_of(member_result, issue_date) =
    issue_date + Day(round(Int, ustrip(u"d", member_result[2].hatch_time)))
mass_mg_of(member_result) = ustrip(u"mg", member_result[2].final_state.egg_mass)

point_summaries = [summarize_point(historical_egg_results[i],
    [forecast_outcomes[i, m] for m in eachindex(members) if forecast_outcomes[i, m] !== nothing], members)
    for i in 1:n]

stats_path = joinpath(egg_dir, "splice_stats_$(run_tag).csv")
open(stats_path, "w") do io
    println(io, "lon,lat,status,n_members,n_hatched,n_died_cold,n_died_heat,n_died_desiccation,n_timeout," *
                 "p25_hatch_date,p25_hatch_mass_mg,median_hatch_date,median_hatch_mass_mg," *
                 "p75_hatch_date,p75_hatch_mass_mg,cv_hatch_days_pct,cv_hatch_mass_pct")
    for i in 1:n
        s = point_summaries[i]
        n_run = s.p25 === nothing && s.status != :no_hatch ? 0 : length(members)
        c(sym) = get(s.counts, sym, 0)
        fields = if s.p25 === nothing
            ("", "", "", "", "", "")
        else
            (string(hatch_date_of(s.p25, issue_date)), string(round(mass_mg_of(s.p25); digits=2)),
             string(hatch_date_of(s.median, issue_date)), string(round(mass_mg_of(s.median); digits=2)),
             string(hatch_date_of(s.p75, issue_date)), string(round(mass_mg_of(s.p75); digits=2)))
        end
        cv_days_str = s.cv_hatch_days === nothing ? "" : string(round(s.cv_hatch_days; digits=2))
        cv_mass_str = s.cv_hatch_mass_mg === nothing ? "" : string(round(s.cv_hatch_mass_mg; digits=2))
        println(io, join((points[i][1], points[i][2], s.status, n_run, c(:hatched), c(:cold), c(:heat),
            c(:desiccation), c(:timeout), fields..., cv_days_str, cv_mass_str), ","))
    end
end
println("Saved $stats_path")

# ── maps: heatmap over the regular grid + basemap (points_australia.jl) ──

using Plots

point_to_index = Dict(p => i for (i, p) in enumerate(all_grid_points))
to_heatmap_z(vals) = permutedims(reshape(vals, length(lon_range), length(lat_range)))

const GI = GeoInterface
function _collect_boundary!(xs, ys, geom)
    trait = GI.geomtrait(geom)
    if trait isa GI.AbstractPointTrait
        push!(xs, GI.x(geom)); push!(ys, GI.y(geom))
    elseif trait isa GI.LinearRingTrait || trait isa GI.LineStringTrait
        for pt in GI.getpoint(geom)
            push!(xs, GI.x(pt)); push!(ys, GI.y(pt))
        end
        push!(xs, NaN); push!(ys, NaN)
    else
        for sub in GI.getgeom(geom)
            _collect_boundary!(xs, ys, sub)
        end
    end
end
au_states = naturalearth("admin_1_states_provinces", 10)
state_xs, state_ys = Float64[], Float64[]
for i in 1:length(au_states)
    au_states[i].iso_a2 == "AU" || continue
    _collect_boundary!(state_xs, state_ys, au_states[i].geometry)
end

map_towns = ["Birdsville, Queensland", "Roma, Queensland", "Charleville, Queensland",
             "Dubbo, New South Wales", "Broken Hill, New South Wales", "Bourke, New South Wales",
             "Mildura, Victoria", "Canberra, Australia"]
towns_cache_path = joinpath(egg_dir, "points_australia_towns.jls")
towns = if isfile(towns_cache_path)
    deserialize(towns_cache_path)
else
    t = map(name -> geocode(name), map_towns)
    serialize(towns_cache_path, t)
    t
end

function add_basemap!(p)
    plot!(p, state_xs, state_ys; color=:black, linewidth=0.75, label=nothing)
    scatter!(p, [t.lon for t in towns], [t.lat for t in towns]; color=:black, markersize=2, label=nothing)
    for t in towns
        annotate!(p, t.lon, t.lat, text(first(split(t.display_name, ",")), 6, :left, :bottom))
    end
    xlims!(p, extrema(lon_range)...)
    ylims!(p, extrema(lat_range)...)
    p
end

historical_hatch_date(r) = oviposition_date + Day(round(Int, ustrip(u"d", r.hatch_time)))

median_dates = fill(NaN, length(all_grid_points))
p25_dates = fill(NaN, length(all_grid_points))
p75_dates = fill(NaN, length(all_grid_points))
median_mass = fill(NaN, length(all_grid_points))
p25_mass = fill(NaN, length(all_grid_points))
p75_mass = fill(NaN, length(all_grid_points))
cv_dates = fill(NaN, length(all_grid_points))
cv_mass = fill(NaN, length(all_grid_points))
frac_hatched = fill(NaN, length(all_grid_points))
frac_cold = fill(NaN, length(all_grid_points))
frac_heat = fill(NaN, length(all_grid_points))
frac_desiccation = fill(NaN, length(all_grid_points))
for i in 1:n
    gi = point_to_index[points[i]]
    s = point_summaries[i]
    if s.status == :hatched_historically
        d = Dates.value(historical_hatch_date(historical_egg_results[i]))
        median_dates[gi] = p25_dates[gi] = p75_dates[gi] = d
        m = ustrip(u"mg", historical_egg_results[i].final_state.egg_mass)
        median_mass[gi] = p25_mass[gi] = p75_mass[gi] = m
    elseif s.p25 !== nothing
        median_dates[gi] = Dates.value(hatch_date_of(s.median, issue_date))
        p25_dates[gi] = Dates.value(hatch_date_of(s.p25, issue_date))
        p75_dates[gi] = Dates.value(hatch_date_of(s.p75, issue_date))
        median_mass[gi] = mass_mg_of(s.median)
        p25_mass[gi] = mass_mg_of(s.p25)
        p75_mass[gi] = mass_mg_of(s.p75)
    end
    s.cv_hatch_days !== nothing && (cv_dates[gi] = s.cv_hatch_days)
    s.cv_hatch_mass_mg !== nothing && (cv_mass[gi] = s.cv_hatch_mass_mg)
    # historically-resolved points are a single outcome (n=1), not an ensemble
    total = s.status == :hatched_historically || startswith(string(s.status), "died_historically") ? 1 : length(members)
    frac_hatched[gi] = get(s.counts, :hatched, 0) / total
    frac_cold[gi] = get(s.counts, :cold, 0) / total
    frac_heat[gi] = get(s.counts, :heat, 0) / total
    frac_desiccation[gi] = get(s.counts, :desiccation, 0) / total
end

# Colourbars can't show date strings (Plots.jl GR limitation) -- a series
# legend is built separately (legend_bar below) and shown once, as its own
# panel, rather than attached to any one subplot here. `lo`/`hi` are shared
# across median/p25/p75 so the same colour always means the same
# date in every subplot.
const HATCH_DATE_GRADIENT = cgrad(:plasma)
function date_heatmap(title, ordinals, lo, hi)
    isempty(filter(!isnan, ordinals)) && return add_basemap!(heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude"))
    p = heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude", color=HATCH_DATE_GRADIENT, clims=(lo, hi), colorbar=false, aspect_ratio = :equal,)
    add_basemap!(p)
end

# A legend-only panel: adjacent colour swatches (a 1-row categorical
# heatmap, no gaps) with the date range as an x-tick label under each one --
# same visual convention as the APLC reference, and immune to Plots.jl's
# legend-engine auto-wrap/clipping quirks since there's no legend involved.
function legend_bar(entries)   # entries :: Vector of (color, label)
    n_entries = length(entries)
    heatmap(1:n_entries, [1], reshape(1:n_entries, 1, n_entries);
        color=cgrad(first.(entries); categorical=true), clims=(0.5, n_entries + 0.5),
        colorbar=false, yticks=false, xticks=(1:n_entries, last.(entries)),
        xlabel="", ylabel="", framestyle=:box, tickfontsize=8)
end

# Fixed colour scheme + date bounds: the APLC forecaster's own 8 dekad bands
# (11 Aug - 31 Oct) kept exactly as published, extended with finer-grained
# dekads on both sides (mid-Jul - 10 Aug early, 1 Nov - 31 Dec late) for
# resolution beyond that reference window. Not derived from this run's own
# data, so genuinely early/late hatch dates land in real dated bands rather
# than a single "too early"/"too late" catch-all.
const HATCH_DATE_BANDS = [
    # early extension -- red tones leading into the APLC's own burnt orange
    (7, 11, "11-20/Jul"), (7, 21, "21-31/Jul"), (8, 1, "01-10/Aug"),
    # unchanged APLC reference bands
    (8, 11, "11-20/Aug"), (8, 21, "21-31/Aug"), (9, 1, "01-10/Sep"), (9, 11, "11-20/Sep"),
    (9, 21, "21-30/Sep"), (10, 1, "01-10/Oct"), (10, 11, "11-20/Oct"), (10, 21, "21-31/Oct"),
    # late extension -- purple/violet tones continuing from the APLC's own navy
    (11, 1, "01-10/Nov"), (11, 11, "11-20/Nov"), (11, 21, "21-30/Nov"),
    (12, 1, "01-10/Dec"), (12, 11, "11-20/Dec"), (12, 21, "21-31/Dec"),
]
const HATCH_DATE_BAND_COLORS = [
    "#7A0C0C", "#A31515", "#C2401A",                                                 # early extension
    "#D2691E", "#F2A900", "#FFFF00", "#8CC63F", "#22B14C", "#2E9E83", "#3B7EA8", "#1B3F8B",  # APLC, unchanged
    "#3B2F8B", "#5B2E99", "#7A2EA6", "#992EB3", "#B82EC0", "#D62ECC",                # late extension
]
const HATCH_DATE_EARLY_COLOR = "#3D0000"   # safety net, beyond even the extended early range
const HATCH_DATE_LATE_COLOR = "#3D003D"    # safety net, beyond even the extended late range
const N_HATCH_BANDS = length(HATCH_DATE_BANDS)
const HATCH_DATE_PALETTE = vcat(HATCH_DATE_EARLY_COLOR, HATCH_DATE_BAND_COLORS, HATCH_DATE_LATE_COLOR)
const HATCH_DATE_LABELS = vcat("before 11/Jul", last.(HATCH_DATE_BANDS), "after 31/Dec")

# band_edges: N_HATCH_BANDS start dates + 1 final (exclusive) end date --
# the last band (21-31/Dec) ends in `year_ref+1`. searchsortedlast on this
# gives 0 (before the first band) through N_HATCH_BANDS+1 (at/after the last edge).
hatch_band_edges(year_ref) = Dates.value.(vcat(
    [Date(year_ref, m, d) for (m, d, _) in HATCH_DATE_BANDS],
    Date(year_ref + 1, 1, 1),   # exclusive end of the last band (21-31/Dec)
))
hatch_band_of(v, band_edges) = searchsortedlast(band_edges, v)

# band_edges is the same fixed reference for every subplot, so a given band
# is always the same colour and date range everywhere. legend_present (the
# legend is built separately (legend_bar) and shown once as its own panel.
function banded_date_heatmap(title, ordinals, band_edges)
    isempty(filter(!isnan, ordinals)) && return add_basemap!(heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude"))
    codes = [isnan(v) ? NaN : Float64(hatch_band_of(v, band_edges)) for v in ordinals]
    p = heatmap(lon_range, lat_range, to_heatmap_z(codes);
        title, xlabel="Longitude", ylabel="Latitude",
        color=cgrad(HATCH_DATE_PALETTE; categorical=true),
        clims=(-0.5, N_HATCH_BANDS + 1.5), colorbar=false,
        aspect_ratio = :equal)
    add_basemap!(p)
end

# `clims`, when given, is shared across a group of subplots so they're all on
# the same colour scale; `colorbar` is only enabled on one of them to avoid
# redundant repeated colorbars.
function mass_heatmap(title, vals; clims=nothing, colorbar=true)
    kw = clims === nothing ? (;) : (; clims)
    p = heatmap(lon_range, lat_range, to_heatmap_z(vals);
        title, xlabel="Longitude", ylabel="Latitude", color=cgrad(:viridis), aspect_ratio = :equal, colorbar, kw...)
    add_basemap!(p)
end
function frac_heatmap(title, vals)
    p = heatmap(lon_range, lat_range, to_heatmap_z(vals);
        title, xlabel="Longitude", ylabel="Latitude", color=cgrad(:viridis), clims=(0, 1))
    add_basemap!(p)
end

plot_title = "Lay date $oviposition_date, issue date $issue_date, diapause=$diapause"

# shared across median/p25/p75 so the same colour (band) always means the
# same date range in every subplot -- otherwise each independently picks its
# own scale and the panels aren't comparable at a glance.
shared_hatch_date_lo, shared_hatch_date_hi = extrema(filter(!isnan, vcat(median_dates, p25_dates, p75_dates)))

# 2x2 map grid + a slim legend-only strip along the bottom (APLC-style)
# instead of attaching a legend to any one subplot. A fresh @layout is built
# for each panel below -- reusing the same layout object across two plot()
# calls corrupts the second render (confirmed: only one subplot renders, in
# the wrong slot, the rest blank).

date_tick_vals = unique(round.(Int, range(shared_hatch_date_lo, shared_hatch_date_hi;
    length=min(6, length(unique(round.(Int, filter(!isnan, vcat(median_dates, p25_dates, p75_dates)))))))))
hatch_date_legend = legend_bar([
    (HATCH_DATE_GRADIENT[(v - shared_hatch_date_lo) / max(shared_hatch_date_hi - shared_hatch_date_lo, 1)], string(Date(Dates.UTD(v))))
    for v in date_tick_vals
])

hatch_date_panel = plot(
    date_heatmap("Median hatch date", median_dates, shared_hatch_date_lo, shared_hatch_date_hi),
    date_heatmap("25th percentile hatch date", p25_dates, shared_hatch_date_lo, shared_hatch_date_hi),
    date_heatmap("75th percentile hatch date", p75_dates, shared_hatch_date_lo, shared_hatch_date_hi),
    mass_heatmap("CV of hatch date (%)", cv_dates),
    hatch_date_legend;
    layout=(@layout [grid(2, 2); b{0.08h}]), size=(1000, 1150), plot_title,
)
hatch_date_path = joinpath(egg_dir, "splice_hatch_dates_$(run_tag).png")
savefig(hatch_date_panel, hatch_date_path)
println("Saved $hatch_date_path")

fixed_band_edges = hatch_band_edges(year(oviposition_date))
hatch_present_overall = sort(unique(Float64(hatch_band_of(v, fixed_band_edges))
    for v in filter(!isnan, vcat(median_dates, p25_dates, p75_dates))))
hatch_date_banded_legend = legend_bar([
    (HATCH_DATE_PALETTE[Int(code)+1], HATCH_DATE_LABELS[Int(code)+1]) for code in hatch_present_overall
])

hatch_date_banded_panel = plot(
    banded_date_heatmap("Median hatch date", median_dates, fixed_band_edges),
    banded_date_heatmap("25th percentile hatch date", p25_dates, fixed_band_edges),
    banded_date_heatmap("75th percentile hatch date", p75_dates, fixed_band_edges),
    mass_heatmap("CV of hatch date (%)", cv_dates),
    hatch_date_banded_legend;
    layout=(@layout [grid(2, 2); b{0.08h}]), size=(1000, 1150), plot_title,
)
hatch_date_banded_path = joinpath(egg_dir, "splice_hatch_dates_banded_$(run_tag).png")
savefig(hatch_date_banded_panel, hatch_date_banded_path)
println("Saved $hatch_date_banded_path")

# shared across median/p25/p75 for the same reason as the hatch-date panels
# above; cv_mass is a different quantity/scale so it keeps its own.
shared_mass_clims = extrema(filter(!isnan, vcat(median_mass, p25_mass, p75_mass)))

mass_panel = plot(
    mass_heatmap("Median egg mass at hatch (mg)", median_mass; clims=shared_mass_clims, colorbar=true),
    mass_heatmap("25th percentile egg mass at hatch (mg)", p25_mass; clims=shared_mass_clims, colorbar=false),
    mass_heatmap("75th percentile egg mass at hatch (mg)", p75_mass; clims=shared_mass_clims, colorbar=false),
    mass_heatmap("CV of egg mass at hatch (%)", cv_mass);
    layout=(2, 2), size=(1000, 1050), plot_title,
)
mass_path = joinpath(egg_dir, "splice_mass_$(run_tag).png")
savefig(mass_panel, mass_path)
println("Saved $mass_path")

mortality_panel = plot(
    frac_heatmap("Fraction hatched", frac_hatched),
    frac_heatmap("Fraction died of heat", frac_heat),
    frac_heatmap("Fraction died of cold", frac_cold),
    frac_heatmap("Fraction died of desiccation", frac_desiccation);
    layout=(2, 2), size=(1100, 1050), plot_title,
)
mortality_path = joinpath(egg_dir, "splice_mortality_$(run_tag).png")
savefig(mortality_panel, mortality_path)
println("Saved $mortality_path")
