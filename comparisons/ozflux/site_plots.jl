# site_plots.jl
#
# Entry point: plots every site's stats from all_stats.csv (written by comparison.jl) with
# CairoMakie, into julia_outputs/site_plots. Needs no solve or observations. For each group of
# variables, one dot per site for each measure; then every site against every variable.
#
# Usage:
#   julia site_plots.jl              ← after comparison.jl
#   include("site_plots.jl")         ← from a Julia REPL

using CSV, DataFrames, Statistics, CairoMakie

outputs_dir = joinpath(@__DIR__, "julia_outputs")
plots_dir = joinpath(outputs_dir, "site_plots")
mkpath(plots_dir)

all_stats = filter(r -> isfinite(r.rmse), CSV.read(joinpath(outputs_dir, "all_stats.csv"), DataFrame))

const INK, MUTED, DOT, REFERENCE = "#1f2933", "#7b8794", "#2f6f8f", "#b0b7bf"

const KIND_ORDER = ["radiation", "flux", "soil_temperature", "soil_moisture", "profile_air_temperature", "profile_wind_speed", "profile_absolute_humidity"]
# The depth or height a variable's name ends in; a range sorts by its top
position_m(variable) = (m = match(r"([\d.]+)(?:-[\d.]+)? m$", variable); m === nothing ? 0.0 : parse(Float64, m[1]))
variables = unique(select(all_stats, :kind, :variable, :units))
variables = variables[sortperm([(findfirst(==(r.kind), KIND_ORDER), position_m(r.variable), r.variable) for r in eachrow(variables)]), :]
site_number = Dict(site => i for (i, site) in enumerate(sort(unique(all_stats.site))))

# `reference`: the perfect value, except normalized RMSE's 1, no better than predicting the
# observed mean. `period`: of a phase shift, which is averaged around its cycle. `reaches`: the
# phase shift's axis limits to choose from, by the sites in the plot (fit_to_sites).
measure(col, title, reference, limits, ticks; scale=identity, period=nothing, reaches=nothing) =
    (; col, title, reference, limits, ticks=(Float64.(ticks), string.(ticks)), scale, period, reaches)
ratio_measure(col, title) = measure(col, title, 1.0, (0.1, 10.0), [0.1, 0.25, 0.5, 1, 2, 4, 10]; scale=log10)
const ANNUAL_PHASE = measure(:annual_phase_shift, "Annual phase shift (days)", 0.0, (-180.0, 180.0), [-180, -90, 0, 90, 180]; period=365.25, reaches=[30, 60, 90, 120, 180])
const HOURLY_FIT = [
    measure(:pearson_correlation, "Pearson correlation", 1.0, (-0.25, 1.0), [0, 0.5, 1]),
    measure(:normalized_bias, "Normalized bias", 0.0, (-3.0, 3.0), [-3, -1, 0, 1, 3]),
    measure(:normalized_rmse, "Normalized RMSE", 1.0, (0.0, 4.0), [0, 1, 2, 4]),
    ratio_measure(:std_ratio, "Std ratio"),
]
const DAILY_AND_ANNUAL = [
    ratio_measure(:daily_range_ratio, "Daily range ratio"),
    measure(:daily_phase_shift, "Daily phase shift (h)", 0.0, (-12.0, 12.0), [-12, -6, 0, 6, 12]; period=24.0, reaches=[3, 6, 12]),
    ANNUAL_PHASE,
]
const STORAGE_AND_ANNUAL = [ratio_measure(:wetting_ratio, "Wetting ratio"), ratio_measure(:drying_ratio, "Drying ratio"), ANNUAL_PHASE]

circular_mean(x, period) = atan(mean(sin.(2π .* x ./ period)), mean(cos.(2π .* x ./ period))) * period / 2π

# A phase shift's measure with its axis reaching only as far as the plot's sites do
function fit_to_sites(m, kinds)
    m.reaches === nothing && return m
    shifts = filter(!isnan, all_stats[in.(all_stats.kind, Ref(kinds)), m.col])
    largest = maximum(abs, shifts; init=0.0)
    reach = m.reaches[something(findfirst(>=(largest), m.reaches), lastindex(m.reaches))]
    ticks = [-reach, -reach ÷ 2, 0, reach ÷ 2, reach]
    return merge(m, (limits=(-Float64(reach), Float64(reach)), ticks=(Float64.(ticks), string.(ticks))))
end

# The axis limits, with room outside `limits` for the sites beyond them
function padded_limits(m)
    m.scale === log10 && return (m.limits[1] / 1.15, m.limits[2] * 1.15)
    pad = 0.04 * (m.limits[2] - m.limits[1])
    return (m.limits[1] - pad, m.limits[2] + pad)
end

function plot_sites(file, heading, kinds, measures)
    group = filter(r -> r.kind in kinds, variables)
    n = nrow(group)
    labels = ["$(r.variable) ($(r.units))" for r in eachrow(group)]
    fig = Figure(size=(280 * length(measures) + 260, 34 * n + 190), fontsize=14, backgroundcolor=:white)
    Label(fig[0, 1:length(measures)], heading; fontsize=20, font=:bold, color=INK, halign=:left)
    for (column, m) in enumerate(fit_to_sites.(measures, Ref(kinds)))
        ax = Axis(fig[1, column]; title=m.title, titlecolor=INK, xscale=m.scale, xticks=m.ticks,
            yticks=(1:n, reverse(labels)), yticklabelsvisible=column == 1, yticksvisible=false,
            ygridcolor=(:black, 0.06), xgridcolor=(:black, 0.05),
            leftspinevisible=false, rightspinevisible=false, topspinevisible=false,
            limits=(padded_limits(m)..., 0.4, n + 0.6))
        vlines!(ax, [m.reference]; color=REFERENCE, linewidth=2)
        for (i, r) in enumerate(eachrow(group))
            y = n + 1 - i
            sites = filter(s -> s.kind == r.kind && s.variable == r.variable && !isnan(s[m.col]), all_stats)
            isempty(sites) && continue
            values = sites[!, m.col]
            ys = [y + 0.26 * sin(2.4 * site_number[site]) for site in sites.site]
            inside = m.limits[1] .<= values .<= m.limits[2]
            below, above = values .< m.limits[1], values .> m.limits[2]
            scatter!(ax, values[inside], ys[inside]; color=(DOT, 0.65), markersize=10)
            any(below) && scatter!(ax, fill(m.limits[1], count(below)), ys[below]; color=INK, marker=:ltriangle, markersize=11)
            any(above) && scatter!(ax, fill(m.limits[2], count(above)), ys[above]; color=INK, marker=:rtriangle, markersize=11)
            length(values) < 2 && continue
            centre = m.period === nothing ? median(values) : circular_mean(values, m.period)
            m.limits[1] <= centre <= m.limits[2] && linesegments!(ax, [centre, centre], [y - 0.4, y + 0.4]; color=INK, linewidth=3)
        end
    end
    Label(fig[2, 1:length(measures)],
        "Each dot is a site, and triangles are sites beyond the axis. The dark bar is the median (phase shifts: the mean around the cycle). " *
        "The grey line is a perfect match; for normalized RMSE it is no better than predicting the observed mean.";
        color=MUTED, fontsize=13, halign=:left, justification=:left, word_wrap=true, tellwidth=false)
    colgap!(fig.layout, 30)
    save(joinpath(plots_dir, file), fig; px_per_unit=2)
end

# Every site against every variable, the best-fitting sites on top
function plot_site_grid(file, heading, col, colormap, colorrange, colorbar_label)
    site_order = sort(combine(groupby(all_stats, :site), :normalized_rmse => (x -> median(filter(isfinite, x))) => :fit), :fit).site
    grid = fill(NaN, nrow(variables), length(site_order))
    for r in eachrow(all_stats)
        i = findfirst(j -> variables.kind[j] == r.kind && variables.variable[j] == r.variable, 1:nrow(variables))
        grid[i, findfirst(==(r.site), site_order)] = r[col]
    end
    fig = Figure(size=(1500, 1250), fontsize=14, backgroundcolor=:white)
    Label(fig[0, 1:2], heading; fontsize=20, font=:bold, color=INK, halign=:left)
    ax = Axis(fig[1, 1]; xticks=(1:nrow(variables), variables.variable), xticklabelrotation=π / 2,
        yticks=(1:length(site_order), site_order), yreversed=true, xgridvisible=false, ygridvisible=false)
    colors = to_colormap(colormap)
    hm = heatmap!(ax, 1:nrow(variables), 1:length(site_order), grid; colormap, colorrange, nan_color=:gray95,
        lowclip=first(colors), highclip=last(colors))
    Colorbar(fig[1, 2], hm; label=colorbar_label)
    save(joinpath(plots_dir, file), fig; px_per_unit=2)
end

const ABOVE_GROUND = ["profile_air_temperature", "profile_wind_speed", "profile_absolute_humidity"]
plot_sites("fluxes_fit.png", "Radiation and energy fluxes: hourly fit", ["radiation", "flux"], HOURLY_FIT)
plot_sites("fluxes_cycles.png", "Radiation and energy fluxes: daily and annual cycles", ["radiation", "flux"], DAILY_AND_ANNUAL)
plot_sites("soil_temperature_fit.png", "Soil temperature by depth: hourly fit", ["soil_temperature"], HOURLY_FIT)
plot_sites("soil_temperature_cycles.png", "Soil temperature by depth: daily and annual cycles", ["soil_temperature"], DAILY_AND_ANNUAL)
plot_sites("soil_moisture_fit.png", "Soil moisture by depth: hourly fit", ["soil_moisture"], HOURLY_FIT)
plot_sites("soil_moisture_cycles.png", "Soil moisture by depth: wetting, drying and annual cycle", ["soil_moisture"], STORAGE_AND_ANNUAL)
plot_sites("above_ground_fit.png", "Air temperature, wind speed and humidity by height: hourly fit", ABOVE_GROUND, HOURLY_FIT)
plot_sites("above_ground_cycles.png", "Air temperature, wind speed and humidity by height: daily and annual cycles", ABOVE_GROUND, DAILY_AND_ANNUAL)
plot_site_grid("sites_normalized_rmse.png", "Normalized RMSE at every site, best-fitting sites on top", :normalized_rmse, :Blues, (0.0, 2.0),
    "Normalized RMSE (1 = no better than predicting the observed mean)")
plot_site_grid("sites_normalized_bias.png", "Normalized bias at every site, best-fitting sites on top", :normalized_bias, Reverse(:RdBu), (-2.0, 2.0),
    "Normalized bias (red: model too high)")
println("Wrote $(length(readdir(plots_dir))) plots to $plots_dir")
