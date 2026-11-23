# Entry point, after comparison.jl: plots every site's stats from the outputs folder's all_stats.csv into its site_plots/

using CSV, DataFrames, Statistics, CairoMakie

outputs_dir = get(ENV, "OZFLUX_OUTPUTS_DIR", joinpath(@__DIR__, "outputs", "default"))
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

# Soil values sit below ground and air values above it, so a quantity keeps its color across the two
const QUANTITIES = [
    (label="Temperature (soil and air)", kinds=["soil_temperature", "profile_air_temperature"], color="#eb6834", dark="#9c3a12", marker=:circle),
    (label="Water (soil moisture and absolute humidity)", kinds=["soil_moisture", "profile_absolute_humidity"], color="#2a78d6", dark="#134a8e", marker=:rect),
    (label="Wind speed", kinds=["profile_wind_speed"], color="#1baf7a", dark="#0b6444", marker=:diamond),
]
const PROFILE_KINDS = reduce(vcat, q.kinds for q in QUANTITIES)
const PROFILE_MEASURES = [HOURLY_FIT[1:3]..., HOURLY_FIT[4], DAILY_AND_ANNUAL[1], STORAGE_AND_ANNUAL[1], STORAGE_AND_ANNUAL[2], DAILY_AND_ANNUAL[2], ANNUAL_PHASE]
const PROFILE_TICKS = [-1.5, -0.5, -0.2, -0.1, -0.05, -0.02, 0, 0.2, 0.5, 1, 2, 5, 10, 20, 45]

# A variable's top and bottom (m) from its name, as heights above ground: depths are negative
function height_range(kind, variable)
    m = match(r"([\d.]+)(?:-([\d.]+))? m$", variable)
    top = parse(Float64, m[1])
    bottom = m[2] === nothing ? top : parse(Float64, m[2])
    return startswith(kind, "soil") ? (-top, -bottom) : (top, bottom)
end

# Every depth and height on one axis, one panel per measure. The axis is logarithmic away from the
# ground, on each side from its own scale: to show centimetres of soil and tens of metres of air
# together, without the empty space between them that a single scale would leave.
const SOIL_SCALE_M, AIR_SCALE_M = 0.01, 0.1
height_axis(y) = y >= 0 ? log10(1 + y / AIR_SCALE_M) : -log10(1 - y / SOIL_SCALE_M)
height_axis_inverse(a) = a >= 0 ? AIR_SCALE_M * (10^a - 1) : -SOIL_SCALE_M * (10^-a - 1)

function plot_profiles(file, heading)
    fig = Figure(size=(1500, 1600), fontsize=14, backgroundcolor=:white)
    Label(fig[0, 1:3], heading; fontsize=20, font=:bold, color=INK, halign=:left)
    profile = filter(r -> r.kind in PROFILE_KINDS, all_stats)
    ranges = height_range.(profile.kind, profile.variable)
    profile.top, profile.bottom = first.(ranges), last.(ranges)
    for (i, m) in enumerate(fit_to_sites.(PROFILE_MEASURES, Ref(PROFILE_KINDS)))
        row, column = fldmod1(i, 3)
        ax = Axis(fig[row + 1, column]; title=m.title, titlecolor=INK, xscale=m.scale, xticks=m.ticks,
            yscale=CairoMakie.Makie.ReversibleScale(height_axis, height_axis_inverse; limits=(-2.0f0, 2.0f0), interval=(-Inf32, Inf32)), yticks=(PROFILE_TICKS, string.(PROFILE_TICKS)), yticklabelsvisible=column == 1,
            ylabel=column == 1 ? "Height above ground (m), negative below" : "", ylabelcolor=MUTED,
            ygridcolor=(:black, 0.06), xgridcolor=(:black, 0.05), rightspinevisible=false, topspinevisible=false,
            limits=(padded_limits(m)..., -2.5, 70.0))
        hlines!(ax, [0.0]; color=INK, linewidth=1)
        vlines!(ax, [m.reference]; color=REFERENCE, linewidth=2)
        for q in QUANTITIES, position in groupby(filter(r -> r.kind in q.kinds && !isnan(r[m.col]), profile), [:kind, :variable])
            values = position[!, m.col]
            shown = clamp.(values, m.limits...)
            ys = (position.top .+ position.bottom) ./ 2
            inside = shown .== values
            # A depth range is its vertical extent, at its site's value
            position.top[1] == position.bottom[1] || linesegments!(ax, repeat(shown; inner=2),
                vec([position.top position.bottom]'); color=(q.color, 0.45), linewidth=2)
            scatter!(ax, shown[inside], ys[inside]; color=(q.color, 0.7), marker=q.marker, markersize=10)
            all(inside) || scatter!(ax, shown[.!inside], ys[.!inside]; color=q.dark,
                marker=[v < m.limits[1] ? :ltriangle : :rtriangle for v in values[.!inside]], markersize=11)
            length(values) < 2 && continue
            centre = m.period === nothing ? median(values) : circular_mean(values, m.period)
            m.limits[1] <= centre <= m.limits[2] && scatter!(ax, [centre], [ys[1]]; color=q.dark, marker=:vline, markersize=20)
        end
    end
    Legend(fig[1, 1:3], [MarkerElement(color=q.color, marker=q.marker, markersize=12) for q in QUANTITIES], [q.label for q in QUANTITIES];
        orientation=:horizontal, framevisible=false, labelcolor=INK, tellwidth=false, halign=:left)
    Label(fig[5, 1:3],
        "Each mark is a site; a vertical line through it is the depth range that sensor averages, and triangles are sites beyond the axis. " *
        "The tick in the darker color is the median (phase shifts: the mean around the cycle). The grey line is a perfect match; for normalized RMSE it is no better than predicting the observed mean. " *
        "Heights are above ground, not relative to each site's canopy.";
        color=MUTED, fontsize=13, halign=:left, justification=:left, word_wrap=true, tellwidth=false)
    colgap!(fig.layout, 24)
    save(joinpath(plots_dir, file), fig; px_per_unit=2)
end

plot_profiles("profiles.png", "Soil and air by depth and height: every measure at every site")
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
