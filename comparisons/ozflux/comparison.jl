# comparison.jl
#
# Batch entry point: runs the canopy model against every requested
# (site, years) combination, writes stats+plots per site-year (report.jl),
# and prints a cross-site summary. All model setup and per-site logic lives
# in pipeline.jl/report.jl (shared with single_site.jl) and config.jl.
#
# Usage:
#   julia comparison.jl              ← runs SITE_YEARS defined below
#   include("comparison.jl")         ← from a Julia REPL

using Microclimate, MicroclimateMapper, Unitful, FluidProperties
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

include("config.jl")
include("utils.jl")
include("pipeline.jl")
include("report.jl")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  RUN CONFIGURATION                                                       ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# Each entry: site name => one or more year-groups to run separately (each
# group becomes one MicroProblem run + one report). :all_years runs every
# available year for that site as a single combined run.
const SITE_YEARS = Dict(
    "CapeTribulation"    => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["CapeTribulation"]],
    "Calperum"           => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["Calperum"]],
    "Whroo"              => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["Whroo"]],
    "Wallaby"            => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["Wallaby"]],
    "GWW"                => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["GWW"]],
    "Longreach"          => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["Longreach"]],
    "TiTreeEast"         => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["TiTreeEast"]],
    "AliceSpringsMulga"  => [[y] for y in discover_site_years(OZFLUX_DATA_DIR)["AliceSpringsMulga"]],
)

# Which forcing mode(s) to run per site-year -- :tower and :silo together
# gives the gridded-vs-point-observed comparison directly; :gapfilled runs a
# third, SILO-gap-filled tower series. Each mode adds a full solve per
# site-year (SILO/gapfilled also mean a live weather fetch, cached to disk).
const FORCING_MODES = [:tower]  # any of :tower, :silo, :gapfilled

# canopy_mode(s) to run alongside each forcing mode -- :full (MultilayerCanopy)
# vs :legacy (NoCanopy + shade/wind/horizon approximation, see pipeline.jl) to
# quantify what the full canopy model buys. :legacy only applies to
# :tower/:gapfilled (not wired up for :silo).
const CANOPY_MODES = [:full]  # any of :full, :legacy

plot_start = nothing
plot_end   = nothing
plot_forcing = true  # false: skip the Ta/RH/Ws/Fsd/Fld/Precip forcing-sanity plots (every site-year, so this adds up in a batch run)
plot_variable = nothing  # nothing plots every variable; or a short key ("Fh", "Ts", "Ta_HMP_2m", ...) -- stats are still computed for all either way

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN LOOP                                                               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

_run_site(mode, site_name, years, canopy_mode) =
    mode == :silo       ? run_site_silo(site_name, years) :
    mode == :gapfilled  ? run_site_gapfilled(site_name, years; canopy_mode) :
    run_site(site_name, years; canopy_mode)

all_stats = DataFrame()
runs = [(site, years) for (site, year_groups) in SITE_YEARS for years in year_groups
        if haskey(SITE_LEAF_AREA_INDEX, site)]
combos = [(mode, cm) for mode in FORCING_MODES for cm in (mode == :silo ? [:full] : CANOPY_MODES)]
println("Running $(length(runs)) site-year(s) x $(length(combos)) mode combo(s).")

for (site_name, years) in runs, (mode, canopy_mode) in combos
    try
        result = _run_site(mode, site_name, years, canopy_mode)
        stats = report_site_results(result, result.output; plot_start, plot_end,
            make_plots=save_outputs, display_plots=false, source="$(mode)_$(canopy_mode)", plot_forcing, plot_variable)
        append!(all_stats, stats)
    catch err
        @warn "Site $site_name years=$years mode=$mode canopy_mode=$canopy_mode failed: $err"
    end
end

# ── Cross-site/year summary (kept per forcing source) ───────────────────────
if nrow(all_stats) > 0
    println("\n" * "="^72)
    println("SUMMARY — $(length(runs)) site-year(s)")
    println("="^72)
    for source in unique(all_stats.source), kind in unique(all_stats.kind)
        sub = filter(r -> r.source == source && r.kind == kind, all_stats)
        nrow(sub) == 0 && continue
        println("\n$kind (forcing=$source):")
        println(rpad("Variable", 30) * "  mean r    mean RMSE   mean bias      n")
        println("-"^72)
        for v in unique(sub.variable)
            rows = filter(r -> r.variable == v, sub)
            valid = filter(r -> isfinite(r.r), rows)
            isempty(valid) && continue
            @printf("  %-28s  %+.3f    %8.4g    %+7.4g   %5d\n", v,
                mean(valid.r), mean(valid.rmse), mean(valid.bias), sum(valid.n))
        end
    end
end
