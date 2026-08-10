# single_site.jl
#
# Flat, single-site debugging entry point — calls `run_site` directly with
# no try/catch, so errors propagate normally with a full stack trace pointing
# at the real failure line. Step through top-to-bottom (Ctrl+Enter /
# Shift+Enter per block in the VS Code Julia extension), or
# `include("single_site.jl")` from the REPL to run it all at once.
#
# All model setup and per-site logic lives in pipeline.jl/report.jl (shared
# with comparison.jl) and config.jl (shared model parameters — change
# leaf_convection_model_choice, soil_hydraulic_model, etc. there).

using Microclimate, MicroclimateMapper, Unitful, FluidProperties
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization, Measures


include("config.jl")
include("utils.jl")
include("pipeline.jl")
include("report.jl")

# ── Pick the site/year(s) to debug here ───────────────────────────────────────
site_name = "Calperum"  # "CapeTribulation", "Calperum" (:legacy), "Whroo", "Wallaby", "GWW" (:legacy), "Longreach" (:legacy), "TiTreeEast" (:legacy), "AliceSpringsMulga"
years = [2025]
max_days = nothing  # set e.g. 30 for a fast iteration loop
forcing_source = :gapfilled  # :tower, :silo (gridded forcing), or :gapfilled (tower + SILO-filled gaps)
canopy_mode = :full  # :full (MultilayerCanopy) or :legacy (NoCanopy + shade/wind/horizon approximation) -- :legacy only wired up for :tower/:gapfilled, not :silo

result = forcing_source == :silo ? run_site_silo(site_name, years; max_days) :
          forcing_source == :gapfilled ? run_site_gapfilled(site_name, years; max_days, canopy_mode) :
          run_site(site_name, years; max_days, canopy_mode)

plot_start = nothing
plot_end   = nothing
plot_start = DateTime(years[1], 1, 1, 0)
plot_end   = DateTime(years[1], 12, 30, 23)
plot_forcing = false  # false: skip the Ta/RH/Ws/Fsd/Fld/Precip forcing-sanity plots
plot_variable = nothing  # nothing plots every variable; or a short key: "Fsu","Flu","Fn","Fh","Fe","Fg","Ts","Sws", or a height-series name like "Ta_HMP_2m" -- stats are still computed/printed for all either way

stats = report_site_results(result, result.output; plot_start, plot_end, display_plots=true,
    source="$(forcing_source)_$(canopy_mode)", plot_forcing, plot_variable)
#display(stats)

# ── Vertical canopy profiles (obs vs model) for a chosen date/hour range ──────
profile_times = DateTime(years[1], 11, 15, 1):Hour(1):DateTime(years[1], 11, 15, 16)
plot_canopy_profiles_all(result, result.output, profile_times)

