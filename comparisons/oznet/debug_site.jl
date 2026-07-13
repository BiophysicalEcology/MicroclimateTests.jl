# debug_site.jl
#
# Flat, single-site debugging entry point — calls `run_site` directly with
# no try/catch, so errors propagate normally with a full stack trace pointing
# at the real failure line. Step through top-to-bottom (Ctrl+Enter /
# Shift+Enter per block in the VS Code Julia extension), or
# `include("debug_site.jl")` from the REPL to run it all at once.
#
# All model setup and per-site logic lives in pipeline.jl (shared with
# comparison.jl) and config.jl (shared model/soil-source/weather-source
# parameters — change `weather_source_choice`, `pedotransfer_model_choice`,
# etc. there to compare variants).

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

include("config.jl")
include("utils.jl")
include("pipeline.jl")

# ── Pick the site to debug here ───────────────────────────────────────────────
site_name = "a1"

sim_start = Date(2007, 1, 1)
sim_end   = Date(2010, 12, 31)
auto_date_range = true

plot_start = nothing
plot_end   = nothing

siteinfo = DataFrame(CSV.File(OZNET_SITEINFO))
row = siteinfo[findfirst(==(site_name), siteinfo.name), :]

weather_cache = Dict{Tuple{Symbol,String,Date,Date}, Any}()

@time site_rows = run_site(row, weather_cache;
    sim_start, sim_end, auto_date_range, max_sim_years, plot_start, plot_end)
