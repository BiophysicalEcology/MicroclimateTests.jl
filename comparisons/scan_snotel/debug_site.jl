# debug_site.jl
#
# Flat, single-site debugging entry point — calls `run_site` directly with
# no try/catch, so errors propagate normally with a full stack trace pointing
# at the real failure line. Step through top-to-bottom (Ctrl+Enter /
# Shift+Enter per block in the VS Code Julia extension), or
# `include("debug_site.jl")` from the REPL to run it all at once.
#
# All model setup and per-site logic lives in pipeline.jl (shared with
# comparison.jl) and config.jl (shared model/soil/snow parameters — change
# `infiltration_algorithm_choice`, `weather_source_choice`, etc. there to
# compare variants).

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

include("config.jl")
include("utils.jl")
include("pipeline.jl")

# ── Pick the site to debug here ───────────────────────────────────────────────
site_num = 2003

sim_start = Date(2013, 1, 1)
sim_end   = Date(2013, 12, 31)
auto_date_range = true
max_sim_years   = 5

plot_start = nothing
plot_end   = nothing

meta_all = DataFrame(CSV.File(joinpath(@__DIR__, "Map metadata export.csv")))
meta_all[!, :ID] = strip.(string.(meta_all[!, :ID]))

weather_cache = Dict{Tuple{Int,Date,Date}, Any}()

micro_model, soil_profile, mapper_model = build_micro_model()

@time site_rows = run_site(site_num, meta_all, weather_cache, micro_model, soil_profile, mapper_model;
    sim_start, sim_end, auto_date_range, max_sim_years, plot_start, plot_end)
