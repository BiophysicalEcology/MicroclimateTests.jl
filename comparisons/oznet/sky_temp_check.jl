# sky_temp_check.jl — TEMPORARY diagnostic: compare Julia's output.sky_temperature
# against NicheMapR's metout$TSKYC for one site, to check whether cloud-cover
# derivation (MicroclimateMapper.jl, from SILO's radiation-based cloud) is
# causing Julia's excess nighttime cooling relative to NMR. Delete when done.

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

include("config.jl")
include("utils.jl")
include("pipeline.jl")

site_name = "a1"
sim_start = Date(2007, 1, 1)
sim_end   = Date(2010, 12, 31)
auto_date_range = true

siteinfo = DataFrame(CSV.File(OZNET_SITEINFO))
row = siteinfo[findfirst(==(site_name), siteinfo.name), :]
weather_cache = Dict{Tuple{Symbol,String,Date,Date}, Any}()

prep = prepare_site(row, weather_cache; sim_start, sim_end, auto_date_range, max_sim_years)
micro_out = Microclimate.solve(prep.problem)

julia_tsky = ustrip.(u"°C", collect(micro_out.sky_temperature))
julia_tair = ustrip.(u"°C", collect(micro_out.reference_temperature))
julia_cloud = collect(micro_out.cloud_cover)

metout = DataFrame(CSV.File(joinpath(nmr_out_dir, site_name, "metout.csv")))
nmr_tsky = Float64.(metout.TSKYC)

n = min(length(julia_tsky), length(nmr_tsky))
println("Julia sky_temperature: n=$(length(julia_tsky)), range=$(extrema(julia_tsky))")
println("NMR   TSKYC:           n=$(length(nmr_tsky)), range=$(extrema(nmr_tsky))")
println("mean(Julia - NMR) over first $n hours: ", mean(julia_tsky[1:n] .- nmr_tsky[1:n]))

# Correlate the zero-Kelvin dips against air_temperature and cloud_cover at
# the exact same hours, to find which quantity is actually degenerate.
bad = findall(<=(-273.0), julia_tsky)
println("\n$(length(bad)) hour(s) at/near absolute zero sky_temperature.")
if !isempty(bad)
    k = first(bad)
    lo, hi = max(1, k - 2), min(n, k + 2)
    println("First occurrence at hour index $k. Context [$lo:$hi]:")
    println("  sky_temperature (°C): ", julia_tsky[lo:hi])
    println("  air_temperature (°C): ", julia_tair[lo:hi])
    println("  cloud_cover:          ", julia_cloud[lo:hi])
end

# Full-period plot, plus a zoomed first-two-weeks plot to see the diurnal shape.
p1 = plot(1:n, julia_tsky[1:n]; label = "Julia", color = :black, lw = 1,
    title = "Sky temperature — $site_name (full period)", ylabel = "°C")
plot!(p1, 1:n, nmr_tsky[1:n]; label = "NicheMapR", color = :blue, lw = 1, alpha = 0.7)

zoom = 1:min(24*14, n)
p2 = plot(zoom, julia_tsky[zoom]; label = "Julia", color = :black, lw = 1.5,
    title = "Sky temperature — $site_name (first 14 days)", ylabel = "°C")
plot!(p2, zoom, nmr_tsky[zoom]; label = "NicheMapR", color = :blue, lw = 1.5, alpha = 0.7)

display(plot(p1, p2; layout = (2, 1), size = (900, 700), link = :x))
