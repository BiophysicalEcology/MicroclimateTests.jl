# soiltexture_nodata_check.jl — small targeted check: does a single SoilGrids
# raster have a declared CRS other than lon/lat, and does a plain lon/lat
# point query find valid data there via Rasters.extract (CRS-aware) vs raw
# X/Y indexing (assumes query coords already match the raster's own CRS)?
#
# Usage: julia --project=. soiltexture_nodata_check.jl [lon] [lat]
# Defaults to SNOTEL/SCAN site 2184 (Ford Dry Lake, CA).

using MicroclimateMapper
using Rasters, RasterDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Statistics

ENV["RASTERDATASOURCES_PATH"] = "Z:"
ENV["GDAL_HTTP_TIMEOUT"] = "15"
ENV["GDAL_HTTP_CONNECTTIMEOUT"] = "10"

args = ARGS
lon = length(args) >= 1 ? parse(Float64, args[1]) : -115.0976
lat = length(args) >= 2 ? parse(Float64, args[2]) : 33.6547

depth_bins = collect(RasterDataSources.depths(SoilGrids))
path = getraster(SoilGrids, :clay; depth = depth_bins, quantile = "mean")[1]
r = Raster(path; name = :clay, lazy = true)

println("CRS: ", crs(r))
println("raster X bounds: ", extrema(dims(r, X)))
println("raster Y bounds: ", extrema(dims(r, Y)))

println("\nraw X/Y Near-indexing with plain lon/lat ($lon, $lat):")
println("  ", r[X(Near(lon)), Y(Near(lat))])

println("\nRasters.extract (CRS-aware) at the same lon/lat:")
println("  ", extract(r, [GIW.Point((lon, lat))]))

# Real production crop path (post-fix), one variable/depth only.
b = 0.05
area = Extent(X = (lon - b, lon + b), Y = (lat - b, lat + b))
var = MicroclimateMapper.texture_variables(SoilGrids)[2]  # :clay
println("\n_texture_values_from_paths([path], Extent(±$(b)°)), clay only:")
println("  ", MicroclimateMapper._texture_values_from_paths([path], area, var))
