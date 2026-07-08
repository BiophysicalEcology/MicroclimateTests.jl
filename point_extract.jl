using RasterDataSources
using Dates
using Rasters
using MicroclimateMapper
using PointDataSources

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

points = [geocode("Daly Waters, NT, Australia")]
lon, lat = points[1].lon, points[1].lat
date = Date(2025, 1)

# === RasterDataSources.jl + Rasters.jl: download the month's file, then index the point out of it ===
#getraster(BARRA{BARRAC2, AUST04}, :hurs; date)
Rasters.DiskArrays.allowscalar(false)
rast = Raster(BARRA{BARRAC2, AUST04}, :hurs; date, lazy=true)

@time raster_series = collect(rast[X(Near(lon)), Y(Near(lat))])
raster_times  = collect(dims(rast, Ti))

# === PointDataSources.jl: fetch the same month directly via OPeNDAP, no local download ===
month_range = (Date(2025, 1, 1), Date(2025, 1, 31))
@time nt = getpoint(BARRA{BARRAC2, AUST04}, :hurs; lon, lat, date=month_range)

point_series = nt.values
point_times  = nt.times
point_units  = nt.units
