#using NCDatasets
using Rasters
using Plots
using Dates

# ds = NCDataset(joinpath(@__DIR__,"data", "CapeTribulation_2018_L3.nc"));
# attributes = ds.attrib
# @show keys(ds)
# time = ds["time"]
# Ta = replace(collect(ds["Ta"][1, 1, :]), -9999 => missing)
# plot(time, Ta)

ds = RasterStack(joinpath(@__DIR__,"data", "CapeTribulation_2018_L3.nc"))

time = lookup(ds, Ti)

air_temperature = replace(vec(collect(ds[:Ta])), -9999.0 => missing)
absolute_humidity = replace(vec(collect(ds[:Ah])), -9999.0 => missing)
downwelling_shortwave_radiation = replace(vec(collect(ds[:Fsd])), -9999.0 => missing)
soil_temperature_10cma = replace(vec(collect(ds[:Ts_10cma])), -9999.0 => missing)
soil_moisture = replace(vec(collect(ds[:Sws])), -9999.0 => missing)
soil_moisture_10cma = replace(vec(collect(ds[:Sws_10cma])), -9999.0 => missing)
soil_moisture_10cmb = replace(vec(collect(ds[:Sws_10cmb])), -9999.0 => missing)
soil_moisture_10cmc = replace(vec(collect(ds[:Sws_10cmc])), -9999.0 => missing)
soil_moisture_75cma = replace(vec(collect(ds[:Sws_75cma])), -9999.0 => missing)
soil_moisture_75cmb = replace(vec(collect(ds[:Sws_75cmb])), -9999.0 => missing)
soil_moisture_75cmc = replace(vec(collect(ds[:Sws_75cmc])), -9999.0 => missing)
soil_moisture_150cma = replace(vec(collect(ds[Symbol("Sws_1.5mc")])), -9999.0 => missing)
soil_moisture_150cmb = replace(vec(collect(ds[Symbol("Sws_1.5mb")])), -9999.0 => missing)
soil_moisture_150cmc = replace(vec(collect(ds[Symbol("Sws_1.5mc")])), -9999.0 => missing)
 
plot(time, absolute_humidity)
plot(time, air_temperature)
plot(time, downwelling_shortwave_radiation)

plot(time, soil_temperature_10cma)

plot(time, soil_moisture_10cma)
plot!(time, soil_moisture_10cmb)
plot!(time, soil_moisture_10cmc)

plot(time, soil_moisture_75cma)
plot!(time, soil_moisture_75cmb)
plot!(time, soil_moisture_75cmc)

plot(time, soil_moisture_150cma)
plot!(time, soil_moisture_150cmb)
plot!(time, soil_moisture_150cmc)


