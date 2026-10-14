# Microclimate forecast from ECMWF SEAS5 seasonal forecasts, continuing on
# from an ERA5-Land history.
#
# History: ERA5-Land, hourly, from `history_start` up to the forecast start.
# Forecast: one microclimate per SEAS5 ensemble member, from the forecast
# start for up to 215 days, each starting from the history's final soil
# temperature and moisture profile so the soil state carries across.
#
# Needs a CDS API key (ENV["CDSAPI_KEY"] or ~/.cdsapirc) for both, and the
# Terms of Use accepted (once) for ERA5-Land and for
# https://cds.climate.copernicus.eu/datasets/seasonal-original-single-levels
# `history` and `forecasts` (one output per member) keep MicroclimateMapper's
# full default outputs, ready to drive the ectotherm/endotherm demos.
#
# The first SEAS5 run queues 7 CDS requests (a few minutes); every member
# comes from those same files, so later members and re-runs are quick.

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model
using Rasters, RasterDataSources
using ZarrDatasets   # ERA5ECMWFLand reads ECMWF's Zarr stores
using Dates, Unitful, Statistics
using Plots

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

site = geocode("Bendigo, Australia")
points = [site]

# Forecast start ("now"): SEAS5 forecasts start on the 1st of each month.
# ECMWF's ERA5-Land store runs ~6 weeks behind real time, and the history must
# reach the day before the forecast starts -- so use the latest forecast
# whose previous day ERA5-Land already covers.
forecast_start = Date(2026, 8, 1)
history_start  = Date(2026, 5, 1)
forecast_days  = 215                 # SEAS5's full horizon; shorten for quicker runs
members        = 0:10                # 0 is the control; up to 0:50

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
soil_profile = example_soil_profile(depths)

plot_depth = 5.0u"cm"
plot_node = argmin(abs.(depths .- plot_depth))

function map_model(weather_source; init_source = nothing, config = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()))
    MicroMapModel(;
        micro_model = MicroModel(;
            depths, heights,
            soil_properties_model = example_soil_properties_model(),
            soil_hydraulic_model  = example_soil_hydraulic_model(),
            snow_model            = NoSnow(),
            config,
        ),
        dem_source              = CRUCL2,
        weather_source,
        init_source,
        surface_albedo_source   = 0.15,
        roughness_height_source = 0.004u"m",
        compute_terrain         = false,
    )
end

# ── history: ERA5-Land, history_start up to the day before the forecast ──

history_dates = history_start:Day(1):(forecast_start - Day(1))
println("Solving history (ERA5-Land): $(first(history_dates)) to $(last(history_dates))...")
# ERA5-Land's own soil state on `history_start` starts the history run.
@time history = solve(MicroVectorProblem(;
    model = map_model(ERA5ECMWFLand; init_source = ERA5ECMWFLand),
    points, dates = history_dates, soil_profile))

# Soil profile at the last history hour -> each forecast member's start.
history_temperature = history.soil_temperature[point = 1]
history_moisture = history.soil_moisture[point = 1]
now_soil_temperature = collect(history_temperature[Ti = size(history_temperature, Ti)])
now_soil_moisture = collect(history_moisture[Ti = size(history_moisture, Ti)])

# ── forecast: one microclimate per SEAS5 ensemble member ──

forecast_dates = forecast_start:Day(1):(forecast_start + Day(forecast_days - 1))
# One solver iteration per day: the soil state already carries over from the
# history, so this trades a little accuracy for a several-fold speed-up
# across members.
forecast_config = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture(),
    convergence = FixedIterationConvergence(1))
forecasts = map(members) do member
    println("Solving SEAS5 member $member: $(first(forecast_dates)) to $(last(forecast_dates))...")
    model = map_model(seas5_forecast(forecast_start; member); config = forecast_config)
    @time solve(MicroVectorProblem(; model, points, dates = forecast_dates, soil_profile,
        init = (; soil_temperature = now_soil_temperature, soil_moisture = now_soil_moisture)))
end

# ── plot daily means at plot_depth: history, each member, ensemble median ──

# Daily means of an hourly (Ti, depth) series at `node`, with their dates.
function daily_at(series, node)
    days = Date.(lookup(series, Ti))
    values = ustrip.(collect(series[:, node]))
    unique_days = unique(days)
    unique_days, [mean(values[days .== d]) for d in unique_days]
end
celsius(series) = uconvert.(u"°C", series)

history_days, history_T = daily_at(celsius(history_temperature), plot_node)
_, history_m = daily_at(history_moisture, plot_node)
member_T = [last(daily_at(celsius(f.soil_temperature[point = 1]), plot_node)) for f in forecasts]
member_m = [last(daily_at(f.soil_moisture[point = 1], plot_node)) for f in forecasts]
forecast_days_axis = first(daily_at(forecasts[1].soil_moisture[point = 1], plot_node))

place = split(site.display_name, ",")[1]
p1 = plot(; ylabel = "soil temperature\nat $(plot_depth) (°C)", legend = :topleft,
    title = "$place: ERA5-Land, then SEAS5 forecast from $forecast_start ($(length(members)) members)")
p2 = plot(; ylabel = "soil moisture\nat $(plot_depth) (m³/m³)", xlabel = "date", legend = false)
for (i, (T, m)) in enumerate(zip(member_T, member_m))
    plot!(p1, forecast_days_axis, T; color = :steelblue, alpha = 0.3, label = i == 1 ? "SEAS5 members" : "")
    plot!(p2, forecast_days_axis, m; color = :steelblue, alpha = 0.3)
end
plot!(p1, forecast_days_axis, vec(median(reduce(hcat, member_T); dims = 2)); color = :navy, lw = 2, label = "ensemble median")
plot!(p2, forecast_days_axis, vec(median(reduce(hcat, member_m); dims = 2)); color = :navy, lw = 2)
plot!(p1, history_days, history_T; color = :black, lw = 1.5, label = "ERA5-Land")
plot!(p2, history_days, history_m; color = :black, lw = 1.5)
vline!(p1, [forecast_start]; color = :grey, ls = :dash, label = "")
vline!(p2, [forecast_start]; color = :grey, ls = :dash)

combined = plot(p1, p2; layout = (2, 1), size = (1000, 700), link = :x)
savefig(combined, joinpath(@__DIR__, "seas5_forecast.png"))
display(combined)
