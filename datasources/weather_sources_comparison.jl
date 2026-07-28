# Pulls a point time series for several variables from SILO, BARRA, NCEP,
# and ERA5 via PointDataSources.jl directly (no MicroclimateMapper solve
# involved) and plots them together for a sanity comparison.

using PointDataSources, RasterDataSources
using Dates, Statistics, Plots
using Zarr, ZarrDatasets, CommonDataModel
using MicroclimateMapper
using Plots

silo_email = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")
ENV["CDS_API_KEY"] = get(ENV, "CDS_API_KEY") do
    error("Set CDS_API_KEY as a user environment variable first")
end
# put your API key in ~/.julia/config/startup.jl
# if !haskey(ENV, "CDS_API_KEY")
#    ENV["CDS_API_KEY"] = "your-api-key-here"
# end

site = geocode("Mildura, Australia")
site_name = split(site.display_name, ",")[1]
lon = site.lon
lat = site.lat
date_range = (Date(2025, 1, 1), Date(2025, 12, 31))

# ERA5 point access
pointlayers(ERA5)
era5_t2m = getpoint(ERA5, :t2m; lon, lat, date=date_range)
plot(Date.(era5_t2m.times), era5_t2m.values .- 273.15; label = "ERA5 t2m", 
    title = "ERA5 t2m at $site_name", xlabel = "Date", ylabel = "Temperature (°C)")

# BARRA{...,Day} (not the Hour default) and NCEP{SurfaceFlux,2} (its own
# native cadence is 6-hourly, aggregated to daily below) for like-for-like
# daily comparison against SILO, rather than aggregating hourly data here.
barra_source = PointDataSources.BARRA{BARRAC2, AUST04, Day}
ncep_source  = PointDataSources.NCEP{SurfaceFlux, 2}
pointlayers(barra_source)
pointlayers(ncep_source)

# ERA5 point access is a package extension -- only active once Rasters.jl +
# ZarrDatasets.jl are loaded. Wrapped so the rest of the comparison still
# runs (minus ERA5) if that fails (a JSON version conflict has blocked this
# before -- see the note in src/MicroclimateMapper.jl).
era5_ok = try
    @eval using Rasters, ZarrDatasets
    true
catch e
    @warn "ERA5 unavailable (Rasters/ZarrDatasets failed to load) -- skipping" exception = e
    false
end

# Group sub-daily values by calendar date and reduce each day with `reducer`
# (`maximum`/`minimum` for instantaneous temperature extremes, `mean` for
# rate variables, `sum` for per-step accumulations like ERA5's tp/ssrd).
function daily_aggregate(times, values, reducer)
    days = Date.(times)
    udays = sort(unique(days))
    return udays, [reducer(values[days .== d]) for d in udays]
end

function compare_plot(title, unit, series...)
    p = plot(; title = "$title ($unit)")
    for (label, dates, values, style) in series
        plot!(p, dates, values; label, linestyle = style)
    end
    p
end

# --- Temperature (max/min) --------------------------------------------------

silo_tmax = getpoint(SILO, :max_temp; lon, lat, date = date_range, username = silo_email)
barra_tmax = getpoint(barra_source, :tasmax; lon, lat, date = date_range)
ncep_tmax_raw = getpoint(ncep_source, :tmax; lon, lat, date = date_range)
ncep_tmax_dates, ncep_tmax_k = daily_aggregate(ncep_tmax_raw.times, ncep_tmax_raw.values, maximum)

silo_tmin = getpoint(SILO, :min_temp; lon, lat, date = date_range, username = silo_email)
barra_tmin = getpoint(barra_source, :tasmin; lon, lat, date = date_range)
ncep_tmin_raw = getpoint(ncep_source, :tmin; lon, lat, date = date_range)
ncep_tmin_dates, ncep_tmin_k = daily_aggregate(ncep_tmin_raw.times, ncep_tmin_raw.values, minimum)

tmax_series = [
    ("SILO", silo_tmax.times, silo_tmax.values, :solid),
    ("BARRA", Date.(barra_tmax.times), barra_tmax.values .- 273.15, :dash),
    ("NCEP", ncep_tmax_dates, ncep_tmax_k .- 273.15, :dot),
]
tmin_series = [
    ("SILO", silo_tmin.times, silo_tmin.values, :solid),
    ("BARRA", Date.(barra_tmin.times), barra_tmin.values .- 273.15, :dash),
    ("NCEP", ncep_tmin_dates, ncep_tmin_k .- 273.15, :dot),
]

if era5_ok
    era5_t2m = getpoint(ERA5, :t2m; lon, lat, date = date_range)
    era5_tmax_dates, era5_tmax_k = daily_aggregate(era5_t2m.times, era5_t2m.values, maximum)
    era5_tmin_dates, era5_tmin_k = daily_aggregate(era5_t2m.times, era5_t2m.values, minimum)
    push!(tmax_series, ("ERA5", era5_tmax_dates, era5_tmax_k .- 273.15, :dashdot))
    push!(tmin_series, ("ERA5", era5_tmin_dates, era5_tmin_k .- 273.15, :dashdot))
end

p1 = compare_plot("$site_name Daily maximum temperature", "°C", tmax_series...)
p2 = compare_plot("$site_name Daily minimum temperature", "°C", tmin_series...)

# --- Rainfall ----------------------------------------------------------------

silo_rain = getpoint(SILO, :daily_rain; lon, lat, date = date_range, username = silo_email)
barra_pr = getpoint(barra_source, :pr; lon, lat, date = date_range)
ncep_prate_raw = getpoint(ncep_source, :prate; lon, lat, date = date_range)
ncep_rain_dates, ncep_prate_mean = daily_aggregate(ncep_prate_raw.times, ncep_prate_raw.values, mean)

rain_series = [
    ("SILO", silo_rain.times, silo_rain.values, :solid),
    ("BARRA", Date.(barra_pr.times), barra_pr.values .* 86400.0, :dash),
    ("NCEP", ncep_rain_dates, ncep_prate_mean .* 86400.0, :dot),
]

if era5_ok
    era5_tp = getpoint(ERA5, :tp; lon, lat, date = date_range)
    era5_rain_dates, era5_tp_sum = daily_aggregate(era5_tp.times, era5_tp.values, sum)
    push!(rain_series, ("ERA5", era5_rain_dates, era5_tp_sum .* 1000.0, :dashdot))
end

p3 = compare_plot("$site_name Daily rainfall", "mm/day", rain_series...)

# --- Radiation ---------------------------------------------------------------

silo_rad = getpoint(SILO, :radiation; lon, lat, date = date_range, username = silo_email)
barra_rsds = getpoint(barra_source, :rsds; lon, lat, date = date_range)
ncep_dswrf_raw = getpoint(ncep_source, :dswrf; lon, lat, date = date_range)
ncep_rad_dates, ncep_dswrf_mean = daily_aggregate(ncep_dswrf_raw.times, ncep_dswrf_raw.values, mean)

rad_series = [
    ("SILO", silo_rad.times, silo_rad.values, :solid),
    ("BARRA", Date.(barra_rsds.times), barra_rsds.values .* 86400.0 ./ 1e6, :dash),
    ("NCEP", ncep_rad_dates, ncep_dswrf_mean .* 86400.0 ./ 1e6, :dot),
]

if era5_ok
    era5_ssrd = getpoint(ERA5, :ssrd; lon, lat, date = date_range)
    era5_rad_dates, era5_ssrd_sum = daily_aggregate(era5_ssrd.times, era5_ssrd.values, sum)
    push!(rad_series, ("ERA5", era5_rad_dates, era5_ssrd_sum ./ 1e6, :dashdot))
end

p4 = compare_plot("$site_name Daily shortwave radiation", "MJ/m²/day", rad_series...)

display(plot(p1, p2, p3, p4; layout = (4, 1), size = (900, 900), link = :x))
