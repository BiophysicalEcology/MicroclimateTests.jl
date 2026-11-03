# barra.jl — hourly forcing from BARRA-R2, the Bureau of Meteorology's ~12 km regional
# reanalysis, for whatever a site's tower can't supply (forcing_selection.jl). Air temperature
# and wind at BARRA_REFERENCE_HEIGHT above the ground,
# specific humidity interpolated to that height between pressure levels, and surface
# radiation, rainfall and pressure; temperature and pressure taken from the grid cell's surface
# height to the site's. Fetched point by point from NCI's THREDDS NetCDF Subset Service, a
# variable-month per request, and cached in BARRA_DIR.

const BARRA_NCSS = "https://thredds.nci.org.au/thredds/ncss/grid/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1"
const BARRA_REFERENCE_HEIGHT = 100.0u"m"
const BARRA_HUMIDITY_LEVELS_HPA = (950, 925, 850)
const DRY_AIR_GAS_CONSTANT = 287.05  # J/kg/K
const STANDARD_LAPSE_RATE = 0.0065   # K/m
const GRAVITY = 9.80665              # m/s²

# (time, value) pairs of a BARRA-R2 variable over a file: hourly `var` for `month`, or the
# fixed field `var` (orog). Times are UTC: instants on the hour, hour means at the half hour.
function _barra_point(var, latitude, longitude; month=nothing)
    file = month === nothing ? "fx/$var/latest/$(var)_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1.nc" :
        "1hr/$var/latest/$(var)_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_1hr_$(Dates.format(month, "yyyymm"))-$(Dates.format(month, "yyyymm")).nc"
    url = "$BARRA_NCSS/$file?var=$var&latitude=$(ustrip(u"°", latitude))&longitude=$(ustrip(u"°", longitude))&temporal=all&accept=csv"
    fetch() = String(take!(Downloads.download(url, IOBuffer())))
    body = retry(fetch; delays=ExponentialBackOff(n=4, first_delay=5.0))()
    return map(split(strip(body), '\n')[2:end]) do line
        fields = split(line, ',')
        (time=DateTime(fields[1][1:19]), value=parse(Float64, fields[end]))
    end
end

# Specific humidity at pressure `p` (hPa) by linear interpolation in log pressure between the
# surface (`ps`, `huss`) and the pressure levels above it
function _humidity_at_pressure(p, ps, huss, level_humidity)
    anchors = sort([(ps, huss); [(level, q) for (level, q) in level_humidity if level < ps]]; rev=true)
    i = something(findfirst(a -> first(a) <= p, anchors), length(anchors))
    i == 1 && return last(anchors[1])
    (p₁, q₁), (p₂, q₂) = anchors[i - 1], anchors[i]
    return q₁ + (q₂ - q₁) * log(p / p₁) / log(p₂ / p₁)
end

# Adds to the tower's `hourly` rows a reanalysis column (forcing_selection.jl's
# reanalysis_column) for each of the `quantities` forcing_sources takes from BARRA-R2
function add_reanalysis_forcing!(hourly, resolved, quantities)
    isempty(quantities) && return hourly
    (; site_name, latitude, longitude, elevation, years) = resolved
    forcing = barra_forcing(site_name, latitude, longitude, elevation, years, site_utc_offset(site_name))
    index = Dict(t => i for (i, t) in enumerate(forcing.DateTime))
    for quantity in quantities
        hourly[!, reanalysis_column(quantity)] = [haskey(index, t) ? forcing[index[t], quantity] : missing for t in hourly.DateTime]
    end
    return hourly
end

# Hourly BARRA-R2 forcing for `years` at a site, in the tower file's units (air temperature °C,
# wind speed m/s, absolute humidity g/m³, radiation W/m², rainfall mm, pressure kPa), on the
# hours the tower's own hourly rows use: local standard time (`utc_offset`), labelled by the
# hour's start
function barra_forcing(site_name, latitude, longitude, elevation, years, utc_offset::Minute)
    path = joinpath(BARRA_DIR, "$(site_name)_$(join(years, '-')).csv")
    isfile(path) && return CSV.read(path, DataFrame)
    mkpath(BARRA_DIR)
    first_hour, last_hour = DateTime(minimum(years), 1, 1) - utc_offset, DateTime(maximum(years) + 1, 1, 1) - utc_offset
    months = Date(Dates.year(first_hour), Dates.month(first_hour)):Month(1):Date(Dates.year(last_hour), Dates.month(last_hour))
    level_names = ["hus$level" for level in BARRA_HUMIDITY_LEVELS_HPA]
    variables = ["ta100m", "ua100m", "va100m", "huss", "ps", "rsds", "rlds", "pr", level_names...]
    requests = [(variable, month) for variable in variables for month in months]
    fetched = asyncmap(((variable, month),) -> _barra_point(variable, latitude, longitude; month), requests; ntasks=16)
    series = Dict(variable => Dict{DateTime,Float64}() for variable in variables)
    for ((variable, _), points) in zip(requests, fetched), point in points
        series[variable][point.time] = point.value
    end
    grid_height = only(_barra_point("orog", latitude, longitude)).value * u"m"
    Δz = ustrip(u"m", elevation - grid_height)

    hours = first_hour:Hour(1):last_hour - Hour(1)
    # States are instants on the hour: an hour's mean is the mean of its two ends
    state(name, t) = (series[name][t] + series[name][t + Hour(1)]) / 2
    # Fluxes are hour means, labelled at the half hour
    flux(name, t) = series[name][t + Minute(30)]
    rows = map(hours) do t
        T_grid = state("ta100m", t)                                  # K, 100 m above the grid cell
        ps_grid = state("ps", t)                                     # Pa
        p_grid = ps_grid * exp(-GRAVITY * ustrip(u"m", BARRA_REFERENCE_HEIGHT) / (DRY_AIR_GAS_CONSTANT * T_grid))
        q = _humidity_at_pressure(p_grid / 100, ps_grid / 100, state("huss", t),
            [(level, state(name, t)) for (level, name) in zip(BARRA_HUMIDITY_LEVELS_HPA, level_names)])
        T = T_grid - STANDARD_LAPSE_RATE * Δz                        # at the same height above the site
        ps = ps_grid * exp(-GRAVITY * Δz / (DRY_AIR_GAS_CONSTANT * T_grid))
        p = ps * exp(-GRAVITY * ustrip(u"m", BARRA_REFERENCE_HEIGHT) / (DRY_AIR_GAS_CONSTANT * T))
        vapour_density = q * p / (DRY_AIR_GAS_CONSTANT * T * (1 + 0.608q))  # kg/m³
        (DateTime = t + utc_offset, air_temperature = T - 273.15, wind_speed = hypot(state("ua100m", t), state("va100m", t)),
         absolute_humidity = 1000vapour_density, downwelling_shortwave = flux("rsds", t), downwelling_longwave = flux("rlds", t),
         rainfall = 3600flux("pr", t), air_pressure = ps / 1000)
    end
    forcing = DataFrame(rows)
    CSV.write(path, forcing)
    return forcing
end
