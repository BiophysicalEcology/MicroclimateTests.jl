# forcing_selection.jl — which measurement drives each forcing variable at a site.
#
# Air temperature and humidity come from the highest height above the canopy where the tower
# measures both; that height is the reference height. Wind comes from the highest anemometer
# above the canopy, scaled to the reference height by the log wind profile. Shortwave and longwave come
# from the highest radiometer above the canopy, rainfall and pressure from any gauge. Among
# sensors at one height a single sensor is preferred to the site's composite. A sensor counts
# only when it has data in the run years and data_corrections.jl doesn't mark it faulty.
# Whatever the tower can't supply comes from the BARRA-R2 reanalysis (barra.jl); pressure falls
# back to the standard atmosphere at the site's elevation (sites.jl's site_pressure).

# The hourly column a forcing variable comes from when it is taken from the reanalysis
reanalysis_column(quantity) = "reanalysis_$quantity"
is_reanalysis_column(name) = startswith(name, "reanalysis_")

# Whether OzFlux variable `name` holds a sensor's readings of the quantity with code `code`: the
# site composite (Ta; "AH" matches Ah too), a named instrument (Ta_HMP_2m) or a numbered
# replicate (Ta_1), not a derived series (Fsd_Diffuse) or statistics of readings (Ta_HMP_2m_Sd)
function measures(name, code)
    parts = split(name, '_')
    lowercase(first(parts)) == lowercase(code) || return false
    length(parts) == 1 && return true
    is_sensor = haskey(OZFLUX_INSTRUMENT_CODES, parts[2]) || all(isdigit, parts[2])
    return is_sensor && !(last(parts) in OZFLUX_STATISTIC_SUFFIXES)
end

is_site_composite(name) = !occursin('_', name)

# Ordering among sensors of one quantity at one height: a single sensor before the site's
# composite, then by name
preference(name) = (is_site_composite(name), name)

# The instrument code in an OzFlux variable name, if it is a known one
function instrument_code(name)
    parts = split(name, '_')
    return length(parts) >= 2 && haskey(OZFLUX_INSTRUMENT_CODES, parts[2]) ? String(parts[2]) : nothing
end

measures_sonic_temperature(name) =
    (code = instrument_code(name); code !== nothing && OZFLUX_INSTRUMENT_CODES[code].measures_sonic_temperature)

# Heights (m) a variable measures at: its own `height` attribute, else the one height of the
# other variables from the same instrument, by name code or `instrument` attribute (a sonic
# anemometer's temperature carries no height, its wind and fluxes do)
function sensor_heights(ds, name)
    heights = attribute_heights(ds, name)
    isempty(heights) || return heights
    code, instrument = instrument_code(name), get(ds[name].attrib, "instrument", "")
    same_instrument(other) = (code !== nothing && instrument_code(other) == code) ||
        (!isempty(instrument) && get(ds[other].attrib, "instrument", "") == instrument)
    sibling_heights = unique(h for other in keys(ds) if other != name && same_instrument(other)
                             for h in attribute_heights(ds, other))
    return length(sibling_heights) == 1 ? sibling_heights : Float64[]
end

# Wind speed at `height` from wind speed measured at `measured_height`, over vegetation
# `canopy_height` tall: the neutral log profile with zero-plane displacement 0.65 h and roughness
# length 0.1 h (Campbell & Norman 1998, eqs. 5.2-5.3)
function wind_at_height(wind_speed, measured_height, height, canopy_height)
    displacement, roughness = 0.65 * canopy_height, 0.1 * canopy_height
    return wind_speed * log((height - displacement) / roughness) / log((measured_height - displacement) / roughness)
end

"""
    forcing_sources(site_name, years, canopy_height)

The hourly column each forcing variable is taken from at a site, the reference height, whether
the air temperature is a sonic temperature and the humidity absolute, and the quantities taken
from the reanalysis.
"""
function forcing_sources(site_name, years, canopy_height)
    corrections = site_corrections(site_name)
    codes = OZFLUX_QUANTITY_CODES
    NCDataset(ozflux_site_path(OZFLUX_DATA_DIR, site_name)) do ds
        rows = findall(t -> Dates.year(t) in years, ds["time"][:])
        has_data(name) = (values = _read_clean_var(ds, name, rows, corrections); values !== nothing && any(!ismissing, values))
        readings(code) = sort(filter(name -> measures(name, code) && !has_correction(corrections, FaultySensor, name) &&
            has_data(name), collect(keys(ds))); by=preference)
        # (name, height) of every reading above the canopy, best first at each height
        above_canopy(code) = [(; name, height) for name in readings(code) for height in sensor_heights(ds, name)
                              if height * u"m" > canopy_height]
        highest(sensors) = isempty(sensors) ? nothing : first(filter(s -> s.height == maximum(last, sensors), sensors)).name
        best_at(sensors, height) = first(filter(s -> s.height == height, sensors)).name

        reanalysis = Symbol[]
        tower_or_reanalysis(name, quantity) = name === nothing ? (push!(reanalysis, quantity); reanalysis_column(quantity)) : name

        temperature, wind = above_canopy(codes.air_temperature), above_canopy(codes.wind_speed)
        relative, absolute = above_canopy(codes.relative_humidity), above_canopy(codes.absolute_humidity)
        shared_heights = intersect(last.(temperature), last.([relative; absolute]))
        if isempty(shared_heights)
            reference_height = wind_height = BARRA_REFERENCE_HEIGHT
            air_temperature = tower_or_reanalysis(nothing, :air_temperature)
            humidity = tower_or_reanalysis(nothing, :absolute_humidity)
            wind_speed = tower_or_reanalysis(nothing, :wind_speed)
            humidity_is_absolute = true
        else
            height = maximum(shared_heights)
            reference_height = height * u"m"
            air_temperature = best_at(temperature, height)
            humidity_is_absolute = !any(s -> s.height == height, relative)
            humidity = best_at(humidity_is_absolute ? absolute : relative, height)
            # Wind from the highest anemometer above the canopy, scaled to the reference height
            wind_speed = tower_or_reanalysis(highest(wind), :wind_speed)
            wind_height = is_reanalysis_column(wind_speed) ? BARRA_REFERENCE_HEIGHT :
                first(s.height for s in wind if s.name == wind_speed) * u"m"
        end
        pressure = readings(codes.air_pressure)
        rainfall = readings(codes.rainfall)
        return (; reference_height, air_temperature, humidity, wind_speed, wind_height,
            downwelling_shortwave = tower_or_reanalysis(highest(above_canopy(codes.downwelling_shortwave)), :downwelling_shortwave),
            downwelling_longwave = tower_or_reanalysis(highest(above_canopy(codes.downwelling_longwave)), :downwelling_longwave),
            rainfall = tower_or_reanalysis(isempty(rainfall) ? nothing : first(rainfall), :rainfall),
            air_pressure = isempty(pressure) ? nothing : first(pressure),
            temperature_is_sonic = measures_sonic_temperature(air_temperature), humidity_is_absolute, reanalysis)
    end
end

# One line per forcing variable and where it comes from, for the run log
function print_forcing_sources(forcing)
    println("  Reference height $(forcing.reference_height), wind measured at $(forcing.wind_height)")
    for quantity in (:air_temperature, :humidity, :wind_speed, :downwelling_shortwave, :downwelling_longwave, :rainfall, :air_pressure)
        source = something(getproperty(forcing, quantity), "standard atmosphere at the site's elevation")
        println("  Forcing $quantity from $source")
    end
end

# The tower variables a site's forcing reads
tower_forcing_names(sources) = String[name for name in (sources.air_temperature, sources.humidity, sources.wind_speed,
    sources.downwelling_shortwave, sources.downwelling_longwave, sources.rainfall, sources.air_pressure)
    if name !== nothing && !is_reanalysis_column(name)]

# Air temperature (°C) from a sonic temperature (°C): sound travels faster in humid air, so a
# sonic anemometer reads warm by the factor (1 + 0.51 q) for specific humidity q (kg/kg)
air_temperature_from_sonic_temperature(sonic_temperature, specific_humidity) =
    (sonic_temperature + 273.15) / (1 + 0.51 * specific_humidity) - 273.15

# Specific humidity (kg/kg) from vapour density (g/m³) at air temperature (°C) and pressure (kPa)
specific_humidity(vapour_density, temperature, pressure) =
    vapour_density / 1000 * DRY_AIR_GAS_CONSTANT * (temperature + 273.15) / (pressure * 1000)
