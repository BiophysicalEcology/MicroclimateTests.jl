using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using DataInterpolations
using Dates
using Unitful

nearest_node(value, nodes) = argmin(abs.(nodes .- value))

# rolling median-of-3: removes isolated single-point spikes (e.g. an
# occasional solver glitch in an upstream microclimate solve) while leaving
# smooth trends untouched, since the median of 3 points including one outlier
# is always one of the two good neighbours.
median3(a, b, c) = max(min(a, b), min(max(a, b), c))
function median_filter3(series)
    n = length(series)
    n < 3 && return copy(series)
    filtered = copy(series)
    for i in 2:(n - 1)
        filtered[i] = median3(series[i-1], series[i], series[i+1])
    end
    filtered
end

# hours from the start of a solved microclimate result to a given oviposition
# date -- lets simulate_egg be re-run from any lay date against one already-
# solved `forcing`/microclimate result, without re-solving the microclimate.
oviposition_offset(oviposition_date::Date, dates) = Dates.value(oviposition_date - first(dates)) * 24.0u"hr"

# Nest-environment forcing at a fixed depth node from MicroclimateMapper output,
# mirroring demos/ectotherm.jl's underground_forcing (zero radiation, near-still
# air, everything read from the soil profile at that depth) -- extended with
# soil_water_potential, which isn't part of HeatExchange's environment vocabulary
# so it's interpolated separately, matching EnvironmentForcing's own convention.
function egg_nest_forcing(result, day_range, depth_node, environment_pars)
    n = length(day_range)
    times = (0:n-1) .* 1.0u"hr" .|> u"s"
    soil_T = result.soil_temperature[day_range, depth_node]
    env_forcing = EnvironmentForcing(times, EnvironmentalVarsVec(;
        air_temperature        = soil_T,
        sky_temperature        = soil_T,
        ground_temperature     = soil_T,
        substrate_temperature  = soil_T,
        relative_humidity      = result.soil_humidity[day_range, depth_node],
        wind_speed              = fill(0.01u"m/s", n),
        atmospheric_pressure    = result.pressure[day_range],
        zenith_angle            = result.solar_radiation.zenith_angle[day_range],
        substrate_conductivity  = result.soil_thermal_conductivity[day_range, depth_node],
        global_radiation        = fill(0.0u"W/m^2", n),
        diffuse_fraction        = result.diffuse_fraction[day_range],
        shade                   = fill(1.0, n),
    ))

    times_s = ustrip.(u"s", times)
    soil_water_potential_series = median_filter3(result.soil_water_potential[day_range, depth_node])
    swp_interp = DataInterpolations.LinearInterpolation(
        soil_water_potential_series, times_s; extrapolation=DataInterpolations.ExtrapolationType.Constant,
    )

    function forcing(t)
        environment = (; environment_pars, environment_vars=env_forcing(t))
        (; environment, soil_water_potential=swp_interp(ustrip(u"s", t)))
    end
    forcing
end
