# mapper.jl — everything that needs MicroclimateMapper: SILO gridded forcing
# (run_site_silo, run_site_gapfilled, run_site_daily) and live SLGA soil
# (soil_source :slga/:slga_uniform). Include after config.jl/utils.jl/
# pipeline.jl, with MicroclimateMapper, Rasters, RasterDataSources,
# PointDataSources, GeoInterface and Serialization loaded. The tower-forced
# run_site path doesn't need it.
#
# MicroclimateMapper still builds Microclimate's old MicroModel/MicroInputs
# (site and soil_profile on MicroInputs), so the fetches here need it ported
# before they run against the current Microclimate.

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", joinpath(homedir(), "Spatial_Data"))

# ── SILO weather source (secondary comparison + forcing gap-fill) ───────────
weather_source_choice = SILO
ENV["SILO_EMAIL"] = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")
use_opendap_points = true  # PointDataSources.jl OPeNDAP point query, not a whole-raster download
dem_source_choice = SRTM
compute_terrain_choice = false

reuse_weather     = true
cache_weather     = true
weather_cache_dir = joinpath(@__DIR__, "weather_cache")

# ── Soil texture (live SLGA fetch + pedotransfer, like oznet) ───────────────
# SLGA has no point-query API -- a uniform profile is derived from a small
# bounding box around each site's coordinates (same approach as oznet).
pedotransfer_model_choice = CosbyMultivariate()
soil_area_buffer_deg = 0.05

# Live SLGA fetch + pedotransfer (same as oznet/pipeline.jl) for
# _fetch_soil_profile (pipeline.jl). Rescales to the site's own metadata
# BulkDensity attribute when present -- SLGA vs a real point measurement, not
# a substitute for it (SLGA is a 90 m grid mean).
function _slga_soil_profile(latitude, longitude, depths, global_attrib)
    b = soil_area_buffer_deg
    lon_dd, lat_dd = ustrip(u"°", longitude), ustrip(u"°", latitude)
    area = Extent(X = (lon_dd - b, lon_dd + b), Y = (lat_dd - b, lat_dd + b))
    soil_profile = build_soil_profile(SLGA, area; depths, pedotransfer_model = pedotransfer_model_choice).soil_profile
    if haskey(global_attrib, "BulkDensity")
        # Metadata value is a real point measurement (assumed near-surface) vs
        # SLGA's 90 m-grid mean -- rescale the whole profile by their ratio at
        # the surface node rather than overwrite it outright, preserving
        # SLGA's own depth trend. Doesn't touch the hydraulic parameters
        # already derived from the pre-scaled bulk density via pedotransfer
        # (texture, not bulk density, dominates those) -- bulk density itself
        # (thermal capacity/conductivity) is the only thing this corrects.
        meta_bd = parse(Float64, string(global_attrib["BulkDensity"])) * u"kg/m^3"
        if meta_bd > 500.0u"kg/m^3" # at least one site has erroneous value
            slga_bd = soil_profile.bulk_density[5]
            scale = meta_bd / slga_bd
            @printf("  BulkDensity: metadata=%.0f kg/m^3  SLGA(surface)=%.0f kg/m^3  (scaling profile by %.2fx)\n",
                ustrip(u"kg/m^3", meta_bd), ustrip(u"kg/m^3", slga_bd), scale)
            soil_profile.bulk_density .*= scale
        end
    end
    return soil_profile
end

# ── SILO gridded-forcing fetch (secondary comparison + gap-fill donor) ──────
# Mirrors oznet/pipeline.jl's weather-fetch block (MicroVectorProblem +
# MicroclimateMapper.init): a single-point fetch/synthesis returning
# .site/.environment_hourly/.environment_daily/.environment_minmax, ready to
# hand straight to MicroInputs. Cached on disk (weather_cache_dir) since a
# fetch is a live network call.
if use_opendap_points
    MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()
end

const _SILO_CACHE = Dict{Tuple{String,Date,Date,Float64},Any}()

# reference_height: the real site model's last(heights) -- MicroclimateMapper's
# assemble_weather! power-law-corrects SILO's native 10m wind to
# maximum(micro_model.heights), so the throwaway model built here must top out
# at the SAME height the actual site model will treat as its reference height,
# or the gap-fill/SILO-run wind ends up scaled to the wrong height (silently
# defaults to 2m otherwise -- see DEFAULT_HEIGHTS_M).
function fetch_silo_forcing(site_name, latitude, longitude, sim_start::Date, sim_end::Date, depths, reference_height)
    ref_m = ustrip(u"m", reference_height)
    wc_key = (site_name, sim_start, sim_end, ref_m)
    disk_path = joinpath(weather_cache_dir, "SILO_$(site_name)_$(sim_start)_$(sim_end)_ref$(ref_m)m.jls")
    if !haskey(_SILO_CACHE, wc_key) && reuse_weather && isfile(disk_path)
        _SILO_CACHE[wc_key] = deserialize(disk_path)
        println("  Loaded SILO forcing from disk cache: $disk_path")
    end
    if !haskey(_SILO_CACHE, wc_key)
        println("  Fetching SILO forcing ($sim_start to $sim_end, wind reference height $(ref_m) m)...")
        micro_model = MicroModel(; hours = collect(0.0:1:23.0), depths, heights = [0.01, ref_m] .* u"m",
            soil_properties_model = soil_properties_model_choice, soil_hydraulic_model = soil_hydraulic_model(site_name),
            snow_model = NoSnow())
        mapper_model = MicroMapModel(; micro_model, dem_source = dem_source_choice, weather_source = weather_source_choice,
            surface_albedo_source = site_albedo(site_name), roughness_height_source = 0.004u"m", compute_terrain = compute_terrain_choice)
        point = GIW.Point((ustrip(u"°", longitude), ustrip(u"°", latitude)))
        # Placeholder soil_profile/init -- only sizes the fetch, doesn't inform
        # it (weather forcing doesn't depend on soil properties).
        vec_problem = MicroVectorProblem(;
            model = mapper_model, points = [point], dates = sim_start:Day(1):sim_end,
            soil_profile = SoilProfile(; bulk_density = fill(1.3u"Mg/m^3", length(depths)),
                mineral_density = fill(2.56u"Mg/m^3", length(depths)), mineral_conductivity = fill(1.25u"W/m/K", length(depths)),
                mineral_heat_capacity = fill(870.0u"J/kg/K", length(depths)),
                hydraulics = CampbellHydraulicProfile(; air_entry_water_potential = fill(1.0u"J/kg", length(depths)),
                    saturated_hydraulic_conductivity = fill(0.001u"kg*s/m^3", length(depths)),
                    campbell_b_parameter = fill(4.5, length(depths)), root_density = fill(1.0u"m/m^3", length(depths)))),
            init = (; soil_moisture = fill(0.2, length(depths))),
        )
        weather_result = MicroclimateMapper.init(vec_problem)
        worker = take!(weather_result.cache_pool)
        fetched = weather_result.init_inputs.build_inputs(worker.scratch, (Dim{:point}(1),))
        put!(weather_result.cache_pool, worker)
        if cache_weather
            mkpath(weather_cache_dir)
            serialize(disk_path, fetched)
        end
        _SILO_CACHE[wc_key] = fetched
    end
    return _SILO_CACHE[wc_key]
end

# Calendar-year-aligned span covering every requested year (MicroclimateMapper
# sizes its annual deep-soil climatology off a full Jan1-Dec31 span).
_year_span(years) = Date(minimum(years), 1, 1), Date(maximum(years), 12, 31)

# SILO-forced alternative to prepare_site: same site metadata, soil profile,
# canopy model and initial conditions (from the tower's own obs, so the two
# runs differ only in forcing), but environment_hourly/daily/minmax come from
# a live SILO fetch instead of the tower file.
#
# canopy_mode defaults to :full, matching prepare_site -- a standalone SILO
# run (gridded-forcing skill comparison) should get the same canopy-resolved
# output the tower-forced run produces. silo_gapfill_donor overrides this to
# :legacy for its own internal use: only reference_temperature/humidity/
# wind_speed/global_radiation/sky_temperature (profile-level, not
# canopy-layer output) feed run_site_gapfilled from there, so building a real
# canopy model would just slow that solve down for no benefit.
function prepare_site_silo(site_name, years; max_days=nothing, canopy_mode=:full)
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, latitude, longitude, height_series) = resolved
    sim_start, sim_end = _year_span(years)

    # Tower obs read first -- needed for observed depths (_build_depths) and
    # initial soil state before the SILO fetch/soil-texture fetch below.
    hourly = _aggregate_hourly(read_ozflux_nc(resolved.path, years; height_series))
    depths = _build_depths(SIM_DEPTHS_M, _observed_depths_m(hourly))

    # heights computed before the fetch (not just below, alongside canopy_model)
    # so fetch_silo_forcing can power-law-correct SILO's 10m wind to the SAME
    # reference height (last(heights)) the real model will assume it's at --
    # otherwise the fetch's own throwaway model silently defaults to 2m.
    extra_heights_m = [h for (h, _) in height_series]
    heights = _build_heights(canopy_height, reference_height, extra_heights_m)

    fetched = fetch_silo_forcing(site_name, latitude, longitude, sim_start, sim_end, depths, last(heights))

    ndays = length(sim_start:Day(1):sim_end)
    max_days !== nothing && (ndays = min(ndays, max_days))
    n = ndays * 24
    println("  $ndays days ($n hours) of SILO forcing")

    mean_air_temperature = mean(skipmissing(hourly.Ta)) * u"°C"
    initial_soil_temperature = initial_temperature_profile(_hourly_depth_series(hourly, "Ts"), depths, mean_air_temperature)

    site = Site(;
        latitude, longitude, elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), sky_view_fraction = 1.0,
        albedo = site_albedo(site_name), roughness_height = 0.004u"m",
        atmospheric_pressure = fetched.site.atmospheric_pressure,
    )
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths, resolved.global_attrib)
    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2; scale=sws_obs_scale(site_name, soil_profile))
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)  # bottom node always saturated (no drainage below the modelled soil column)

    # canopy_mode==:legacy here just means NoCanopy() -- the tower path's
    # fuller :legacy approximation (shade/wind-knockdown/horizon) isn't
    # wired up for SILO.
    canopy = canopy_model(site_name, resolved, heights, canopy_mode)
    model = micro_model(site_name, site, soil_profile, depths, heights, canopy, DailyRainfall())  # SILO rainfall is a daily total

    # max_days truncates the solve to `ndays` days, but fetch_silo_forcing
    # above is deliberately full-year -- truncate environment_minmax/
    # environment_hourly's own arrays to match (a no-op when max_days is
    # nothing).
    environment_minmax = max_days === nothing ? fetched.environment_minmax :
        truncate_environment_minmax(fetched.environment_minmax, ndays)
    environment_hourly = max_days === nothing ? fetched.environment_hourly :
        truncate_environment_hourly(fetched.environment_hourly, n)

    inputs = MicroInputs(;
        environment_minmax, environment_daily = fetched.environment_daily,
        environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    # Model spans the full sim_start:sim_end calendar range (SILO forcing);
    # `hourly` (tower obs, for initial conditions + comparison) may cover a
    # different/shorter span -- report.jl aligns the two by DateTime, not position.
    t_model = [DateTime(sim_start) + Day(d) + Hour(h) for d in 0:(ndays - 1) for h in 0:23]
    problem = MicroProblem(model, inputs; days = days_of_year(t_model), time_mode = ConsecutiveDayMode(; spinup_first_day = false))
    # sim_start/environment_daily kept (full sim_start:sim_end span, not
    # truncated by max_days) for silo_gapfill_donor's daily rainfall donor.
    return (; site_name, years, resolved, hourly, t_model, ndays, n, heights, depths, height_series, canopy_mode,
              sim_start, environment_daily = fetched.environment_daily, problem)
end

function run_site_silo(site_name, years; max_days=nothing, canopy_mode=:full)
    println("\n" * "="^72)
    println("Site $site_name  years=$years (SILO forcing, canopy_mode=$canopy_mode)")
    println("="^72)
    prep = prepare_site_silo(site_name, years; max_days, canopy_mode)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output)
end

# The SAME real+gap-filled hourly forcing prepare_site/run_site_gapfilled
# consumes directly (_gapfilled_hourly_forcing -- needs gap_fill_donor, or
# variables the tower doesn't measure at all, e.g. Whroo's Precip, come out
# all-zero rather than SILO-backed), reduced to a daily min/max/mean/total
# envelope; cloud_cover inferred from that same gap-filled Fsd's daily mean
# (clear-sky-ratio method, same as MicroclimateMapper's own SILO cloud
# inference -- MicroclimateMapper.jl/src/mesoclimate/cloud.jl -- and the
# same min=0.5x/max=2x spread its derive!(::CloudCover{Minimum/Maximum})
# uses), fed through Microclimate.jl's DailyMinMaxEnvironment/minmax_forcings
# -- the SAME diel-synthesis pathway a SILO run uses. Isolates whether
# MultilayerCanopy's response to raw hourly forcing *noise* (as opposed to
# real-vs-gridded *values*) is behind run_site_gapfilled's damped
# soil-temperature fluctuations relative to run_site_silo.
function _cloud_from_daily_fsd(daily_fsd_mean_Wm2, days_of_year, terrain)
    ndays = length(daily_fsd_mean_Wm2)
    hours = collect(0.0:1.0:23.0)
    solar_model = SolarRadiation.SolarProblem()
    solar_out = SolarRadiation.allocate_output_arrays(ndays * length(hours), ndays, solar_model.wavelength_count)
    solar_buffers = SolarRadiation.allocate_buffers(solar_model.wavelength_count, solar_model.diffuse_model)
    cloud = zeros(ndays)
    MicroclimateMapper.cloud_from_solar_radiation!(cloud, solar_out, daily_fsd_mean_Wm2 .* u"W/m^2",
        terrain, days_of_year, solar_buffers; solar_model, hours)
    return clamp.(cloud, 0.0, 1.0)
end

function prepare_site_daily(site_name, years; max_days=nothing, canopy_mode=:full, gap_fill_donor=nothing)
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, forcing_vars, latitude, longitude, leaf_area_index, height_series) = resolved

    println("  Reading $(resolved.path) for $site_name $years...")
    hourly = _aggregate_hourly(read_ozflux_nc(resolved.path, years; height_series))
    ndays_full = nrow(hourly) ÷ 24
    ndays = max_days === nothing ? ndays_full : min(ndays_full, max_days)
    hourly = hourly[1:(ndays * 24), :]
    n = ndays * 24
    println("  $ndays days ($n hours) after 30-min -> hourly aggregation, daily-envelope-forced")

    # Same real+gap-filled hourly forcing prepare_site/run_site_gapfilled
    # uses (SILO donor backs every variable's real gaps, incl. Precip --
    # e.g. Whroo's own tower Precip channel is 100% missing, so without this
    # the daily envelope would be built from real tower rain that doesn't
    # exist), degraded to a daily min/max/mean/total envelope instead of
    # consumed hourly.
    (; reference_temperature, reference_humidity, reference_wind_speed, global_radiation, rainfall_hourly) =
        _gapfilled_hourly_forcing(site_name, hourly, forcing_vars, gap_fill_donor, canopy_mode)
    frame = DataFrame(date = Date.(hourly.DateTime),
        ta = ustrip.(u"°C", reference_temperature), ws = ustrip.(u"m/s", reference_wind_speed),
        rh = reference_humidity .* 100.0, fsd = ustrip.(u"W/m^2", global_radiation),
        precip = ustrip.(u"kg/m^2", rainfall_hourly))
    daily = DataFrames.combine(DataFrames.groupby(frame, :date, sort=true),
        :ta => minimum => :ta_min, :ta => maximum => :ta_max,
        :ws => minimum => :ws_min, :ws => maximum => :ws_max,
        :rh => minimum => :rh_min, :rh => maximum => :rh_max,
        :fsd => mean => :fsd_mean, :precip => sum => :precip_total,
    )
    ndays_actual = nrow(daily)
    days_of_year = Dates.dayofyear.(daily.date)

    depths = _build_depths(SIM_DEPTHS_M, _observed_depths_m(hourly))
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths, resolved.global_attrib)

    solar_terrain = SolarRadiation.SolarTerrain(;
        elevation = 0.0u"m", horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
        albedo = site_albedo(site_name), atmospheric_pressure = atmospheric_pressure(0.0u"m"),
        latitude, longitude,
    )
    cloud_mean = _cloud_from_daily_fsd(daily.fsd_mean, days_of_year, solar_terrain)
    forcings = minmax_forcings(;
        reference_temperature_min = daily.ta_min .* u"°C", reference_temperature_max = daily.ta_max .* u"°C",
        reference_wind_speed_min = max.(daily.ws_min, 0.1) .* u"m/s", reference_wind_speed_max = max.(daily.ws_max, 0.1) .* u"m/s",
        reference_humidity_min = clamp.(daily.rh_min ./ 100.0, 0.0, 1.0), reference_humidity_max = clamp.(daily.rh_max ./ 100.0, 0.0, 1.0),
        cloud_cover_min = clamp.(cloud_mean .* 0.5, 0.0, 1.0), cloud_cover_max = clamp.(cloud_mean .* 2.0, 0.0, 1.0),
    )
    environment_minmax = DailyMinMaxEnvironment(; forcings)

    mean_air_temperature = mean((daily.ta_min .+ daily.ta_max) ./ 2) * u"°C"
    lp = legacy_params(site_name)
    shade_level = canopy_mode == :legacy ? _legacy_shade(leaf_area_index, lp.extinction_coefficient) : 0.0
    environment_daily = DailyTimeseries(;
        shade = fill(shade_level, ndays_actual),
        soil_wetness = fill(0.0, ndays_actual),
        surface_emissivity = fill(emissivity, ndays_actual),
        cloud_emissivity = fill(emissivity, ndays_actual),
        rainfall = daily.precip_total .* u"kg/m^2",
        deep_soil_temperature = fill(mean_air_temperature, ndays_actual),
        leaf_area_index = fill(leaf_area_index, ndays_actual),
    )

    horizon = canopy_mode == :legacy ? fill(lp.horizon_angle, 24) : fill(0.0u"°", 24)
    roughness_height = canopy_mode == :legacy ? lp.roughness_height : 0.004u"m"
    site = Site(;
        latitude, longitude, elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = horizon, sky_view_fraction = 1.0,
        albedo = site_albedo(site_name), roughness_height,
        atmospheric_pressure = atmospheric_pressure(0.0u"m"),
    )

    extra_heights_m = [h for (h, _) in height_series]
    heights = _build_heights(canopy_height, reference_height, extra_heights_m)

    canopy = canopy_model(site_name, resolved, heights, canopy_mode)
    model = micro_model(site_name, site, soil_profile, depths, heights, canopy, DailyRainfall())

    initial_soil_temperature = initial_temperature_profile(_hourly_depth_series(hourly, "Ts"), depths, mean_air_temperature)
    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2; scale=sws_obs_scale(site_name, soil_profile))
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)

    # A minmax-driven solve still reads environment_hourly.pressure directly
    # (interpolate_minmax!, Microclimate.jl/src/simulation.jl:386) -- nothing
    # else on it is used. Same stub shape MicroclimateMapper's own SILO fetch
    # uses (_environment_hourly_stub, weather.jl): every other field nothing.
    environment_hourly = HourlyTimeseries(;
        pressure = fill(atmospheric_pressure(0.0u"m"), ndays_actual * 24),
        reference_temperature = nothing, reference_humidity = nothing,
        reference_wind_speed = nothing, global_radiation = nothing,
        longwave_radiation = nothing, cloud_cover = nothing,
        rainfall = nothing, zenith_angle = nothing,
    )
    inputs = MicroInputs(;
        environment_minmax, environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    problem = MicroProblem(model, inputs; days = days_of_year, time_mode = ConsecutiveDayMode(; spinup_first_day = false))

    sim_start = daily.date[1]
    t_model = [DateTime(sim_start) + Day(d) + Hour(h) for d in 0:(ndays_actual - 1) for h in 0:23]

    return (; site_name, years, resolved, hourly, t_model, ndays = ndays_actual, n = ndays_actual * 24,
              heights, depths, height_series, canopy_mode, sim_start, environment_daily, problem)
end

function run_site_daily(site_name, years; max_days=nothing, canopy_mode=:full)
    gap_fill_donor = silo_gapfill_donor(site_name, years; max_days)
    println("\n" * "="^72)
    println("Site $site_name  years=$years (tower-daily-envelope forcing, canopy_mode=$canopy_mode)")
    println("="^72)
    prep = prepare_site_daily(site_name, years; max_days, canopy_mode, gap_fill_donor)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output)
end

# Real per-hour forcing donor for gap-filling, built from a SILO-forced solve.
# SILO's own `environment_hourly` is a stub (nothing fields; native resolution
# is daily min/max, not hourly) -- MicroResult's reference_*/global_radiation
# columns are the actual per-hour values the solver used (diel-synthesized
# from SILO's daily min/max + CRUCL2 wind climatology), recorded forcing
# echoes rather than simulated fluxes, so this is still "SILO data", not
# model predictions, feeding back into another run's forcing.
# Daily rainfall total isn't part of the per-hour donor above (SILO has no
# hourly rain at all, just a daily total -- see _gap_fill_precip_daily).
function silo_gapfill_donor(site_name, years; max_days=nothing)
    silo = run_site_silo(site_name, years; max_days, canopy_mode=:legacy)
    return (;
        t_model = silo.t_model,
        daily_dates = collect(silo.sim_start:Day(1):(silo.sim_start + Day(length(silo.environment_daily.rainfall) - 1))),
        daily_rainfall_mm = ustrip.(u"kg/m^2", silo.environment_daily.rainfall),
        reference_temperature_C = ustrip.(u"°C", silo.output.reference_temperature),
        reference_humidity_pct = silo.output.reference_humidity .* 100.0,
        reference_wind_speed_ms = ustrip.(u"m/s", silo.output.reference_wind_speed),
        global_radiation_Wm2 = ustrip.(u"W/m^2", silo.output.global_radiation),
        # SILO has no downward-longwave layer at all -- Fld is backed out of
        # the solver's own sky_temperature (Stefan-Boltzmann: incoming
        # longwave = σ*T_sky^4, the same relation view_factor.jl derives
        # sky_temperature FROM), itself computed from the SILO run's cloud
        # cover/humidity/temperature. A physically-grounded estimate, not a
        # straight-line interpolation across a multi-month tower gap.
        longwave_radiation_Wm2 = ustrip.(u"W/m^2", FluidProperties.σ .* silo.output.sky_temperature .^ 4),
    )
end

function run_site_gapfilled(site_name, years; max_days=nothing, canopy_mode=:full)
    gap_fill_donor = silo_gapfill_donor(site_name, years; max_days)
    println("\n" * "="^72)
    println("Site $site_name  years=$years (tower forcing, SILO-gap-filled)")
    println("="^72)
    prep = prepare_site(site_name, years; max_days, gap_fill_donor, canopy_mode)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output, solve_time)
end
