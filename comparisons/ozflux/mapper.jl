# mapper.jl — everything that needs MicroclimateMapper: SILO gridded forcing
# (run_site_silo) and live SLGA soil
# (soil_source :slga/:slga_uniform). Include after config.jl/utils.jl/
# pipeline.jl, with MicroclimateMapper, Rasters, RasterDataSources,
# PointDataSources, GeoInterface and Serialization loaded. The tower-forced
# run_site path doesn't need it.
#
# MicroclimateMapper still builds Microclimate's old MicroModel/MicroInputs
# (site and soil_profile on MicroInputs), so the fetches here need it ported
# before they run against the current Microclimate.

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", joinpath(homedir(), "Spatial_Data"))

# ── SILO weather source (secondary comparison) ──────────────────────────────
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
# MicroclimateMapper still carries its own pedotransfer types until it uses Microclimate.jl's
slga_pedotransfer_model_choice = MicroclimateMapper.CosbyMultivariate()
soil_area_buffer_deg = 0.05

# Live SLGA fetch + pedotransfer (same as oznet/pipeline.jl) for
# _fetch_soil_profile (pipeline.jl). Rescales to the site's own metadata
# BulkDensity attribute when present -- SLGA vs a real point measurement, not
# a substitute for it (SLGA is a 90 m grid mean).
function _slga_soil_profile(latitude, longitude, depths, global_attrib)
    b = soil_area_buffer_deg
    lon_dd, lat_dd = ustrip(u"°", longitude), ustrip(u"°", latitude)
    area = Extent(X = (lon_dd - b, lon_dd + b), Y = (lat_dd - b, lat_dd + b))
    soil_profile = build_soil_profile(SLGA, area; depths, pedotransfer_model = slga_pedotransfer_model_choice).soil_profile
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

# ── SILO gridded-forcing fetch (secondary comparison) ───────────────────────
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
# or the SILO-run wind ends up scaled to the wrong height (silently
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
                mineral_density = fill(2.56u"Mg/m^3", length(depths)), mineral_conductivity = fill(2.5u"W/m/K", length(depths)),
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
function prepare_site_silo(site_name, years; max_days=nothing)
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, latitude, longitude, height_series) = resolved
    sim_start, sim_end = _year_span(years)

    # Tower obs read first -- needed for observed depths (_build_depths) and
    # initial soil state before the SILO fetch/soil-texture fetch below.
    hourly = _aggregate_hourly(read_ozflux_nc(resolved.path, years, site_corrections(site_name); height_series))
    depths = _build_depths(SIM_DEPTHS_M, _observed_depths_m(hourly))

    # heights computed before the fetch (not just below, alongside canopy_model)
    # so fetch_silo_forcing can power-law-correct SILO's 10m wind to the SAME
    # reference height (last(heights)) the real model will assume it's at --
    # otherwise the fetch's own throwaway model silently defaults to 2m.
    heights = _build_heights(canopy_height, reference_height)

    fetched = fetch_silo_forcing(site_name, latitude, longitude, sim_start, sim_end, depths, last(heights))

    ndays = length(sim_start:Day(1):sim_end)
    max_days !== nothing && (ndays = min(ndays, max_days))
    n = ndays * 24
    println("  $ndays days ($n hours) of SILO forcing")

    (; reference_temperature_min, reference_temperature_max) = fetched.environment_minmax
    mean_air_temperature = mean((ustrip.(u"°C", reference_temperature_min) .+ ustrip.(u"°C", reference_temperature_max)) ./ 2) * u"°C"
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

    canopy = canopy_model(site_name, resolved, heights)
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
    return (; site_name, years, resolved, hourly, t_model, ndays, n, heights, depths, height_series,
              sim_start, environment_daily = fetched.environment_daily, problem)
end

function run_site_silo(site_name, years; max_days=nothing)
    println("\n" * "="^72)
    println("Site $site_name  years=$years (SILO forcing)")
    println("="^72)
    prep = prepare_site_silo(site_name, years; max_days)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output)
end

