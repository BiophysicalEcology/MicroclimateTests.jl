# mapper.jl — everything that needs MicroclimateMapper: SILO gridded forcing (run_site_silo)
# and SLGA soil (soil_source :slga/:slga_uniform). ozflux.jl includes it last.

# ── SILO weather source (secondary comparison) ──────────────────────────────
weather_source_choice = SILO
use_opendap_points = true  # PointDataSources.jl OPeNDAP point query, not a whole-raster download
dem_source_choice = SRTM
compute_terrain_choice = false

reuse_weather     = true
cache_weather     = true
weather_cache_dir = joinpath(@__DIR__, "weather_cache")

# ── SLGA soil at the site's cell, for _fetch_soil_profile (pipeline.jl) ─────
_slga_soil_profile(latitude, longitude, depths) =
    build_soil_profile(SLGA, (ustrip(u"°", longitude), ustrip(u"°", latitude));
        depths, pedotransfer_model = pedotransfer_model_choice).soil_profile

# ── SILO gridded-forcing fetch (secondary comparison) ───────────────────────
# A single-point MicroVectorProblem on the site's own model: MicroclimateMapper builds the
# point's MicroProblem, and its forcing (environment_minmax/daily/hourly) is kept. Cached on
# disk (weather_cache_dir) since a fetch is a live network call.
if use_opendap_points
    MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()
end

const _SILO_CACHE = Dict{Tuple{String,Date,Date,Float64},Any}()

# The forcing depends on the model only through its reference height: MicroclimateMapper
# power-law-corrects SILO's 10 m wind to maximum(model.atmosphere.heights)
function fetch_silo_forcing(site_name, model, latitude, longitude, sim_start::Date, sim_end::Date, initial_soil_moisture)
    ref_m = ustrip(u"m", maximum(model.atmosphere.heights))
    wc_key = (site_name, sim_start, sim_end, ref_m)
    disk_path = joinpath(weather_cache_dir, "SILO_$(site_name)_$(sim_start)_$(sim_end)_ref$(ref_m)m.jls")
    if !haskey(_SILO_CACHE, wc_key) && reuse_weather && isfile(disk_path)
        _SILO_CACHE[wc_key] = deserialize(disk_path)
        println("  Loaded SILO forcing from disk cache: $disk_path")
    end
    if !haskey(_SILO_CACHE, wc_key)
        println("  Fetching SILO forcing ($sim_start to $sim_end, wind reference height $(ref_m) m)...")
        mapper_model = MicroMapModel(; micro_model = model, dem_source = dem_source_choice,
            weather_source = weather_source_choice, compute_terrain = compute_terrain_choice,
            surface_albedo_source = site_albedo(site_name), roughness_height_source = site_ground_roughness_height(site_name))
        point = GIW.Point((ustrip(u"°", longitude), ustrip(u"°", latitude)))
        vec_problem = MicroVectorProblem(; model = mapper_model, points = [point], dates = sim_start:Day(1):sim_end,
            init = (; soil_moisture = initial_soil_moisture))
        cache = MicroclimateMapper.init(vec_problem)
        worker = take!(cache.cache_pool)
        problem = cache.init_inputs.build_problem(worker.scratch, (Dim{:point}(1),))
        # The inputs alias the worker's scratch buffers: copy before handing the worker back
        (; environment_minmax, environment_daily, environment_hourly) = problem.inputs
        fetched = deepcopy((; environment_minmax, environment_daily, environment_hourly))
        put!(cache.cache_pool, worker)
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

# SILO-forced alternative to prepare_site: the same site, model and initial conditions (from
# the tower's own obs, so the two runs differ only in forcing), with environment_minmax/daily/
# hourly from a live SILO fetch instead of the tower file.
function prepare_site_silo(site_name, years; max_days=nothing, litter_model=site_litter(site_name))
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, latitude, longitude, height_series) = resolved
    sim_start, sim_end = _year_span(years)

    hourly = _aggregate_hourly(read_ozflux_nc(resolved.path, years, site_corrections(site_name); height_series))
    depths = _build_depths(SIM_DEPTHS_M, vcat(_observed_depths_m(hourly), extra_soil_depths_m))
    heights = _build_heights(canopy_height, reference_height, first.(height_series))

    site = Site(;
        latitude, longitude, resolved.elevation, slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), sky_view_fraction = 1.0,
    )
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths)
    canopy = canopy_model(site_name, resolved, heights)
    model = micro_model(site_name, site, soil_profile, depths, heights, canopy, DailyRainfall(); litter_model)  # SILO rainfall is a daily total

    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2, soil_profile; scale=sws_obs_scale(site_name, soil_profile))
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)  # bottom node always saturated (no drainage below the modelled soil column)

    fetched = fetch_silo_forcing(site_name, model, latitude, longitude, sim_start, sim_end, initial_soil_moisture)

    ndays = length(sim_start:Day(1):sim_end)
    max_days !== nothing && (ndays = min(ndays, max_days))
    n = ndays * 24
    println("  $ndays days ($n hours) of SILO forcing")

    (; reference_temperature_min, reference_temperature_max) = fetched.environment_minmax
    mean_air_temperature = mean((ustrip.(u"°C", reference_temperature_min) .+ ustrip.(u"°C", reference_temperature_max)) ./ 2) * u"°C"
    initial_soil_temperature = initial_temperature_profile(_hourly_depth_series(hourly, "Ts"), depths, mean_air_temperature)

    # max_days truncates the solve to `ndays` days, but fetch_silo_forcing
    # above is deliberately full-year -- truncate the forcing to match (a no-op
    # when max_days is nothing).
    environment_minmax = max_days === nothing ? fetched.environment_minmax :
        truncate_environment_minmax(fetched.environment_minmax, ndays)
    environment_hourly = max_days === nothing ? fetched.environment_hourly :
        truncate_environment_hourly(fetched.environment_hourly, n)

    inputs = MicroInputs(;
        environment_minmax, fetched.environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    # Model spans the full sim_start:sim_end calendar range (SILO forcing);
    # `hourly` (tower obs, for initial conditions + comparison) may cover a
    # different/shorter span -- report.jl aligns the two by DateTime, not position.
    t_model = [DateTime(sim_start) + Day(d) + Hour(h) for d in 0:(ndays - 1) for h in 0:23]
    problem = MicroProblem(model, inputs; days = days_of_year(t_model), time_mode = ConsecutiveDayMode(; spinup_first_day = false))
    return (; site_name, years, resolved, hourly, t_model, ndays, n, heights, depths, height_series,
              sim_start, fetched.environment_daily, problem)
end

function run_site_silo(site_name, years; max_days=nothing)
    haskey(ENV, "SILO_EMAIL") || error("SILO forcing needs your email in ENV[\"SILO_EMAIL\"], as SILO's API asks")
    println("\n" * "="^72)
    println("Site $site_name  years=$years (SILO forcing)")
    println("="^72)
    prep = prepare_site_silo(site_name, years; max_days)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output)
end

