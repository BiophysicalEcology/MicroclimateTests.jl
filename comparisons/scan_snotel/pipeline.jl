# pipeline.jl — shared model setup + per-site pipeline for comparison.jl
# (multi-site, phase-batched) and debug_site.jl (single-site, flat, no
# try/catch — for stepping through in the REPL/debugger).
#
# Pipeline phases, each usable standalone or via `run_site` (sequential,
# single-site composition of the same pieces):
#   prefetch_weather_batch!  — one MicroVectorProblem covering many sites
#   prepare_site             — metadata, weather (from cache), obs, initial
#                              conditions, builds the final MicroProblem
#   run_julia_batch!         — Threads.@threads over many prepared sites
#   run_nmr_batch!            — concurrent Rscript processes, capped
#   report_site_results      — NMR read, stats, save, plots for one site
#
# Single source of truth: change the model in build_micro_model() once, both
# scripts pick it up. Requires config.jl and utils.jl to already be included.

# ── Site-independent model setup ─────────────────────────────────────────────
# Built once per script run (not per-site — nothing here depends on site_num),
# using the MODEL ALGORITHM CHOICES section of config.jl. That's the one
# place to change e.g. `infiltration_algorithm_choice` or `dem_source_choice`.
function build_micro_model()
    soil_properties_model = CampbelldeVriesSoilProperties(;
        de_vries_shape_factor,
        recirculation_power,
        return_flow_threshold,
    )

    soil_profile = SoilProfile(;
        bulk_density    = fill(bulk_density, length(depths)),
        mineral_density = fill(mineral_density, length(depths)),
        mineral_conductivity,
        mineral_heat_capacity,
        hydraulics = CampbellHydraulicProfile(;
            air_entry_water_potential        = fill(air_entry_potential, length(depths)),
            saturated_hydraulic_conductivity = fill(sat_hydraulic_cond, length(depths)),
            campbell_b_parameter             = fill(campbell_b, length(depths)),
            root_density,
        ),
    )
    soil_hydraulic_model = CampbellSoilHydraulics(;
        root_resistance,
        stomatal_closure_potential,
        leaf_resistance,
        stomatal_stability_parameter,
        root_radius,
        infiltration_algorithm = infiltration_algorithm_choice,
    )

    snow_model = SnowModel(;
        snow_temperature_threshold = snow_temp_threshold,
        snow_density               = snow_density_init,
        snow_melt_factor           = snow_melt_factor,
        undercatch                 = undercatch,
        rain_multiplier            = rain_multiplier,
        rain_melt_factor           = rain_melt_factor,
        density_function           = density_function,
        snow_conductivity          = snow_conductivity,
        canopy_interception        = canopy_interception,
    )

    micro_config = MicroConfig(;
        convergence            = convergence_choice,
        rainfall_schedule      = rainfall_schedule_choice,
        soil_moisture_strategy = soil_moisture_strategy_choice,
    )

    # `radiation`/`boundary_layer_model` are left at MicroModel's defaults
    # (RadiationModel() / MoninObukhov()) — they already match the values
    # this script used to set explicitly, so no `using SolarRadiation` needed.
    micro_model = MicroModel(;
        hours   = collect(0.0:1:23.0),
        depths,
        heights = profile_heights,
        soil_properties_model,
        soil_hydraulic_model,
        snow_model,
        config = micro_config,
    )

    # NOTE: `compute_terrain_choice = false` keeps the flat-terrain assumption
    # (zero slope/aspect, flat horizon), matching NicheMapR's point-run
    # convention. Site elevation comes from `dem_source_choice` (SRTM) at the
    # site's lon/lat rather than the surveyed station elevation in the
    # metadata CSV; for GRIDMET, gridMET's own elevation grid is additionally
    # used for lapse-rate correction via MicroclimateMapper's
    # `weather_grid_elevation(GRIDMET)` — ERA5/NCEP do their own equivalent
    # internally.
    mapper_model = MicroMapModel(;
        micro_model,
        dem_source              = dem_source_choice,
        weather_source           = weather_source_choice,
        surface_albedo_source    = albedo,
        roughness_height_source  = roughness_height,
        compute_terrain          = compute_terrain_choice,
    )

    return micro_model, soil_profile, mapper_model
end

# ── Site metadata + date-range resolution ────────────────────────────────────
# Shared by every phase below so they all agree on exactly which sites are
# valid and what date range each one gets. Returns `nothing` for expected skip
# conditions (not in metadata, outside weather-source coverage, no usable
# obs/weather date overlap) — these are routine, not errors.
function resolve_site(site_num, meta_all; sim_start, sim_end, auto_date_range, max_sim_years)
    obsfile = joinpath(@__DIR__, "observations", "$(site_num).csv")

    site_row = filter(r -> strip(string(r.ID)) == string(site_num), meta_all)
    if nrow(site_row) == 0
        @warn "Site $site_num not found in metadata CSV — lat/lon/elevation unknown."
        lat_dd = NaN;  lon_dd = NaN;  elev_m = NaN
        site_name = string(site_num);  site_state = "";  site_network = ""
    else
        r           = site_row[1, :]
        lat_dd      = parse(Float64, string(r.Latitude))
        lon_dd      = parse(Float64, string(r.Longitude))   # negative = W
        elev_m      = parse(Float64, string(r.Elevation_ft)) / 3.28084
        site_name   = string(r.Name)
        site_state  = string(r.State)
        site_network= string(r.Network)
    end
    println("  $site_name ($site_state, $site_network) — " *
            "$(round(lat_dd, digits=4))°N  $(round(lon_dd, digits=4))°E  " *
            "$(round(elev_m, digits=1)) m asl")

    # ── Weather-source coverage check ─────────────────────────────────────────
    # Only GRIDMET (CONUS ~4 km) needs a coverage check — ERA5/NCEP are global.
    if weather_source_choice == GRIDMET
        if site_state in GRIDMET_EXCLUDED_STATES
            println("  Skipping — state $site_state is outside gridMET coverage.")
            return nothing
        end
        if isnan(lat_dd) || isnan(lon_dd) ||
           lat_dd < GRIDMET_LAT_MIN || lat_dd > GRIDMET_LAT_MAX ||
           lon_dd < GRIDMET_LON_MIN || lon_dd > GRIDMET_LON_MAX
            println("  Skipping — coordinates ($(round(lat_dd,digits=3))°N, $(round(lon_dd,digits=3))°E) outside gridMET domain.")
            return nothing
        end
        grid_elevation = GRIDMET_ELEV_RASTER[X(Near(lon_dd)), Y(Near(lat_dd))] * u"m"
        println("  GridMET grid elevation: $(round(ustrip(grid_elevation), digits=1)) m")
    elseif isnan(lat_dd) || isnan(lon_dd)
        println("  Skipping — coordinates unknown.")
        return nothing
    end

    # ── Auto date range ───────────────────────────────────────────────────────
    _sim_start = sim_start
    _sim_end   = sim_end
    weather_start = WEATHER_SOURCE_START[weather_source_choice]
    if auto_date_range && isfile(obsfile)
        _obs_dates_raw = DataFrame(CSV.File(obsfile; select=["DateTime"],
            normalizenames=true, silencewarnings=true))
        _obs_dates = DateTime.(_obs_dates_raw[!, :DateTime],
                               dateformat"yyyy-mm-dd HH:MM:SS")
        _obs_data  = DataFrame(CSV.File(obsfile; normalizenames=true,
            missingstring=["NA",""], silencewarnings=true))
        _snow_cols = [:SNWD_I, :WTEQ_I]   # prefer anchoring to snow obs availability
        _soil_cols = [:STO_I_2, :STO_I_4, :STO_I_8, :STO_I_20, :STO_I_40,
                      :SMS_I_2, :SMS_I_4, :SMS_I_8, :SMS_I_20, :SMS_I_40]
        _has_data_for(cols) = begin
            mask = falses(nrow(_obs_data))
            for c in cols
                c in propertynames(_obs_data) &&
                    (mask .|= .!ismissing.(_obs_data[!, c]))
            end
            mask
        end
        _has_data = _has_data_for(_snow_cols)
        _valid_dt = _obs_dates[_has_data]
        if length(_valid_dt) < 2
            # No usable snow record at this site (e.g. a non-snow SCAN site) —
            # fall back to soil temp/moisture obs so it still gets a real,
            # data-driven multi-year window and contributes soil-only stats,
            # rather than silently defaulting to whatever sim_start/sim_end
            # happens to be hardcoded in the entry-point script.
            _has_data = _has_data_for(_soil_cols)
            _valid_dt = _obs_dates[_has_data]
        end
        if length(_valid_dt) >= 2
            _sim_end   = min(Date(maximum(_valid_dt)), today() - Day(1))
            _sim_start = max(_sim_end - Year(max_sim_years) + Day(1),
                             Date(minimum(_valid_dt)), weather_start)
        end
        if _sim_start >= _sim_end
            @warn "Site $site_num: no usable obs/weather overlap — skipping."
            return nothing
        end
        println("  Auto date range: $_sim_start to $_sim_end")
    end

    # Snap to full calendar-year boundaries (Jan 1 .. Dec 31). MicroclimateMapper.jl
    # sizes its internal weather buffers — and computes annual-mean deep-soil
    # climatology — off the full Jan1-Dec31 span of the touched year(s), but
    # slices the loaded weather stack down to just the requested sub-range
    # before assembling those buffers. Any request that isn't itself full-year
    # aligned leaves the tail of those buffers unfilled with garbage, corrupting
    # downstream derivations (surfaces as a DomainError deep in vapour_pressure()).
    # Requesting full years sidesteps the mismatch entirely, at the cost of
    # simulating a few extra months of data at each end. Clamped so the end
    # doesn't run into the weather source's actual availability (today - 1
    # day) — a request ending in the current (incomplete) calendar year still
    # can't be fully year-aligned; that residual case isn't covered by this
    # workaround.
    _sim_start = max(Date(Dates.year(_sim_start), 1, 1), weather_start)
    _sim_end   = min(Date(Dates.year(_sim_end), 12, 31), today() - Day(1))

    return (; lat_dd, lon_dd, elev_m, site_name, site_state, site_network,
              sim_start = _sim_start, sim_end = _sim_end)
end

# Disk path for one site's cached weather forcing (see config.jl's
# reuse_weather/cache_weather/weather_cache_dir). Namespaced by weather source
# so switching weather_source_choice can't silently reuse another source's cache.
_weather_cache_path(site_num, sim_start, sim_end) =
    joinpath(weather_cache_dir, "$(nameof(weather_source_choice))_$(site_num)_$(sim_start)_$(sim_end).jls")

# ── Batch-prefetch weather forcing for many sites in one MicroVectorProblem ──
# Reading each (variable, year) file's cropped region separately per site
# means redundantly re-decompressing overlapping/nearby chunks once per site.
# Passing all sites' points to a single MicroVectorProblem call reads each
# region once, then extracts every site's own series from the shared result.
# Only makes sense when sites share (or nearly share) a date range — this
# groups sites by their resolved (sim_start, sim_end) and does one batched
# call per group, so it degrades gracefully even if ranges do vary.
# Populates `weather_cache` in place; `prepare_site`/`run_site` already check
# the cache first, so no other changes are needed for them to pick up the
# pre-fetched data.
function prefetch_weather_batch!(weather_cache, meta_all, sites, mapper_model, soil_profile;
                                  sim_start, sim_end, auto_date_range, max_sim_years)
    groups = Dict{Tuple{Date,Date}, Vector{Tuple{Int,Float64,Float64}}}()
    for site_num in sites
        println("\nResolving site $site_num for batch fetch...")
        resolved = resolve_site(site_num, meta_all; sim_start, sim_end, auto_date_range, max_sim_years)
        isnothing(resolved) && continue
        wc_key = (site_num, resolved.sim_start, resolved.sim_end)
        haskey(weather_cache, wc_key) && continue   # already fetched (e.g. from a previous run)
        _cache_file = _weather_cache_path(site_num, resolved.sim_start, resolved.sim_end)
        if reuse_weather && isfile(_cache_file)
            weather_cache[wc_key] = deserialize(_cache_file)
            println("  Loaded from disk cache: $_cache_file")
            continue
        end
        key = (resolved.sim_start, resolved.sim_end)
        push!(get!(groups, key, Tuple{Int,Float64,Float64}[]), (site_num, resolved.lon_dd, resolved.lat_dd))
    end

    println("\n$(length(groups)) distinct date range(s) across $(sum(length(v) for v in values(groups); init=0)) site(s) to fetch.")

    for ((_sim_start, _sim_end), entries) in groups
        # Sub-partition into spatial tiles: the SRTM DEM load covers the
        # bounding box of every point handed to one MicroVectorProblem at
        # once, so grouping geographically scattered sites (e.g. a
        # nationwide site_subset) into one call blows the DEM read up to a
        # multi-GB, whole-region raster. Binning by a fixed-size lon/lat grid
        # cell keeps each batch's bounding box close to what a single site's
        # own buffer already needed — only sites that are genuinely close
        # together get batched.
        tiles = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Float64,Float64}}}()
        for e in entries
            (_, lon_dd, lat_dd) = e
            tile_key = (floor(Int, lon_dd / PREFETCH_TILE_DEG), floor(Int, lat_dd / PREFETCH_TILE_DEG))
            push!(get!(tiles, tile_key, Tuple{Int,Float64,Float64}[]), e)
        end
        # Further cap tile size (a single dense tile could still exceed the
        # secondary per-tile site cap, e.g. many stations in one small area).
        batches = Vector{Tuple{Int,Float64,Float64}}[]
        for tile_entries in values(tiles)
            for i in 1:PREFETCH_MAX_SITES_PER_TILE:length(tile_entries)
                push!(batches, tile_entries[i:min(i+PREFETCH_MAX_SITES_PER_TILE-1, length(tile_entries))])
            end
        end
        println("\n$(length(entries)) site(s) for $_sim_start to $_sim_end split into $(length(batches)) spatial batch(es) (tile=$(PREFETCH_TILE_DEG)°).")

        for batch in batches
            println("\nBatch-fetching weather forcing for $(length(batch)) site(s) ($_sim_start to $_sim_end)...")
            points = [GIW.Point((lon_dd, lat_dd)) for (_, lon_dd, lat_dd) in batch]
            dates_range = _sim_start:Day(1):_sim_end
            vec_problem = MicroVectorProblem(;
                model = mapper_model,
                points,
                dates  = dates_range,
                soil_profile,
                init = (; soil_moisture = fill(0.3, length(depths))),
            )
            init_time = @elapsed weather_result = MicroclimateMapper.init(vec_problem)
            @printf("  MicroclimateMapper.init (DEM + weather fetch, %d points): %.2f s\n", length(batch), init_time)
            worker = take!(weather_result.cache_pool)
            for (i, (site_num, _, _)) in enumerate(batch)
                inputs = weather_result.init_inputs.build_inputs(worker.scratch, (Dim{:point}(i),))
                weather_cache[(site_num, _sim_start, _sim_end)] = inputs
                if cache_weather
                    mkpath(weather_cache_dir)
                    serialize(_weather_cache_path(site_num, _sim_start, _sim_end), inputs)
                end
            end
            put!(weather_result.cache_pool, worker)
        end
    end
    return weather_cache
end

# ── Prepare one site: metadata, weather, obs, initial conditions, problem ───
# Does everything up through building the final `MicroProblem`, but doesn't
# run anything (no NMR launch, no Julia solve) — that's `run_nmr_batch!`/
# `run_julia_batch!`, or `run_site`'s synchronous single-site composition.
# Returns `nothing` for the same routine skip conditions as `resolve_site`.
function prepare_site(site_num, meta_all, weather_cache, micro_model, soil_profile;
                       sim_start, sim_end, auto_date_range, max_sim_years)
    obsfile = joinpath(@__DIR__, "observations", "$(site_num).csv")

    println("\n" * "="^72)
    println("Site $site_num")
    println("="^72)

    println("\nReading site metadata...")
    resolved = resolve_site(site_num, meta_all; sim_start, sim_end, auto_date_range, max_sim_years)
    isnothing(resolved) && return nothing
    (; lat_dd, lon_dd, elev_m) = resolved
    _sim_start = resolved.sim_start
    _sim_end   = resolved.sim_end

    # ── Daily forcing from the chosen weather source (via MicroclimateMapper) ─
    dates_range = _sim_start:Day(1):_sim_end
    wc_key = (site_num, _sim_start, _sim_end)
    _weather_disk_cache = _weather_cache_path(site_num, _sim_start, _sim_end)
    if !haskey(weather_cache, wc_key) && reuse_weather && isfile(_weather_disk_cache)
        weather_cache[wc_key] = deserialize(_weather_disk_cache)
        println("\nLoaded weather forcing from disk cache: $_weather_disk_cache")
    end
    if !haskey(weather_cache, wc_key)
        # Not pre-fetched via prefetch_weather_batch! — fetch just this one site.
        println("\nFetching weather forcing ($_sim_start to $_sim_end)...")
        point = GIW.Point((lon_dd, lat_dd))
        vec_problem = MicroVectorProblem(;
            model = MicroMapModel(; micro_model,
                dem_source = dem_source_choice, weather_source = weather_source_choice,
                surface_albedo_source = albedo, roughness_height_source = roughness_height,
                compute_terrain = compute_terrain_choice),
            points = [point],
            dates  = dates_range,
            soil_profile,
            # Most weather sources don't provide soil_moisture natively, so
            # `init()` needs a placeholder here to build `base_inputs` — this
            # value is discarded; the real initial soil moisture (obs-derived
            # `initial_sm`) is set on the final `MicroInputs` below.
            init = (; soil_moisture = fill(0.3, length(depths))),
        )
        init_time = @elapsed weather_result = MicroclimateMapper.init(vec_problem)
        @printf("  MicroclimateMapper.init (DEM + weather fetch): %.2f s\n", init_time)
        worker = take!(weather_result.cache_pool)
        inputs = weather_result.init_inputs.build_inputs(worker.scratch, (Dim{:point}(1),))
        put!(weather_result.cache_pool, worker)
        weather_cache[wc_key] = inputs
        if cache_weather
            mkpath(weather_cache_dir)
            serialize(_weather_disk_cache, inputs)
        end
    else
        println("\nUsing cached weather forcing ($_sim_start to $_sim_end).")
    end
    base_inputs = weather_cache[wc_key]

    dates_vec = MicroclimateMapper._normalise_dates(dates_range)  # Vector{Date}, real calendar (leap years included)
    f_doy    = dayofyear.(dates_vec)   # real calendar DOY
    ndays    = length(f_doy)
    f_tminn  = base_inputs.environment_minmax.reference_temperature_min  # Unitful K
    f_tmaxx  = base_inputs.environment_minmax.reference_temperature_max  # Unitful K
    f_rhminn = base_inputs.environment_minmax.reference_humidity_min     # 0–1 fraction
    f_rhmaxx = base_inputs.environment_minmax.reference_humidity_max     # 0–1 fraction
    f_ccmaxx = base_inputs.environment_minmax.cloud_max                 # 0–1 fraction
    f_ccminn = base_inputs.environment_minmax.cloud_min                 # 0–1 fraction
    f_wnminn = base_inputs.environment_minmax.reference_wind_min        # Unitful m/s
    f_wnmaxx = base_inputs.environment_minmax.reference_wind_max        # Unitful m/s
    f_rain   = base_inputs.environment_daily.rainfall                  # Unitful kg/m²
    f_tannul = base_inputs.environment_daily.deep_soil_temperature      # Unitful K
    println("  $ndays days from $(nameof(weather_source_choice))")

    # ── Observations ───────────────────────────────────────────────────────
    println("\nReading observations from $obsfile...")
    obs_raw = DataFrame(CSV.File(obsfile, normalizenames=true,
        missingstring=["NA",""], silencewarnings=true))
    obs_raw[!, :DateTime] = DateTime.(string.(obs_raw[!, :DateTime]),
        dateformat"yyyy-mm-dd HH:MM:SS")
    omask   = (obs_raw.DateTime .>= DateTime(_sim_start)) .&
              (obs_raw.DateTime .<= DateTime(_sim_end, Time(23,0,0)))
    obs_sub = obs_raw[omask, :]

    obs = select(obs_sub, :DateTime)
    # Column mapping: obs name (after normalizenames) → output name, unit conversion
    col_map = [
        (:SNWD_I,   :SNWD_cm,   x -> ismissing(x) ? missing : Float64(x) * 2.54),
        (:WTEQ_I,   :WTEQ_cm,   x -> ismissing(x) ? missing : Float64(x) * 2.54),
        (:STO_I_2,  :STO_5cm,   x -> x),
        (:STO_I_4,  :STO_10cm,  x -> x),
        (:STO_I_8,  :STO_20cm,  x -> x),
        (:STO_I_20, :STO_50cm,  x -> x),
        (:STO_I_40, :STO_100cm, x -> x),
        (:SMS_I_2,  :SMS_5cm,   x -> x),
        (:SMS_I_4,  :SMS_10cm,  x -> x),
        (:SMS_I_8,  :SMS_20cm,  x -> x),
        (:SMS_I_20, :SMS_50cm,  x -> x),
        (:SMS_I_40, :SMS_100cm, x -> x),
    ]
    for (src, dst, tfm) in col_map
        src in propertynames(obs_sub) && (obs[!, dst] = tfm.(obs_sub[!, src]))
    end
    println("  $(nrow(obs)) obs rows; columns: $(join(names(obs), ", "))")
    println("  De-spiking observations...")
    despiked!(obs, obs_despike_specs)

    # ── Initial conditions ─────────────────────────────────────────────────
    # Node depths (cm): 0,1.25,2.5,3.75,5,7.5,10,12.5,15,17.5,20,25,30,40,50,75,100,150,200
    # Obs available at 5→idx5, 10→idx7, 20→idx11, 50→idx15, 100→idx17
    node_depths_cm = ustrip.(depths) .* 100.0   # depths is in m, convert to cm

    t_air_mean = f_tannul[1]  # Unitful K — mean annual temperature
    t_air_day1 = (f_tminn[1] + f_tmaxx[1]) / 2  # daily mean for first day

    # Find first obs row on or after sim_start
    first_obs_rows = obs[obs.DateTime .>= DateTime(_sim_start), :]
    first_obs = nrow(first_obs_rows) > 0 ? first_obs_rows[1, :] : nothing

    # --- Soil temperature ---
    # Surface anchor: use daily mean air temp, but clamp up to the shallowest
    # available obs if the soil is warmer (snow insulation keeps soil near 0°C
    # even when air is well below freezing — the extreme gradient otherwise
    # destabilises the ODE solver at hour 1).
    st_obs_map = [(:STO_5cm, 5), (:STO_10cm, 7), (:STO_20cm, 11),
                  (:STO_50cm, 15), (:STO_100cm, 17)]
    surface_anchor = t_air_day1
    if !isnothing(first_obs)
        for (col, _) in st_obs_map
            if col in propertynames(first_obs) && !ismissing(first_obs[col])
                obs_T = (Float64(first_obs[col]) + 273.15) * u"K"
                surface_anchor = max(surface_anchor, obs_T)
                break  # use only the shallowest available obs
            end
        end
    end
    st_known_d = Float64[0.0]
    st_known_T = typeof(surface_anchor)[surface_anchor]
    if !isnothing(first_obs)
        for (col, idx) in st_obs_map
            if col in propertynames(first_obs) && !ismissing(first_obs[col])
                push!(st_known_d, node_depths_cm[idx])
                push!(st_known_T, (Float64(first_obs[col]) + 273.15) * u"K")
            end
        end
    end
    push!(st_known_d, 200.0)   # deep anchor
    push!(st_known_T, t_air_mean)

    initial_st = [_linterp(st_known_d, st_known_T, d) for d in node_depths_cm]
    println("  Initial ST (node 1→19): " *
        join([@sprintf("%.1f", ustrip(uconvert(u"°C", T))) for T in initial_st], ", ") * " °C")

    # --- Soil moisture ---
    sat_sm = ustrip(saturation_moisture)
    sm_obs_map = [(:SMS_5cm, 5), (:SMS_10cm, 7), (:SMS_20cm, 11),
                  (:SMS_50cm, 15), (:SMS_100cm, 17)]
    sm_known_d = Float64[];  sm_known_v = Float64[]
    if !isnothing(first_obs)
        for (col, idx) in sm_obs_map
            if col in propertynames(first_obs) && !ismissing(first_obs[col])
                v = Float64(first_obs[col])
                v > 1.0 && (v /= 100.0)   # guard: convert % → fraction if needed
                push!(sm_known_d, node_depths_cm[idx])
                push!(sm_known_v, clamp(v, 0.01, sat_sm))
            end
        end
    end
    push!(sm_known_d, 200.0)   # deep anchor: saturated at 200 cm
    push!(sm_known_v, sat_sm)
    # Clamp lower bound to 0.01 — SMS sensors read 0 when frozen, but Campbell
    # matric potential ψ = ψe*(θ/θs)^(-b) → ±Inf at θ=0, causing NaN in ODE RHS.
    initial_sm = if length(sm_known_d) >= 2
        clamp.([_linterp(sm_known_d, sm_known_v, d) for d in node_depths_cm], 0.01, sat_sm)
    else
        fill(sat_sm, 19)
    end
    println("  Initial SM (node 1→19): " *
        join([@sprintf("%.3f", v) for v in initial_sm], ", "))

    # --- Snow depth ---
    initial_snow_depth = if init_snow_from_obs && !isnothing(first_obs) &&
                            :SNWD_cm in propertynames(first_obs) &&
                            !ismissing(first_obs[:SNWD_cm])
        Float64(first_obs[:SNWD_cm]) * u"cm"
    else
        0.0u"cm"
    end
    println("  Initial snow depth: $(round(ustrip(initial_snow_depth), digits=1)) cm")

    # ── Build the final MicroProblem ──────────────────────────────────────────
    # `micro_model`/`soil_profile` were already built (site-independent);
    # `base_inputs.site`/`.environment_minmax`/`.environment_daily`/
    # `.environment_hourly` came from the MicroclimateMapper weather fetch —
    # only the obs-derived initial conditions are added here.
    inputs = MicroInputs(;
        site               = base_inputs.site,
        soil_profile,
        environment_minmax = base_inputs.environment_minmax,
        environment_daily  = base_inputs.environment_daily,
        environment_hourly = base_inputs.environment_hourly,
        initial_soil_temperature = initial_st,
        initial_soil_moisture    = Vector{Float64}(initial_sm),
        initial_snow_depth       = initial_snow_depth,
        initial_snow_temperature = u"K"(0.0u"°C"),
    )

    problem = MicroProblem(micro_model, inputs;
        days      = f_doy,
        time_mode = ConsecutiveDayMode(; spinup_first_day=false),
    )

    return (; site_num, lat_dd, lon_dd, elev_m,
              sim_start = _sim_start, sim_end = _sim_end,
              dates_vec, f_doy, ndays,
              f_tminn, f_tmaxx, f_rhminn, f_rhmaxx, f_ccmaxx, f_ccminn,
              f_wnminn, f_wnmaxx, f_rain, f_tannul,
              obs, initial_st, initial_sm, initial_snow_depth,
              problem)
end

# ── Write NMR input CSVs for one prepared site ───────────────────────────────
# Separated from launching Rscript so `run_nmr_batch!` can write all sites'
# inputs up front, then launch/wait on the R processes concurrently.
# Returns (nmr_site_dir, needs_run::Bool).
function write_nmr_inputs(prep)
    (; site_num, lat_dd, lon_dd, elev_m, sim_start, sim_end, f_doy,
       f_tminn, f_tmaxx, f_rhminn, f_rhmaxx, f_wnminn, f_rain, f_ccmaxx, f_ccminn,
       f_tannul, initial_st, initial_sm) = prep
    _sim_start, _sim_end = sim_start, sim_end

    _nmr_site_dir = joinpath(nmr_out_dir, string(site_num))
    _nmr_metout   = joinpath(_nmr_site_dir, "metout.csv")
    needs_run = run_nmr && (!reuse_nmr || !isfile(_nmr_metout))
    if needs_run
        mkpath(_nmr_site_dir)
        # NMR DEP nodes are subset of Julia's 19 at indices [1,3,5,7,9,11,13,15,17,19]
        _nmr_idx  = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
        _mc  = ustrip.(mineral_conductivity)   # W/m/K
        _mhc = ustrip.(mineral_heat_capacity)  # J/kg/K
        _rd  = ustrip.(root_density)           # m/m³
        # nmr_forcing.csv
        CSV.write(joinpath(_nmr_site_dir, "nmr_forcing.csv"), DataFrame(
            doy       = f_doy,   # real calendar DOY (leap years included) — matches f_tminn etc. length
            Tmin_C    = ustrip.(uconvert.(u"°C", f_tminn)),
            Tmax_C    = ustrip.(uconvert.(u"°C", f_tmaxx)),
            RHmin_pct = f_rhminn .* 100.0,
            RHmax_pct = f_rhmaxx .* 100.0,
            Wind_ms   = ustrip.(f_wnminn),
            Rain_mm   = ustrip.(f_rain),
            CCmax_pct = f_ccmaxx .* 100.0,
            CCmin_pct = f_ccminn .* 100.0,
            tannul_C  = ustrip.(uconvert.(u"°C", f_tannul)),
        ))
        # nmr_params.csv (column names match what run_nmr.R expects)
        CSV.write(joinpath(_nmr_site_dir, "nmr_params.csv"), DataFrame(
            latitude    = [lat_dd],
            longitude   = [lon_dd],
            elevation_m = [elev_m],
            dstart      = [Dates.format(_sim_start, "dd/mm/yyyy")],
            dfinish     = [Dates.format(_sim_end,   "dd/mm/yyyy")],
            nyears      = [Dates.year(_sim_end) - Dates.year(_sim_start) + 1],
            snowtemp    = [ustrip(uconvert(u"°C", snow_temp_threshold))],
            snowdens    = [ustrip(uconvert(u"g/cm^3", snow_density_init))],
            snowmelt    = [snow_melt_factor],
            undercatch  = [undercatch],
            rainmelt    = [rain_melt_factor],
            densfun1    = [density_function[1]],
            densfun2    = [density_function[2]],
            densfun3    = [density_function[3]],
            densfun4    = [density_function[4]],
            snowcond    = [ustrip(snow_conductivity)],
            intercept   = [canopy_interception],
            REFL        = [albedo],
            SLE         = [emissivity],
            RUF         = [ustrip(roughness_height)],
        ))
        # nmr_soil.csv (10 NMR depth nodes — thermal props)
        CSV.write(joinpath(_nmr_site_dir, "nmr_soil.csv"), DataFrame(
            DEP            = [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200],
            Thcond         = _mc[_nmr_idx],
            SpecHeat       = _mhc[_nmr_idx],
            BulkDensity    = fill(ustrip(bulk_density), 10),
            MineralDensity = fill(ustrip(mineral_density), 10),
        ))
        # nmr_soil19.csv (19-node vectors for soil moisture model)
        CSV.write(joinpath(_nmr_site_dir, "nmr_soil19.csv"), DataFrame(
            node = 1:19,
            PE   = fill(ustrip(air_entry_potential), 19),
            KS   = fill(ustrip(sat_hydraulic_cond),  19),
            BB   = fill(campbell_b,                   19),
            BD   = fill(ustrip(bulk_density),         19),
            DD   = fill(ustrip(mineral_density),      19),
            L    = _rd,
        ))
        # nmr_initial.csv (initial conditions at 10 NMR nodes)
        CSV.write(joinpath(_nmr_site_dir, "nmr_initial.csv"), DataFrame(
            reshape(
                vcat(ustrip.(uconvert.(u"°C", initial_st[_nmr_idx])),
                     initial_sm[_nmr_idx]),
                1, 20),
            vcat(["ST$i" for i in 1:10], ["SM$i" for i in 1:10]),
        ))
    elseif reuse_nmr && isfile(_nmr_metout)
        println("  Reusing existing NMR outputs in $_nmr_site_dir")
    end
    return _nmr_site_dir, needs_run
end

# ── Run NicheMapR for many prepared sites concurrently ───────────────────────
# Each site's `Rscript run_nmr.R` only touches its own input/output directory,
# so they're safe to run as independent concurrent OS processes. Capped at
# `max_concurrent` (default: all available cores) since each is a real process
# with its own memory footprint — launching dozens at once would thrash.
function run_nmr_batch!(preps; max_concurrent = Sys.CPU_THREADS)
    to_run = Tuple{Int,String}[]
    for prep in preps
        nmr_dir, needs_run = write_nmr_inputs(prep)
        needs_run && push!(to_run, (prep.site_num, nmr_dir))
    end

    if isempty(to_run)
        println("\nNo NMR runs needed — all sites already have cached outputs.")
        return nothing
    end

    println("\nLaunching NMR runs for $(length(to_run)) site(s), up to $max_concurrent concurrent...")
    rscript = joinpath(@__DIR__, "run_nmr.R")
    queue = copy(to_run)
    running = Dict{Int, Base.Process}()
    print_lock = ReentrantLock()

    batch_time = @elapsed begin
        while !isempty(queue) || !isempty(running)
            while !isempty(queue) && length(running) < max_concurrent
                site_num, nmr_dir = popfirst!(queue)
                lock(print_lock) do
                    println("  Launching NMR for site $site_num...")
                end
                running[site_num] = run(`Rscript $rscript $nmr_dir`; wait=false)
            end
            sleep(0.5)
            for (site_num, proc) in collect(running)
                if process_exited(proc)
                    lock(print_lock) do
                        status = success(proc) ? "ok" : "FAILED (exit $(proc.exitcode))"
                        println("  NMR finished for site $site_num ($status)")
                    end
                    delete!(running, site_num)
                end
            end
        end
    end
    @printf("NMR batch total: %.2f s (%d site(s))\n", batch_time, length(to_run))
    return nothing
end

# ── Run the Julia model for many prepared sites in parallel ─────────────────
# `Microclimate.solve` allocates its own working state per call and shares
# only read-only parameters (`micro_model`) across problems — the same
# sharing pattern MicroclimateMapper's own internal parallel solve already
# relies on — so independent problems are safe to solve concurrently.
# Requires Julia started with multiple threads (`-t auto` or `-t N`) to
# actually run concurrently; falls back to sequential automatically otherwise.
# Returns a Dict{Int,Any} of site_num => (micro_out, julia_solve_time).
function run_julia_batch!(preps)
    n = length(preps)
    results = Vector{Tuple{Any,Float64}}(undef, n)
    print_lock = ReentrantLock()
    println("\nRunning $n Julia model(s) across $(Threads.nthreads()) thread(s)...")
    batch_time = @elapsed Threads.@threads for i in 1:n
        prep = preps[i]
        _cache_file = joinpath(sim_cache_dir, "$(prep.site_num)_$(prep.sim_start)_$(prep.sim_end).jls")
        if reuse_simulation && isfile(_cache_file)
            results[i] = (deserialize(_cache_file), NaN)
        else
            t = @elapsed result = Microclimate.solve(prep.problem)
            if cache_simulation
                mkpath(sim_cache_dir)
                serialize(_cache_file, result)
            end
            results[i] = (result, t)
        end
        lock(print_lock) do
            @printf("  [site %d] solved in %.2f s\n", prep.site_num, results[i][2])
        end
    end
    @printf("Julia batch total: %.2f s (%d site(s))\n", batch_time, n)
    return Dict(preps[i].site_num => results[i] for i in 1:n)
end

# ── Report results for one prepared site: NMR read, stats, save, plots ──────
# `make_plots`/`display_plots` let the multi-site batch skip/suppress plotting
# when not saving to disk; debug_site.jl always wants both.
function report_site_results(prep, micro_out, julia_solve_time;
                              plot_start = nothing, plot_end = nothing,
                              make_plots = true, display_plots = true)
    (; site_num, sim_start, sim_end, dates_vec, ndays, obs) = prep
    _sim_start, _sim_end = sim_start, sim_end

    # ── NicheMapR outputs (optional) ───────────────────────────────────────
    println("\nLooking for NicheMapR outputs in $nmr_out_dir ...")
    metout_df = try_read(joinpath(nmr_out_dir, string(site_num), "metout.csv"),    "metout.csv")
    soil_df   = try_read(joinpath(nmr_out_dir, string(site_num), "soil.csv"),      "soil.csv")
    smoist_df = try_read(joinpath(nmr_out_dir, string(site_num), "soilmoist.csv"), "soilmoist.csv")

    has_nmr = !isnothing(metout_df) && !isnothing(soil_df) && !isnothing(smoist_df)
    nmr = nothing
    if has_nmr
        nhours_nmr = ndays * 24
        nmr = DataFrame(
            DOY      = Int.(metout_df.DOY[1:nhours_nmr]),
            HOUR     = Int.(round.(metout_df.TIME[1:nhours_nmr] ./ 60.0)),
            SNOWDEP  = Float64.(metout_df.SNOWDEP[1:nhours_nmr]),
            SNOWDENS = Float64.(metout_df.SNOWDENS[1:nhours_nmr]),
            D5cm     = Float64.(soil_df.D5cm[1:nhours_nmr]),
            D20cm    = Float64.(soil_df.D20cm[1:nhours_nmr]),
            D50cm    = Float64.(soil_df.D50cm[1:nhours_nmr]),
            WC5cm    = Float64.(smoist_df.WC5cm[1:nhours_nmr]),
            WC20cm   = Float64.(smoist_df.WC20cm[1:nhours_nmr]),
            WC50cm   = Float64.(smoist_df.WC50cm[1:nhours_nmr]),
        )
        for (src_df, col) in [(soil_df, :D10cm), (soil_df, :D100cm),
                               (smoist_df, :WC10cm), (smoist_df, :WC100cm)]
            col in propertynames(src_df) &&
                (nmr[!, col] = Float64.(src_df[!, col][1:nhours_nmr]))
        end
        println("  NMR data loaded ($nhours_nmr rows)")
    end
    _nmr_timing_file = joinpath(nmr_out_dir, string(site_num), "nmr_timing.csv")
    nmr_solve_time = if isfile(_nmr_timing_file)
        Float64(DataFrame(CSV.File(_nmr_timing_file))[1, :elapsed_s])
    else
        NaN
    end

    # ── Extract Julia outputs ───────────────────────────────────────────────
    # Fine-node indices (1-based) for comparison depths:
    # depths (cm): 0, 1.25, 2.5, 3.75, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 40, 50, 75, 100, 150, 200
    #  index:       1  2     3    4     5  6     7   8     9   10    11  12  13  14  15  16  17   18   19
    di = Dict("5cm" => 5, "10cm" => 7, "20cm" => 11, "50cm" => 15, "100cm" => 17)

    nhours  = ndays * 24
    # Built from `dates_vec` (the model's own real-calendar day sequence,
    # leap years included), not `_sim_start`+offset — keeps alignment exact.
    t_model = [DateTime(d) + Hour(h) for d in dates_vec for h in 0:23]

    julia_D5cm    = ustrip.(u"°C", micro_out.soil_temperature[:, di["5cm"]])
    julia_D20cm   = ustrip.(u"°C", micro_out.soil_temperature[:, di["20cm"]])
    julia_D50cm   = ustrip.(u"°C", micro_out.soil_temperature[:, di["50cm"]])
    julia_D100cm  = ustrip.(u"°C", micro_out.soil_temperature[:, di["100cm"]])
    julia_WC5cm   = micro_out.soil_moisture[:, di["5cm"]]
    julia_WC20cm  = micro_out.soil_moisture[:, di["20cm"]]
    julia_WC50cm  = micro_out.soil_moisture[:, di["50cm"]]
    julia_WC100cm = micro_out.soil_moisture[:, di["100cm"]]
    julia_snow    = ustrip.(u"cm",     micro_out.snow_depth)
    julia_dens    = ustrip.(u"g/cm^3", micro_out.snow_density)
    julia_swe     = julia_snow .* julia_dens   # cm water equivalent

    # ── Statistics: Julia vs NicheMapR ────────────────────────────────────
    if has_nmr
        nr = nmr[1:nhours, :]
        println("\n── Julia vs NicheMapR (hourly, full period) ──")
        for (label, jv, nv) in [
                ("Snow depth (cm)",      julia_snow,   Float64.(nr.SNOWDEP)),
                ("SWE (cm)",             julia_swe,    Float64.(nr.SNOWDEP) .* Float64.(nr.SNOWDENS)),
                ("Soil temp 5 cm (°C)",   julia_D5cm,    Float64.(nr.D5cm)),
                ("Soil temp 20 cm (°C)",  julia_D20cm,   Float64.(nr.D20cm)),
                ("Soil temp 50 cm (°C)",  julia_D50cm,   Float64.(nr.D50cm)),
                ("Soil temp 100 cm (°C)", julia_D100cm,  :D100cm in propertynames(nr) ? Float64.(nr.D100cm) : fill(NaN, nhours)),
                ("Soil moist 5 cm",       julia_WC5cm,   Float64.(nr.WC5cm)),
                ("Soil moist 20 cm",      julia_WC20cm,  Float64.(nr.WC20cm)),
                ("Soil moist 50 cm",      julia_WC50cm,  Float64.(nr.WC50cm)),
                ("Soil moist 100 cm",     julia_WC100cm, :WC100cm in propertynames(nr) ? Float64.(nr.WC100cm) : fill(NaN, nhours)),
            ]
            s = compute_stats(nv, jv)
            println("  $(rpad(label, 22))  $(fmt_stat(s))")
        end
    end

    # ── Align models to obs timestamps and compute stats vs obs ───────────
    julia_cols = Dict{Symbol, Vector{Float64}}(
        :snow    => julia_snow,
        :swe     => julia_swe,
        :D5cm    => julia_D5cm,
        :D20cm   => julia_D20cm,
        :D50cm   => julia_D50cm,
        :D100cm  => julia_D100cm,
        :WC5cm   => julia_WC5cm   .* 100.0,   # m³/m³ → % to match obs
        :WC20cm  => julia_WC20cm  .* 100.0,
        :WC50cm  => julia_WC50cm  .* 100.0,
        :WC100cm => julia_WC100cm .* 100.0,
    )
    jal = align_model_to_obs(obs.DateTime, t_model, julia_cols)

    nal = nothing
    if has_nmr
        nr = nmr[1:nhours, :]
        nmr_cols = Dict{Symbol, Vector{Float64}}(
            :snow    => Float64.(nr.SNOWDEP),
            :swe     => Float64.(nr.SNOWDEP) .* Float64.(nr.SNOWDENS),
            :D5cm    => Float64.(nr.D5cm),
            :D20cm   => Float64.(nr.D20cm),
            :D50cm   => Float64.(nr.D50cm),
            :D100cm  => :D100cm in propertynames(nr) ? Float64.(nr.D100cm) : fill(NaN, nhours),
            :WC5cm   => Float64.(nr.WC5cm)   .* 100.0,
            :WC20cm  => Float64.(nr.WC20cm)  .* 100.0,
            :WC50cm  => Float64.(nr.WC50cm)  .* 100.0,
            :WC100cm => :WC100cm in propertynames(nr) ? Float64.(nr.WC100cm) .* 100.0 : fill(NaN, nhours),
        )
        nal = align_model_to_obs(obs.DateTime, t_model, nmr_cols)
    end

    obs_col_map = [
        ("Snow depth (cm)",      :snow,   :SNWD_cm),
        ("SWE (cm)",             :swe,    :WTEQ_cm),
        ("Soil temp 5 cm (°C)",   :D5cm,   :STO_5cm),
        ("Soil temp 20 cm (°C)",  :D20cm,  :STO_20cm),
        ("Soil temp 50 cm (°C)",  :D50cm,  :STO_50cm),
        ("Soil temp 100 cm (°C)", :D100cm, :STO_100cm),
        ("Soil moist 5 cm (%)",   :WC5cm,  :SMS_5cm),
        ("Soil moist 20 cm (%)",  :WC20cm, :SMS_20cm),
        ("Soil moist 50 cm (%)",  :WC50cm, :SMS_50cm),
        ("Soil moist 100 cm (%)", :WC100cm,:SMS_100cm),
    ]

    println("\n── Model vs observations ──")
    hdr = rpad("Variable", 24) * "  " * rpad("Julia", 54)
    !isnothing(nal) && (hdr *= "  NMR")
    println(hdr);  println("-"^length(hdr))

    site_rows = []
    for (label, sym, obscol) in obs_col_map
        ov = obscol in propertynames(obs) ? obs[!, obscol] : fill(missing, nrow(obs))
        sj = compute_stats(ov, jal[sym])
        sn = !isnothing(nal) ? compute_stats(ov, nal[sym]) : ModelStats(NaN, NaN, NaN, 0)
        row_str = rpad(label, 24) * "  " * rpad(fmt_stat(sj), 54)
        !isnothing(nal) && (row_str *= "  " * fmt_stat(sn))
        println(row_str)
        push!(site_rows, (site=site_num, variable=label,
            j_r=sj.r, j_rmse=sj.rmse, j_bias=sj.bias, j_n=sj.n,
            nmr_r=sn.r, nmr_rmse=sn.rmse, nmr_bias=sn.bias, nmr_n=sn.n))
    end
    if isfinite(julia_solve_time) || isfinite(nmr_solve_time)
        @printf("\n  Solver runtime — Julia: %.2f s   NMR: %s\n",
            julia_solve_time,
            isfinite(nmr_solve_time) ? @sprintf("%.2f s", nmr_solve_time) : "—")
    end

    if save_outputs
        mkpath(joinpath(outputs_dir, "stats"))
        CSV.write(joinpath(outputs_dir, "stats", "$(site_num)_$(_sim_start)_$(_sim_end).csv"),
                  DataFrame(site_rows))

        # ── Combined time-series CSV (obs + Julia + NMR at hourly resolution) ─
        mkpath(joinpath(outputs_dir, "data"))
        _data = DataFrame(DateTime = t_model)
        _data[!, :julia_snow_cm]    = julia_snow
        _data[!, :julia_swe_cm]     = julia_swe
        _data[!, :julia_D5cm_C]     = julia_D5cm
        _data[!, :julia_D20cm_C]    = julia_D20cm
        _data[!, :julia_D50cm_C]    = julia_D50cm
        _data[!, :julia_D100cm_C]   = julia_D100cm
        _data[!, :julia_WC5cm]      = julia_WC5cm
        _data[!, :julia_WC20cm]     = julia_WC20cm
        _data[!, :julia_WC50cm]     = julia_WC50cm
        _data[!, :julia_WC100cm]    = julia_WC100cm
        if has_nmr
            _data[!, :nmr_snow_cm]  = Float64.(nmr[1:nhours, :SNOWDEP])
            _data[!, :nmr_swe_cm]   = Float64.(nmr[1:nhours, :SNOWDEP]) .* Float64.(nmr[1:nhours, :SNOWDENS])
            _data[!, :nmr_D5cm_C]   = Float64.(nmr[1:nhours, :D5cm])
            _data[!, :nmr_D20cm_C]  = Float64.(nmr[1:nhours, :D20cm])
            _data[!, :nmr_D50cm_C]  = Float64.(nmr[1:nhours, :D50cm])
            _data[!, :nmr_D100cm_C] = :D100cm in propertynames(nmr) ? Float64.(nmr[1:nhours, :D100cm]) : fill(missing, nhours)
            _data[!, :nmr_WC5cm]    = Float64.(nmr[1:nhours, :WC5cm])
            _data[!, :nmr_WC20cm]   = Float64.(nmr[1:nhours, :WC20cm])
            _data[!, :nmr_WC50cm]   = Float64.(nmr[1:nhours, :WC50cm])
            _data[!, :nmr_WC100cm]  = :WC100cm in propertynames(nmr) ? Float64.(nmr[1:nhours, :WC100cm]) : fill(missing, nhours)
        end
        # Map obs to hourly model grid (obs are sparse; most hours will be missing)
        _obs_col_map = [(:SNWD_cm, :obs_snow_cm), (:WTEQ_cm, :obs_swe_cm),
                        (:STO_5cm, :obs_D5cm_C),  (:STO_20cm, :obs_D20cm_C),
                        (:STO_50cm, :obs_D50cm_C), (:STO_100cm, :obs_D100cm_C),
                        (:SMS_5cm, :obs_WC5cm),   (:SMS_20cm, :obs_WC20cm),
                        (:SMS_50cm, :obs_WC50cm), (:SMS_100cm, :obs_WC100cm)]
        _t_idx = Dict(t => i for (i, t) in enumerate(t_model))
        for (obscol, datacol) in _obs_col_map
            _v = Vector{Union{Float64,Missing}}(missing, nhours)
            if obscol in propertynames(obs)
                for r in eachrow(obs)
                    idx = get(_t_idx, round(r.DateTime, Hour(1)), 0)
                    idx > 0 && !ismissing(r[obscol]) && (_v[idx] = Float64(r[obscol]))
                end
            end
            _data[!, datacol] = _v
        end
        CSV.write(joinpath(outputs_dir, "data", "$(site_num).csv"), _data)

        # ── Cumulative timings CSV (one row per site, appended across runs) ──
        _timing_file = joinpath(outputs_dir, "timings.csv")
        _timing_row  = DataFrame(
            site          = [site_num],
            weather_source = [string(nameof(weather_source_choice))],
            sim_start     = [string(_sim_start)],
            sim_end       = [string(_sim_end)],
            ndays         = [ndays],
            julia_solve_s = [isfinite(julia_solve_time) ? julia_solve_time : missing],
            nmr_solve_s   = [isfinite(nmr_solve_time)   ? nmr_solve_time   : missing],
            cached        = [!isfinite(julia_solve_time)],
            run_at        = [string(now())],
        )
        if isfile(_timing_file)
            CSV.write(_timing_file, _timing_row; append=true)
        else
            CSV.write(_timing_file, _timing_row)
        end
    end

    # ── Plots ─────────────────────────────────────────────────────────────
    if make_plots
        t_ps = isnothing(plot_start) ? DateTime(_sim_start) : DateTime(plot_start)
        t_pe = isnothing(plot_end)   ? DateTime(_sim_end)   : DateTime(plot_end)
        hm   = findall(t -> t_ps <= t <= t_pe, t_model)
        t_w  = t_model[hm]
        om   = (obs.DateTime .>= t_ps) .& (obs.DateTime .<= t_pe)
        t_o  = obs.DateTime[om]
        sfx  = " — site $site_num ($(Dates.format(t_ps,"yyyy-mm-dd")) to $(Dates.format(t_pe,"yyyy-mm-dd")))"

        obs_vals(col) =   # safely extract a plotting vector (empty if column absent)
            Symbol(col) in propertynames(obs) ?
                replace(Float64.(coalesce.(obs[om, Symbol(col)], NaN)), NaN => NaN) :
                Float64[]

        # ── Snow ──────────────────────────────────────────────────────────────
        snow_ymax = max(maximum(julia_snow[hm]; init=0.0), has_nmr ? maximum(Float64.(nmr[hm, :SNOWDEP]); init=0.0) : 0.0) * 1.5
        snow_ymax = max(snow_ymax, 1.0)
        pa = plot(t_w, julia_snow[hm]; label="Julia", color=:black, lw=1.5,
                  title="Snow depth (cm)$sfx", ylabel="cm", ylims=(0, snow_ymax))
        has_nmr && plot!(pa, t_w, Float64.(nmr[hm, :SNOWDEP]);
            label="NicheMapR", color=:blue, lw=1, alpha=0.7)
        !isempty(obs_vals("SNWD_cm")) &&
            plot!(pa, t_o, obs_vals("SNWD_cm"); label="obs", color=:red, lw=1)

        swe_ymax = max(maximum(julia_swe[hm]; init=0.0), has_nmr ? maximum(Float64.(nmr[hm, :SNOWDEP]) .* Float64.(nmr[hm, :SNOWDENS]); init=0.0) : 0.0) * 1.5
        swe_ymax = max(swe_ymax, 1.0)
        pb = plot(t_w, julia_swe[hm]; label="Julia", color=:black, lw=1.5,
                  title="Snow water equiv. (cm w.e.)$sfx", ylabel="cm", ylims=(0, swe_ymax))
        has_nmr && plot!(pb, t_w, Float64.(nmr[hm, :SNOWDEP]) .* Float64.(nmr[hm, :SNOWDENS]);
            label="NicheMapR", color=:blue, lw=1, alpha=0.7)
        !isempty(obs_vals("WTEQ_cm")) &&
            plot!(pb, t_o, obs_vals("WTEQ_cm"); label="obs", color=:red, lw=1)

        snow_fig = plot(pa, pb; layout=(2,1), size=(900,500), link=:x,
            plot_title="Snow — site $site_num")
        display_plots && display(snow_fig)
        if save_outputs
            mkpath(joinpath(outputs_dir, "snow"))
            savefig(snow_fig, joinpath(outputs_dir, "snow", "$(site_num)_$(_sim_start)_$(_sim_end).$figure_format"))
        end

        # ── Soil temperature ──────────────────────────────────────────────────
        function make_temp_panel(julia_vec, nmr_col, obs_col, depth_label)
            lo, hi = extrema(julia_vec[hm])
            pad = max((hi - lo) * 0.5, 2.0)
            p = plot(t_w, julia_vec[hm]; label="Julia", color=:black, lw=1.5,
                     title="Soil temp $depth_label (°C)$sfx", ylabel="°C",
                     ylims=(lo - pad, hi + pad))
            has_nmr && nmr_col in names(nmr) &&
                plot!(p, t_w, Float64.(nmr[hm, Symbol(nmr_col)]);
                      label="NicheMapR", color=:blue, lw=1, alpha=0.7)
            !isempty(obs_vals(obs_col)) &&
                plot!(p, t_o, obs_vals(obs_col); label="obs", color=:red, lw=1)
            p
        end

        tp = [make_temp_panel(julia_D5cm,   "D5cm",   "STO_5cm",   "5 cm"),
              make_temp_panel(julia_D20cm,  "D20cm",  "STO_20cm",  "20 cm"),
              make_temp_panel(julia_D50cm,  "D50cm",  "STO_50cm",  "50 cm"),
              make_temp_panel(julia_D100cm, "D100cm", "STO_100cm", "100 cm")]
        temp_fig = plot(tp...; layout=(4,1), size=(900,900), link=:x,
            plot_title="Soil temperature — site $site_num")
        display_plots && display(temp_fig)
        if save_outputs
            mkpath(joinpath(outputs_dir, "soil_temperature"))
            savefig(temp_fig, joinpath(outputs_dir, "soil_temperature", "$(site_num)_$(_sim_start)_$(_sim_end).$figure_format"))
        end

        # ── Soil moisture ─────────────────────────────────────────────────────
        function make_moist_panel(julia_vec, nmr_col, obs_col, depth_label)
            jv_pct = julia_vec[hm] .* 100.0
            moist_ymax = max(maximum(jv_pct; init=0.0) * 1.5, 5.0)
            p = plot(t_w, jv_pct; label="Julia (%)", color=:black, lw=1.5,
                     title="Soil moisture $depth_label (%)$sfx", ylabel="%",
                     ylims=(0, moist_ymax))
            has_nmr && nmr_col in names(nmr) &&
                plot!(p, t_w, Float64.(nmr[hm, Symbol(nmr_col)]) .* 100.0;
                      label="NicheMapR (%)", color=:blue, lw=1, alpha=0.7)
            !isempty(obs_vals(obs_col)) &&
                plot!(p, t_o, obs_vals(obs_col); label="obs (%)", color=:red, lw=1)
            p
        end

        mp = [make_moist_panel(julia_WC5cm,   "WC5cm",   "SMS_5cm",   "5 cm"),
              make_moist_panel(julia_WC20cm,  "WC20cm",  "SMS_20cm",  "20 cm"),
              make_moist_panel(julia_WC50cm,  "WC50cm",  "SMS_50cm",  "50 cm"),
              make_moist_panel(julia_WC100cm, "WC100cm", "SMS_100cm", "100 cm")]
        moist_fig = plot(mp...; layout=(4,1), size=(900,900), link=:x,
            plot_title="Soil moisture — site $site_num")
        display_plots && display(moist_fig)
        if save_outputs
            mkpath(joinpath(outputs_dir, "soil_moisture"))
            savefig(moist_fig, joinpath(outputs_dir, "soil_moisture", "$(site_num)_$(_sim_start)_$(_sim_end).$figure_format"))
        end
    end

    return site_rows
end

# ── Single-site run: sequential composition of the phases above ─────────────
# Used by debug_site.jl. Real problems (bad data, solver failure, etc.) raise
# normally; the caller decides whether to catch them (comparison.jl wraps
# this per-site in try/catch) or let them propagate with a full stack trace
# (debug_site.jl calls this with nothing in between).
function run_site(site_num, meta_all, weather_cache, micro_model, soil_profile, mapper_model;
                   sim_start, sim_end, auto_date_range, max_sim_years,
                   plot_start = nothing, plot_end = nothing,
                   make_plots = true, display_plots = true)
    prep = prepare_site(site_num, meta_all, weather_cache, micro_model, soil_profile;
        sim_start, sim_end, auto_date_range, max_sim_years)
    isnothing(prep) && return nothing

    nmr_dir, needs_run = write_nmr_inputs(prep)
    if needs_run
        println("\nCalling run_nmr.R for site $site_num...")
        rscript = joinpath(@__DIR__, "run_nmr.R")
        run(`Rscript $rscript $nmr_dir`)
    end

    println("\nRunning Julia microclimate model ($(prep.ndays) days)...")
    _cache_file = joinpath(sim_cache_dir, "$(site_num)_$(prep.sim_start)_$(prep.sim_end).jls")
    julia_solve_time = NaN
    micro_out = if reuse_simulation && isfile(_cache_file)
        println("\nLoading cached simulation from $(_cache_file)...")
        deserialize(_cache_file)
    else
        result = nothing
        julia_solve_time = @elapsed result = Microclimate.solve(prep.problem)
        @printf("  Julia solver: %.2f s\n", julia_solve_time)
        if cache_simulation
            mkpath(sim_cache_dir)
            serialize(_cache_file, result)
            println("  Saved to $(_cache_file)")
        end
        result
    end

    return report_site_results(prep, micro_out, julia_solve_time;
        plot_start, plot_end, make_plots, display_plots)
end
