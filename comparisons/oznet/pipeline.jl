# pipeline.jl — shared model setup + per-site pipeline for comparison.jl
# (multi-site) and debug_site.jl (single-site, flat, no try/catch — for
# stepping through in the REPL/debugger).
#
# `depths` is fixed across sites (NMR_DEP19_CM, shared with NicheMapR's own
# 19-node hydraulic scheme -- see its definition), but soil hydraulic
# properties still come from a live per-site SLGA fetch + pedotransfer model
# rather than a hand-specified literature profile, so there's no
# site-independent `soil_profile` to build once up front.
#
# Pipeline phases, each usable standalone or via `run_site` (sequential,
# single-site composition of the same pieces):
#   resolve_site   — metadata, date-range resolution
#   prepare_site   — weather (fetched per-site; no cross-site batching here,
#                    unlike scan_snotel — OzNet is 38 sites, not 600+), SLGA
#                    soil profile, obs, initial conditions, builds the
#                    MicroProblem
#   run_julia_batch!  — Threads.@threads over many prepared sites
#   run_nmr_batch!    — concurrent Rscript processes, capped
#   report_site_results — NMR read, stats, save, plots for one site
#
# Requires config.jl and utils.jl to already be included.

# ── OPeNDAP point-query weather fetch for SILO / BARRA ───────────────────────
# MicroclimateMapper.PointQuery fetches via PointDataSources.jl's `getpoint`
# (both area- and points-mode) instead of file-download-based getraster.
if use_opendap_points
    MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()
    MicroclimateMapper.loader(::Type{<:BARRA}) = MicroclimateMapper.PointQuery()
end

# ── Site-independent model setup (except depths, which vary per site) ───────
function build_micro_model(depths)
    soil_properties_model = Microclimate.example_soil_properties_model()
    soil_hydraulic_model  = Microclimate.example_soil_hydraulic_model(;
        infiltration_algorithm = infiltration_algorithm_choice)

    micro_config = MicroConfig(;
        convergence            = convergence_choice,
        rainfall_schedule      = rainfall_schedule_choice,
        soil_moisture_strategy = soil_moisture_strategy_choice,
    )

    micro_model = MicroModel(;
        hours   = collect(0.0:1:23.0),
        depths,
        heights = profile_heights,
        soil_properties_model,
        soil_hydraulic_model,
        snow_model = NoSnow(),   # OzNet (Murrumbidgee catchment) doesn't see persistent snow
        config = micro_config,
    )

    mapper_model = MicroMapModel(;
        micro_model,
        dem_source              = dem_source_choice,
        weather_source           = weather_source_choice,
        surface_albedo_source    = albedo,
        roughness_height_source  = roughness_height,
        compute_terrain          = compute_terrain_choice,
    )

    return micro_model, mapper_model
end

# Parses a quoted comma list like "\"4,15,45,75\"" (as stored in
# oznetsiteinfo2.csv) into a Vector{Float64}.
_parse_depth_list(s::AbstractString) = parse.(Float64, split(strip(s, '"'), ","))

# NicheMapR's Fortran temperature solver hardcodes exactly 10 depth nodes
# (Numtyps, soilinit). Its soil hydraulic parameters (PE/KS/BB/BD/DD/L) use a
# finer 19-node scheme instead: the 10 real nodes plus a midpoint between
# each consecutive pair (odd 1-indexed positions = the real nodes, matching
# NicheMapR's own convention -- see e.g. micro_silo.R's
# `BulkDensity <- BD[seq(1,19,2)]`). Julia's MicroModel has no fixed node
# count, so both models now share this same 19-node scheme directly --
# nmr_soil.csv needs no separate interpolation step.
const NMR_DEP_CM   = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]
const NMR_DEP19_CM = let d = NMR_DEP_CM
    vcat([d[1]], reduce(vcat, [[(d[i] + d[i+1]) / 2, d[i+1]] for i in 1:(length(d) - 1)]))
end

# ── Site metadata + date-range resolution ────────────────────────────────────
# Returns `nothing` for expected skip conditions (no usable obs/weather date
# overlap) — routine, not an error.
function resolve_site(row; sim_start, sim_end, auto_date_range, max_sim_years)
    site_name = row.name
    lon_dd, lat_dd = Float64(row.lon), Float64(row.lat)
    obsfile = joinpath(OZNET_DATA_DIR, row.file)
    temp_depths_cm  = _parse_depth_list(row["temp.depths.cm"])
    moist_depths_cm = _parse_depth_list(row["soil.moisture.depths.cm"])
    sim_depths_cm   = NMR_DEP19_CM

    println("  $site_name — $(round(lat_dd, digits=4))°N  $(round(lon_dd, digits=4))°E  " *
            "($(length(sim_depths_cm)) depth nodes)")

    _sim_start, _sim_end = sim_start, sim_end
    weather_start = WEATHER_SOURCE_START[weather_source_choice]
    if auto_date_range
        isfile(obsfile) || (println("  Skipping — obs file not found: $obsfile"); return nothing)
        first_line_dates = String[]
        open(obsfile) do io
            readline(io); readline(io)   # skip site name + header lines
            for line in eachline(io)
                push!(first_line_dates, first(split(line)))
            end
        end
        isempty(first_line_dates) && (println("  Skipping — no data rows."); return nothing)
        obs_dates = Date.(first_line_dates, dateformat"d/m/y")
        _sim_end   = min(maximum(obs_dates), today() - Day(1))
        _sim_start = max(_sim_end - Year(max_sim_years) + Day(1), minimum(obs_dates), weather_start)
        if _sim_start >= _sim_end
            @warn "Site $site_name: no usable obs/weather overlap — skipping."
            return nothing
        end
        println("  Auto date range: $_sim_start to $_sim_end")
    end

    # Snap to full calendar years -- see scan_snotel/pipeline.jl for why
    # (MicroclimateMapper.jl sizes its annual-mean deep-soil climatology off
    # the full Jan1-Dec31 span of the touched year(s)).
    _sim_start = max(Date(Dates.year(_sim_start), 1, 1), weather_start)
    _sim_end   = min(Date(Dates.year(_sim_end), 12, 31), today() - Day(1))

    return (; site_name, lon_dd, lat_dd, obsfile, temp_depths_cm, moist_depths_cm, sim_depths_cm,
              sim_start = _sim_start, sim_end = _sim_end)
end

# Disk path for one site's cached weather forcing, namespaced by weather
# source so switching weather_source_choice can't silently reuse another
# source's cached forcing.
_weather_cache_path(site_name, sim_start, sim_end) =
    joinpath(weather_cache_dir, "$(nameof(weather_source_choice))_$(site_name)_$(sim_start)_$(sim_end).jls")

# ── Prepare one site: weather, SLGA soil profile, obs, initial conditions ───
# Does everything up through building the final `MicroProblem`, but doesn't
# run anything (no NMR launch, no Julia solve) — that's `run_nmr_batch!`/
# `run_julia_batch!`, or `run_site`'s synchronous single-site composition.
# Returns `nothing` for the same routine skip conditions as `resolve_site`.
function prepare_site(row, weather_cache; sim_start, sim_end, auto_date_range, max_sim_years)
    println("\n" * "="^72)
    println("Site $(row.name)")
    println("="^72)

    resolved = resolve_site(row; sim_start, sim_end, auto_date_range, max_sim_years)
    isnothing(resolved) && return nothing
    (; site_name, lon_dd, lat_dd, obsfile, temp_depths_cm, moist_depths_cm, sim_depths_cm) = resolved
    _sim_start = resolved.sim_start
    _sim_end   = resolved.sim_end

    depths = (sim_depths_cm ./ 100.0) .* u"m"
    micro_model, mapper_model = build_micro_model(depths)

    # ── Weather forcing ───────────────────────────────────────────────────
    dates_range = _sim_start:Day(1):_sim_end
    wc_key = (nameof(weather_source_choice), site_name, _sim_start, _sim_end)
    _weather_disk_cache = _weather_cache_path(site_name, _sim_start, _sim_end)
    if !haskey(weather_cache, wc_key) && reuse_weather && isfile(_weather_disk_cache)
        weather_cache[wc_key] = deserialize(_weather_disk_cache)
        println("\nLoaded weather forcing from disk cache: $_weather_disk_cache")
    end
    if !haskey(weather_cache, wc_key)
        println("\nFetching weather forcing ($_sim_start to $_sim_end)...")
        point = GIW.Point((lon_dd, lat_dd))
        vec_problem = MicroVectorProblem(;
            model  = mapper_model,
            points = [point],
            dates  = dates_range,
            soil_profile = SoilProfile(; bulk_density = fill(1.3u"Mg/m^3", length(depths)),
                mineral_density = fill(2.56u"Mg/m^3", length(depths)),
                mineral_conductivity = fill(1.25u"W/m/K", length(depths)),
                mineral_heat_capacity = fill(870.0u"J/kg/K", length(depths)),
                hydraulics = CampbellHydraulicProfile(;
                    air_entry_water_potential = fill(1.0u"J/kg", length(depths)),
                    saturated_hydraulic_conductivity = fill(0.001u"kg*s/m^3", length(depths)),
                    campbell_b_parameter = fill(4.5, length(depths)),
                    root_density = fill(1.0u"m/m^3", length(depths)))),
            # Placeholder soil_profile/init: only used to size the fetch, not
            # to inform it -- weather forcing doesn't depend on soil
            # properties. The real, SLGA-derived soil_profile is built below.
            init = (; soil_moisture = fill(0.2, length(depths))),
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

    dates_vec = MicroclimateMapper._normalise_dates(dates_range)
    f_doy    = dayofyear.(dates_vec)
    ndays    = length(f_doy)
    em = base_inputs.environment_minmax
    f_tminn  = _minmax_series(em, :reference_temperature, :reference_temperature_min)
    f_tmaxx  = _minmax_series(em, :reference_temperature, :reference_temperature_max)
    f_rhminn = _minmax_series(em, :reference_humidity, :reference_humidity_min)
    f_rhmaxx = _minmax_series(em, :reference_humidity, :reference_humidity_max)
    f_ccmaxx = _minmax_series(em, :cloud_cover, :cloud_cover_max)
    f_ccminn = _minmax_series(em, :cloud_cover, :cloud_cover_min)
    f_wnminn = _minmax_series(em, :reference_wind_speed, :reference_wind_speed_min)
    f_wnmaxx = _minmax_series(em, :reference_wind_speed, :reference_wind_speed_max)
    f_rain   = base_inputs.environment_daily.rainfall
    f_tannul = base_inputs.environment_daily.deep_soil_temperature
    println("  $ndays days from $(nameof(weather_source_choice))")

    # ── SLGA soil profile (live fetch + Cosby pedotransfer) ─────────────────
    println("\nFetching SLGA soil texture...")
    b = soil_area_buffer_deg
    area = Extent(X = (lon_dd - b, lon_dd + b), Y = (lat_dd - b, lat_dd + b))
    soil_result = build_soil_profile(SLGA, area; depths, pedotransfer_model = pedotransfer_model_choice)
    soil_profile = soil_result.soil_profile

    # Organic litter/crust cap on the top ~2.5 cm (nodes 1:3 of NMR_DEP19_CM:
    # 0, 1.25, 2.5 cm), matching NicheMapR's own default `cap=1` behaviour
    # (see micro_global.R) -- overrides the SLGA-derived mineral thermal
    # properties there with low-conductivity, high-heat-capacity litter
    # values, damping the diurnal amplitude near the surface.
    soil_profile.mineral_conductivity[1:3] .= 0.2u"W/m/K"
    soil_profile.mineral_heat_capacity[1:3] .= 1920.0u"J/kg/K"

    # ── Observations ───────────────────────────────────────────────────────
    println("\nReading observations from $obsfile...")
    obs = read_oznet_obs(obsfile, temp_depths_cm, moist_depths_cm)
    omask = (obs.DateTime .>= DateTime(_sim_start)) .& (obs.DateTime .<= DateTime(_sim_end, Time(23, 0, 0)))
    obs = obs[omask, :]
    println("  $(nrow(obs)) obs rows")

    # ── Initial conditions: interpolate from the shallowest obs + deep-soil
    #    annual mean anchor (same approach as scan_snotel) ───────────────────
    t_air_mean = f_tannul[1]
    t_air_day1 = (f_tminn[1] + f_tmaxx[1]) / 2
    first_obs_rows = obs[obs.DateTime .>= DateTime(_sim_start), :]
    first_obs = nrow(first_obs_rows) > 0 ? first_obs_rows[1, :] : nothing

    st_known_d = Float64[0.0]
    st_known_T = typeof(t_air_day1)[t_air_day1]
    if !isnothing(first_obs)
        for d in temp_depths_cm
            col = Symbol("Temp_$(d)cm")
            v = first_obs[col]
            if !ismissing(v)
                push!(st_known_d, d)
                push!(st_known_T, (Float64(v) + 273.15) * u"K")
            end
        end
    end
    push!(st_known_d, sim_depths_cm[end])
    push!(st_known_T, t_air_mean)
    initial_st = [_linterp(st_known_d, st_known_T, d) for d in sim_depths_cm]

    sm_known_d = Float64[]; sm_known_v = Float64[]
    if !isnothing(first_obs)
        for d in moist_depths_cm
            col = Symbol("SM_$(d)cm")
            v = first_obs[col]
            if !ismissing(v)
                vv = Float64(v)
                vv > 1.0 && (vv /= 100.0)   # guard: % -> fraction if needed
                push!(sm_known_d, d)
                push!(sm_known_v, clamp(vv, 0.02, 0.5))
            end
        end
    end
    push!(sm_known_d, sim_depths_cm[end])
    push!(sm_known_v, 0.3)
    initial_sm = length(sm_known_d) >= 2 ?
        clamp.([_linterp(sm_known_d, sm_known_v, d) for d in sim_depths_cm], 0.02, 0.5) :
        fill(0.3, length(sim_depths_cm))

    println("  Initial ST: " * join([@sprintf("%.1f", ustrip(uconvert(u"°C", T))) for T in initial_st], ", ") * " °C")
    println("  Initial SM: " * join([@sprintf("%.3f", v) for v in initial_sm], ", "))

    inputs = MicroInputs(;
        site               = base_inputs.site,
        soil_profile,
        environment_minmax = base_inputs.environment_minmax,
        environment_daily  = base_inputs.environment_daily,
        environment_hourly = base_inputs.environment_hourly,
        initial_soil_temperature = initial_st,
        initial_soil_moisture    = initial_sm,
    )

    problem = MicroProblem(micro_model, inputs;
        days      = f_doy,
        time_mode = ConsecutiveDayMode(; spinup_first_day = false),
    )

    return (; site_name, lon_dd, lat_dd, elev_m = ustrip(u"m", base_inputs.site.elevation),
              sim_start = _sim_start, sim_end = _sim_end,
              temp_depths_cm, moist_depths_cm, sim_depths_cm,
              dates_vec, f_doy, ndays,
              f_tminn, f_tmaxx, f_rhminn, f_rhmaxx, f_ccmaxx, f_ccminn,
              f_wnminn, f_wnmaxx, f_rain, f_tannul,
              obs, initial_st, initial_sm, soil_result,
              problem)
end

# ── Write NMR input CSVs for one prepared site ───────────────────────────────
# Returns (nmr_site_dir, needs_run::Bool).
function write_nmr_inputs(prep)
    (; site_name, lon_dd, lat_dd, elev_m, sim_start, sim_end, f_doy,
       f_tminn, f_tmaxx, f_rhminn, f_rhmaxx, f_wnmaxx, f_rain, f_ccmaxx, f_ccminn,
       f_tannul, initial_st, initial_sm, soil_result) = prep
    _sim_start, _sim_end = sim_start, sim_end

    _nmr_site_dir = joinpath(nmr_out_dir, site_name)
    _nmr_metout   = joinpath(_nmr_site_dir, "metout.csv")
    needs_run = run_nmr && (!reuse_nmr || !isfile(_nmr_metout))
    if needs_run
        mkpath(_nmr_site_dir)
        sp = soil_result.soil_profile

        CSV.write(joinpath(_nmr_site_dir, "nmr_forcing.csv"), DataFrame(
            doy       = f_doy,
            Tmin_C    = ustrip.(uconvert.(u"°C", f_tminn)),
            Tmax_C    = ustrip.(uconvert.(u"°C", f_tmaxx)),
            RHmin_pct = f_rhminn .* 100.0,
            RHmax_pct = f_rhmaxx .* 100.0,
            Wind_ms   = ustrip.(f_wnmaxx),
            Rain_mm   = ustrip.(f_rain),
            CCmax_pct = f_ccmaxx .* 100.0,
            CCmin_pct = f_ccminn .* 100.0,
            tannul_C  = ustrip.(uconvert.(u"°C", f_tannul)),
        ))
        CSV.write(joinpath(_nmr_site_dir, "nmr_params.csv"), DataFrame(
            latitude    = [lat_dd],
            longitude   = [lon_dd],
            elevation_m = [elev_m],
            dstart      = [Dates.format(_sim_start, "dd/mm/yyyy")],
            dfinish     = [Dates.format(_sim_end,   "dd/mm/yyyy")],
            nyears      = [Dates.year(_sim_end) - Dates.year(_sim_start) + 1],
            REFL        = [albedo],
            SLE         = [emissivity],
            RUF         = [ustrip(roughness_height)],
        ))
        # 19 rows (NMR_DEP19_CM), straight from Julia's own soil profile --
        # it's already built at this exact scheme, no interpolation needed.
        # run_nmr.R subsets the odd (real-node) positions for the 10-node
        # temperature grid, matching NicheMapR's own pedotransfer() convention.
        CSV.write(joinpath(_nmr_site_dir, "nmr_soil.csv"), DataFrame(
            DEP            = NMR_DEP19_CM,
            BulkDensity    = ustrip.(u"Mg/m^3", sp.bulk_density),
            MineralDensity = ustrip.(u"Mg/m^3", sp.mineral_density),
            Thcond         = ustrip.(u"W/m/K", sp.mineral_conductivity),
            SpecHeat       = ustrip.(u"J/kg/K", sp.mineral_heat_capacity),
            PE             = ustrip.(u"J/kg", sp.hydraulics.air_entry_water_potential),
            KS             = ustrip.(u"kg*s/m^3", sp.hydraulics.saturated_hydraulic_conductivity),
            BB             = sp.hydraulics.campbell_b_parameter,
            L              = ustrip.(u"m/m^3", sp.hydraulics.root_density),
        ))
        # Odd positions of the 19-node initial_st/initial_sm are the 10 real
        # nodes (same convention as nmr_soil.csv) -- NMR's temperature/
        # moisture initial state is still 10-node.
        st10 = ustrip.(uconvert.(u"°C", initial_st))[1:2:end]
        sm10 = initial_sm[1:2:end]
        CSV.write(joinpath(_nmr_site_dir, "nmr_initial.csv"), DataFrame(
            reshape(vcat(st10, sm10), 1, 20),
            vcat(["ST$i" for i in 1:10], ["SM$i" for i in 1:10]),
        ))
    elseif reuse_nmr && isfile(_nmr_metout)
        println("  Reusing existing NMR outputs in $_nmr_site_dir")
    end
    return _nmr_site_dir, needs_run
end

# ── Run NicheMapR for many prepared sites concurrently ───────────────────────
function run_nmr_batch!(preps; max_concurrent = Sys.CPU_THREADS)
    to_run = Tuple{String,String}[]
    for prep in preps
        nmr_dir, needs_run = write_nmr_inputs(prep)
        needs_run && push!(to_run, (prep.site_name, nmr_dir))
    end

    if isempty(to_run)
        println("\nNo NMR runs needed — all sites already have cached outputs.")
        return nothing
    end

    println("\nLaunching NMR runs for $(length(to_run)) site(s), up to $max_concurrent concurrent...")
    rscript = joinpath(@__DIR__, "run_nmr.R")
    queue = copy(to_run)
    running = Dict{String, Base.Process}()
    print_lock = ReentrantLock()

    batch_time = @elapsed begin
        while !isempty(queue) || !isempty(running)
            while !isempty(queue) && length(running) < max_concurrent
                site_name, nmr_dir = popfirst!(queue)
                lock(print_lock) do
                    println("  Launching NMR for site $site_name...")
                end
                running[site_name] = run(`Rscript $rscript $nmr_dir`; wait=false)
            end
            sleep(0.5)
            for (site_name, proc) in collect(running)
                if process_exited(proc)
                    lock(print_lock) do
                        status = success(proc) ? "ok" : "FAILED (exit $(proc.exitcode))"
                        println("  NMR finished for site $site_name ($status)")
                    end
                    delete!(running, site_name)
                end
            end
        end
    end
    @printf("NMR batch total: %.2f s (%d site(s))\n", batch_time, length(to_run))
    return nothing
end

# ── Run the Julia model for many prepared sites in parallel ─────────────────
function run_julia_batch!(preps)
    n = length(preps)
    results = Vector{Tuple{Any,Float64}}(undef, n)
    print_lock = ReentrantLock()
    println("\nRunning $n Julia model(s) across $(Threads.nthreads()) thread(s)...")
    batch_time = @elapsed Threads.@threads for i in 1:n
        prep = preps[i]
        _cache_file = joinpath(sim_cache_dir, "$(prep.site_name)_$(prep.sim_start)_$(prep.sim_end).jls")
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
            @printf("  [site %s] solved in %.2f s\n", prep.site_name, results[i][2])
        end
    end
    @printf("Julia batch total: %.2f s (%d site(s))\n", batch_time, n)
    return Dict(preps[i].site_name => results[i] for i in 1:n)
end

# ── Report results for one prepared site: NMR read, stats, save, plots ──────
function report_site_results(prep, micro_out, julia_solve_time;
                              plot_start = nothing, plot_end = nothing,
                              make_plots = true, display_plots = true)
    (; site_name, sim_start, sim_end, dates_vec, ndays, obs,
       temp_depths_cm, moist_depths_cm, sim_depths_cm) = prep
    _sim_start, _sim_end = sim_start, sim_end

    println("\nLooking for NicheMapR outputs in $nmr_out_dir ...")
    soil_df   = try_read(joinpath(nmr_out_dir, site_name, "soil.csv"),      "soil.csv")
    smoist_df = try_read(joinpath(nmr_out_dir, site_name, "soilmoist.csv"), "soilmoist.csv")
    has_nmr = !isnothing(soil_df) && !isnothing(smoist_df)

    nhours  = ndays * 24
    t_model = [DateTime(d) + Hour(h) for d in dates_vec for h in 0:23]

    # Plot window -- stats below always use the full period; only the saved
    # figures are restricted to [plot_start, plot_end] when given.
    t_ps = isnothing(plot_start) ? DateTime(_sim_start) : DateTime(plot_start)
    t_pe = isnothing(plot_end)   ? DateTime(_sim_end)   : DateTime(plot_end)
    hm   = findall(t -> t_ps <= t <= t_pe, t_model)
    t_w  = t_model[hm]
    om   = (obs.DateTime .>= t_ps) .& (obs.DateTime .<= t_pe)

    # Julia now shares NMR's 19-node hydraulic scheme (NMR_DEP19_CM), but obs
    # depths (e.g. 4, 15, 45, 75 cm) don't generally land exactly on either
    # model's nodes, so both use nearest-node lookup. NMR's own output stays
    # on its 10-node temperature grid (NMR_DEP_CM), not the 19-node scheme --
    # see write_nmr_inputs. NMR's soil.csv/soilmoist.csv columns are
    # [DOY, TIME, <one per DEP node, in DEP order>] (NicheMapR::microrun.R)
    # -- indexed positionally rather than by parsed column name, since
    # name-based lookup would need to exactly reproduce R's
    # `paste0("D", DEP, "cm")` number-to-string formatting.
    node_idx(d) = argmin(abs.(sim_depths_cm .- d))
    nmr_node_idx(d) = argmin(abs.(NMR_DEP_CM .- d))

    println("\n── Model vs observations ──")
    hdr = rpad("Variable", 26) * "  " * rpad("Julia", 54)
    has_nmr && (hdr *= "  NMR")
    println(hdr); println("-"^length(hdr))

    site_rows = []
    make_plots && (temp_figs = Plots.Plot[]; moist_figs = Plots.Plot[])

    for d in temp_depths_cm
        i = node_idx(d)
        julia_vec = ustrip.(u"°C", micro_out.soil_temperature[:, i])
        obscol = Symbol("Temp_$(d)cm")
        jal = align_model_to_obs(obs.DateTime, t_model, Dict(:v => julia_vec))
        sj = compute_stats(obs[!, obscol], jal[:v])
        sn = ModelStats(NaN, NaN, NaN, 0)
        nmr_vec = Float64[]
        if has_nmr && nrow(soil_df) >= nhours
            nmr_vec = Float64.(soil_df[1:nhours, 2 + nmr_node_idx(d)])
            nal = align_model_to_obs(obs.DateTime, t_model, Dict(:v => nmr_vec))
            sn = compute_stats(obs[!, obscol], nal[:v])
        end
        label = "Soil temp $(d) cm (°C)"
        row_str = rpad(label, 26) * "  " * rpad(fmt_stat(sj), 54)
        has_nmr && (row_str *= "  " * fmt_stat(sn))
        println(row_str)
        push!(site_rows, (site=site_name, variable=label, depth_cm=d, kind="temperature",
            j_r=sj.r, j_rmse=sj.rmse, j_bias=sj.bias, j_n=sj.n,
            nmr_r=sn.r, nmr_rmse=sn.rmse, nmr_bias=sn.bias, nmr_n=sn.n))
        if make_plots
            p = plot(t_w, julia_vec[hm]; label="Julia", color=:black, lw=1.5, title=label, ylabel="°C")
            !isempty(nmr_vec) && plot!(p, t_w, nmr_vec[hm]; label="NicheMapR", color=:blue, lw=1, alpha=0.7)
            obsmask = om .& .!ismissing.(obs[!, obscol])
            any(obsmask) && plot!(p, obs.DateTime[obsmask], Float64.(obs[obsmask, obscol]); label="obs", color=:red, lw=1)
            push!(temp_figs, p)
        end
    end

    for d in moist_depths_cm
        i = node_idx(d)
        julia_vec = micro_out.soil_moisture[:, i] .* 100.0
        obscol = Symbol("SM_$(d)cm")
        jal = align_model_to_obs(obs.DateTime, t_model, Dict(:v => julia_vec))
        sj = compute_stats(obs[!, obscol], jal[:v])
        sn = ModelStats(NaN, NaN, NaN, 0)
        nmr_vec = Float64[]
        if has_nmr && nrow(smoist_df) >= nhours
            nmr_vec = Float64.(smoist_df[1:nhours, 2 + nmr_node_idx(d)]) .* 100.0
            nal = align_model_to_obs(obs.DateTime, t_model, Dict(:v => nmr_vec))
            sn = compute_stats(obs[!, obscol], nal[:v])
        end
        label = "Soil moisture $(d) cm (%vol)"
        row_str = rpad(label, 26) * "  " * rpad(fmt_stat(sj), 54)
        has_nmr && (row_str *= "  " * fmt_stat(sn))
        println(row_str)
        push!(site_rows, (site=site_name, variable=label, depth_cm=d, kind="moisture",
            j_r=sj.r, j_rmse=sj.rmse, j_bias=sj.bias, j_n=sj.n,
            nmr_r=sn.r, nmr_rmse=sn.rmse, nmr_bias=sn.bias, nmr_n=sn.n))
        if make_plots
            p = plot(t_w, julia_vec[hm]; label="Julia (%)", color=:black, lw=1.5, title=label, ylabel="%")
            !isempty(nmr_vec) && plot!(p, t_w, nmr_vec[hm]; label="NicheMapR (%)", color=:blue, lw=1, alpha=0.7)
            obsmask = om .& .!ismissing.(obs[!, obscol])
            any(obsmask) && plot!(p, obs.DateTime[obsmask], Float64.(obs[obsmask, obscol]); label="obs (%)", color=:red, lw=1)
            push!(moist_figs, p)
        end
    end

    if isfinite(julia_solve_time)
        @printf("\n  Julia solver runtime: %.2f s\n", julia_solve_time)
    end

    if save_outputs
        mkpath(joinpath(outputs_dir, "stats"))
        CSV.write(joinpath(outputs_dir, "stats", "$(site_name)_$(_sim_start)_$(_sim_end).csv"),
                  DataFrame(site_rows))
    end

    if make_plots
        temp_fig = plot(temp_figs...; layout=(length(temp_figs), 1), ylims=(-10, 60),
            size=(900, 300 * length(temp_figs)), link=:x, plot_title="Soil temperature — $site_name")
        moist_fig = plot(moist_figs...; layout=(length(moist_figs), 1), ylims=(0, 50),
            size=(900, 300 * length(moist_figs)), link=:x, plot_title="Soil moisture — $site_name")
        display_plots && display(temp_fig)
        display_plots && display(moist_fig)
        if save_outputs
            mkpath(joinpath(outputs_dir, "soil_temperature"))
            mkpath(joinpath(outputs_dir, "soil_moisture"))
            savefig(temp_fig, joinpath(outputs_dir, "soil_temperature", "$(site_name)_$(_sim_start)_$(_sim_end).$figure_format"))
            savefig(moist_fig, joinpath(outputs_dir, "soil_moisture", "$(site_name)_$(_sim_start)_$(_sim_end).$figure_format"))
        end
    end

    return site_rows
end

# ── Single-site run: sequential composition of the phases above ─────────────
function run_site(row, weather_cache;
                   sim_start, sim_end, auto_date_range, max_sim_years,
                   plot_start = nothing, plot_end = nothing,
                   make_plots = true, display_plots = true)
    prep = prepare_site(row, weather_cache; sim_start, sim_end, auto_date_range, max_sim_years)
    isnothing(prep) && return nothing

    nmr_dir, needs_run = write_nmr_inputs(prep)
    if needs_run
        println("\nCalling run_nmr.R for site $(prep.site_name)...")
        rscript = joinpath(@__DIR__, "run_nmr.R")
        run(`Rscript $rscript $nmr_dir`)
    end

    println("\nRunning Julia microclimate model ($(prep.ndays) days)...")
    _cache_file = joinpath(sim_cache_dir, "$(prep.site_name)_$(prep.sim_start)_$(prep.sim_end).jls")
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
        end
        result
    end

    return report_site_results(prep, micro_out, julia_solve_time;
        plot_start, plot_end, make_plots, display_plots)
end
