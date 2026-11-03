# pipeline.jl — shared model setup + per-site-year pipeline for comparison.jl
# (multi-site/multi-year) and single_site.jl (single site-year, flat, no
# try/catch -- for stepping through in the REPL/debugger).
#
# Forcing comes straight from each site's own tower file (30-min, aggregated
# to hourly, QC-filtered during aggregation), not a fetched grid -- unlike
# oznet there's no weather-source fetch phase here.
#
# Phases: resolve_site (metadata) -> prepare_site (read/aggregate/build
# MicroProblem) -> run_site (solve, pair against obs).
#
# Requires config.jl and utils.jl to already be included.

# Absolute humidity (g/m^3) + air temperature -> relative humidity (%), via
# FluidProperties' saturation vapour pressure (ideal-gas vapour density ->
# vapour pressure; no real-gas correction -- fine for a cross-check/secondary
# gap-fill source, not the primary model driver).
const _MW_WATER = 0.018015u"kg/mol"
function _rh_from_ah(ah_gm3, ta_C; vapour_pressure_equation=GoffGratch())
    T = u"K"(ta_C * u"°C")
    p_v = (ah_gm3 * 1e-3u"kg/m^3") * Unitful.R * T / _MW_WATER
    p_sat = FluidProperties.vapour_pressure(vapour_pressure_equation, T)
    return ustrip(u"Pa", p_v) / ustrip(u"Pa", p_sat) * 100.0
end

function resolve_site(site_name, years)
    path = ozflux_site_path(OZFLUX_DATA_DIR, site_name)
    available = ozflux_years(path)
    missing_years = setdiff(years, available)
    isempty(missing_years) || error("$path has no data for year(s) $missing_years (has $available)")

    global_attrib = site_attributes(site_name)
    canopy_height = site_canopy_height(site_name)
    latitude  = parse(Float64, string(global_attrib["latitude"])) * u"°"
    longitude = parse(Float64, string(global_attrib["longitude"])) * u"°"

    forcing = forcing_sources(site_name, years, canopy_height)
    reference_height = forcing.reference_height
    canopy_height < reference_height || error("$site_name's canopy ($canopy_height) reaches the reanalysis " *
        "forcing height ($reference_height)")

    return (; site_name, years, path, canopy_height, reference_height, forcing,
              forcing_backups = forcing_backups(site_name),
              latitude, longitude, global_attrib,
              elevation = site_elevation(site_name, latitude, longitude),
              leaf_area_index = site_leaf_area_index(site_name, latitude, longitude, years),
              leaf_reflectance = site_setting(site_name, "leaf_reflectance"),
              leaf_transmittance = site_setting(site_name, "leaf_transmittance"),
              height_series = discover_height_series(site_name, years, forcing, forcing_backups(site_name)))
end

# 30-min -> hourly: mean of the QC-good sub-hourly readings; `missing` if none
# were good. QC-filtering happens here so downstream code only ever checks
# for `missing`, not a separate flag column. OzFlux is END-of-period labeled
# (an HH:30 reading covers [HH:00,HH:30), an HH:00 reading covers
# [HH-1:30,HH:00)) -- but the STARTING PHASE differs per file, confirmed by
# direct inspection: CapeTribulation's raw series starts exactly at 00:00,
# Whroo/Calperum's starts at 00:30. A fixed positional stride from row 1
# (this function's previous approach) mislabels every hour for whichever
# phase it wasn't tested against -- silently, since it still produces a
# same-length, monotonic DateTime column, just offset by half an hour from
# the true hour-start label. Hour labels are computed per-row instead and
# grouped, which is correct regardless of starting phase.
_hour_group_start(t) = Dates.minute(t) == 30 ? floor(t, Dates.Hour) : floor(t, Dates.Hour) - Dates.Hour(1)

function _aggregate_hourly(df)
    n_per_hour = Minute(60) ÷ DataFrames.metadata(df, "time_step")
    groups = Dict{DateTime,Vector{Int}}()
    for (i, t) in enumerate(df.DateTime)
        push!(get!(groups, _hour_group_start(t), Int[]), i)
    end
    # Drop any partial hour at the edges (fewer than n_per_hour sub-steps) --
    # an incomplete average would silently bias that hour rather than
    # reflecting a real QC gap; only affects the very first/last hour of the
    # requested span, or a gap in the raw time index.
    hour_starts = filter(h -> length(groups[h]) == n_per_hour, sort(collect(keys(groups))))
    out = DataFrame(DateTime = hour_starts)
    value_cols = [n for n in names(df) if n != "DateTime" && !endswith(n, "_QCFlag")]
    for name in value_cols
        col = df[!, name]
        qc_name = name * "_QCFlag"
        has_qc = qc_name in names(df)
        hourly = Vector{Union{Float64,Missing}}(missing, length(hour_starts))
        @inbounds for (k, h) in enumerate(hour_starts)
            vals = Float64[]
            for r in groups[h]
                v = col[r]
                (ismissing(v) || (has_qc && !(df[r, qc_name] in OZFLUX_GOOD_QC_FLAGS))) && continue
                push!(vals, v)
            end
            isempty(vals) || (hourly[k] = mean(vals))
        end
        out[!, name] = hourly
    end
    return out
end

# Linear interpolation across `missing` runs (edges held at the nearest known
# value); O(n). Forcing series must be gap-free for the solver -- comparison
# targets are left with `missing` instead (compute_stats skips those pairs).
function fill_gaps(v::AbstractVector{<:Union{T,Missing}}) where T
    known_idx = findall(!ismissing, v)
    isempty(known_idx) && error("all values missing -- cannot gap-fill")
    known_val = T[v[i] for i in known_idx]
    n = length(v)
    out = Vector{T}(undef, n)
    k = 1
    for i in 1:n
        if !ismissing(v[i])
            out[i] = v[i]
            continue
        end
        while k < length(known_idx) && known_idx[k+1] < i
            k += 1
        end
        if i < known_idx[1]
            out[i] = known_val[1]
        elseif k == length(known_idx)
            out[i] = known_val[end]
        else
            lo_i, hi_i = known_idx[k], known_idx[k+1]
            t = (i - lo_i) / (hi_i - lo_i)
            out[i] = known_val[k] + t * (known_val[k+1] - known_val[k])
        end
    end
    return out
end

# _build_heights/_build_depths/plant_area_index_profile: see utils.jl
# (generic grid + PAI construction, not site-specific).

# Both ends of every Ts and Sws depth range: the depths _build_depths lands nodes on
_observed_depths_m(hourly) = [d for prefix in ("Ts", "Sws") for r in _hourly_depth_series(hourly, prefix) for d in (r.top_m, r.bottom_m)]

# Soil profile from soil_source(site_name): SLGA at the site's cell (:slga, _slga_soil_profile
# in mapper.jl), measured horizons (:horizons), peat (:peat, tables/peat.toml), or a fixed
# CAMPBELL_NORMAN_TEXTURES texture. :slga_uniform is the SLGA profile collapsed to one
# depth-weighted-mean slab (flatten_soil_profile, utils.jl), so micropoint can run the same
# flat numbers.
#
# Cached per (site_name, depths, source): prepare_site/prepare_site_silo fetch the same site's
# same depth grid.
const _SOIL_CACHE = Dict{Tuple{String,Vector{Float64},Symbol},Any}()

function _fetch_soil_profile(site_name, latitude, longitude, depths)
    source = soil_source(site_name)
    cache_key = (site_name, ustrip.(u"m", depths), source)
    haskey(_SOIL_CACHE, cache_key) && return _SOIL_CACHE[cache_key]
    soil_profile = if source === :slga || source === :slga_uniform
        _slga_soil_profile(latitude, longitude, depths)
    elseif source === :horizons
        soil_profile_from_horizons(soil_horizons(site_name), depths, pedotransfer_model_choice)
    elseif source === :peat
        soil_profile_from_peat(PEAT, depths)
    else
        soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[source], depths)
    end
    if source === :slga_uniform
        soil_profile = flatten_soil_profile(soil_profile, depths)
    end
    _SOIL_CACHE[cache_key] = soil_profile
    return soil_profile
end

# Fills `primary`'s missing hours from `backup`, mapped onto the primary
# instrument by a least-squares line fitted on hours where both report.
function _fill_from_backup(primary, backup)
    both = .!ismissing.(primary) .& .!ismissing.(backup)
    x = Float64.(backup[both])
    y = Float64.(primary[both])
    slope = cov(x, y) / var(x)
    intercept = mean(y) - slope * mean(x)
    return [ismissing(p) && !ismissing(b) ? intercept + slope * b : p for (p, b) in zip(primary, backup)]
end

# Forcing columns with gaps filled from the site's forcing_backups (sites.toml), before
# interpolation, so the gap masks downstream mark only hours no tower instrument covers. Backup variables are summed (missing if any is).
function _fill_forcing_from_backups(hourly, forcing, forcing_backups)
    for (quantity, backup_names) in pairs(forcing_backups)
        name = getproperty(forcing, quantity)
        backup = sum(hourly[!, b] for b in backup_names)
        hourly[!, name] = _fill_from_backup(hourly[!, name], backup)
    end
    return hourly
end

# Gap-free hourly forcing from the columns forcing_selection.jl's forcing_sources chose
# (tower sensors, or reanalysis), interpolated across gaps (fill_gaps)
function _gapfilled_hourly_forcing(resolved, hourly)
    (; site_name, forcing, elevation, canopy_height) = resolved
    pressure = site_pressure(forcing.air_pressure === nothing ? fill(missing, nrow(hourly)) : hourly[!, forcing.air_pressure], elevation)
    pressure_kPa = ustrip.(u"kPa", pressure)
    ta_col = hourly[!, forcing.air_temperature]
    humidity_col = hourly[!, forcing.humidity]
    if forcing.temperature_is_sonic
        forcing.humidity_is_absolute || error("$site_name: a sonic temperature needs absolute humidity to convert")
        ta_col = [ismissing(t) || ismissing(ah) ? missing : air_temperature_from_sonic_temperature(t, specific_humidity(ah, t, p))
                  for (t, ah, p) in zip(ta_col, humidity_col, pressure_kPa)]
    end
    ws_col = wind_at_height.(hourly[!, forcing.wind_speed], forcing.wind_height, forcing.reference_height, canopy_height)
    reference_temperature = fill_gaps(ta_col) .* u"°C"

    # Relative humidity: measured, or from absolute humidity at the air temperature
    ta_C = ustrip.(u"°C", reference_temperature)
    rh_raw = forcing.humidity_is_absolute ?
        [ismissing(ah) ? missing : _rh_from_ah(ah, t) for (ah, t) in zip(humidity_col, ta_C)] : humidity_col
    reference_humidity = fill_gaps(rh_raw) ./ 100.0
    shortwave_col, longwave_col = hourly[!, forcing.downwelling_shortwave], hourly[!, forcing.downwelling_longwave]
    reference_wind_speed  = fill_gaps(ws_col) .* u"m/s"
    global_radiation      = max.(fill_gaps(shortwave_col), 0.0) .* u"W/m^2"
    longwave_radiation    = max.(fill_gaps(longwave_col), 0.0) .* u"W/m^2"
    precip = hourly[!, forcing.rainfall]
    rainfall_hourly       = coalesce.(precip, 0.0) .* u"kg/m^2"  # any remaining missing (partial-coverage days) assumed zero, not extrapolated

    # What actually went into the model (post gap-fill) + which hours were
    # originally missing -- report.jl's forcing-sanity plots color these separately.
    forcing_used = (;
        Ta = ustrip.(u"°C", reference_temperature), RH = reference_humidity .* 100.0,
        Ws = ustrip.(u"m/s", reference_wind_speed), Fsd = ustrip.(u"W/m^2", global_radiation),
        Fld = ustrip.(u"W/m^2", longwave_radiation), Precip = ustrip.(u"kg/m^2", rainfall_hourly),
    )
    filled_mask = (;
        Ta = ismissing.(ta_col), RH = ismissing.(rh_raw), Ws = ismissing.(ws_col),
        Fsd = ismissing.(shortwave_col), Fld = ismissing.(longwave_col), Precip = ismissing.(precip),
    )

    return (; reference_temperature, reference_humidity, reference_wind_speed, global_radiation,
              longwave_radiation, rainfall_hourly, pressure, forcing_used, filled_mask)
end

# Hourly tower measurements of a resolved site's years: the validation series and the forcing
# sensors forcing_sources picked, their gaps filled from the site's backup sensors, plus the
# reanalysis columns for whatever the tower can't supply
function read_hourly_site_data(resolved)
    (; site_name, years, path, forcing, forcing_backups, height_series) = resolved
    println("  Reading $path for $site_name $years...")
    print_forcing_sources(forcing)
    extra_vars = String[tower_forcing_names(forcing)..., Iterators.flatten(values(forcing_backups))...]
    tower = read_ozflux_nc(path, years, site_corrections(site_name); height_series, extra_vars)
    hourly = _fill_forcing_from_backups(_aggregate_hourly(tower), forcing, forcing_backups)
    add_reanalysis_forcing!(hourly, resolved, forcing.reanalysis)
    return hourly
end

# Reads + aggregates all requested years, builds the MicroProblem. Doesn't
# run anything -- that's run_site. `max_days` truncates the run for fast
# debug-loop iteration (unlike oznet, ozflux is keyed on whole-year files,
# not a continuous date range, so there's no sim_start/sim_end to shorten).
function prepare_site(site_name, years; max_days=nothing, litter_model=site_litter(site_name))
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, forcing, latitude, longitude, leaf_area_index, height_series) = resolved

    hourly = read_hourly_site_data(resolved)
    # fill_gaps would hold a leading/trailing forcing gap flat at the nearest observation for
    # days (sun-up Fsd all night, fixed Ta) -- start and end on whole observed days instead.
    hourly = _trim_unforced_days(hourly, forcing)
    ndays = nrow(hourly) ÷ 24
    max_days !== nothing && (ndays = min(ndays, max_days))
    hourly = hourly[1:(ndays * 24), :]
    n = ndays * 24
    println("  $ndays days ($n hours) from $(first(hourly.DateTime)) after 30-min -> hourly aggregation")

    (; reference_temperature, reference_humidity, reference_wind_speed, global_radiation,
       longwave_radiation, rainfall_hourly, pressure, forcing_used, filled_mask) =
        _gapfilled_hourly_forcing(resolved, hourly)
    mean_air_temperature = mean(ustrip.(u"°C", reference_temperature)) * u"°C"

    environment_hourly = HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity, reference_wind_speed,
        global_radiation, longwave_radiation,
        cloud_cover = fill(0.0, n),  # unused once global_radiation/longwave_radiation are both supplied directly
        rainfall = rainfall_hourly,
        zenith_angle = nothing,
    )
    daily_rainfall = [sum(ustrip.(u"kg/m^2", rainfall_hourly[(d-1)*24+1 : d*24])) for d in 1:ndays] .* u"kg/m^2"
    environment_daily = DailyTimeseries(;
        shade = fill(0.0, ndays),
        soil_wetness = fill(0.0, ndays),  # unused -- DynamicSoilMoisture simulates it instead
        surface_emissivity = fill(emissivity, ndays),
        cloud_emissivity = fill(emissivity, ndays),
        rainfall = daily_rainfall,
        deep_soil_temperature = fill(mean_air_temperature, ndays),  # no borehole climatology available; site-year mean air temp as a stand-in
        leaf_area_index = fill(leaf_area_index, ndays),  # only read by NoCanopy's ET partitioning; MultilayerCanopy uses its own plant_area_index
        albedo = fill(site_albedo(site_name), ndays),
        roughness_height = fill(site_ground_roughness_height(site_name), ndays),
    )

    site = Site(;
        latitude, longitude, resolved.elevation, slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), sky_view_fraction = 1.0,
    )

    depths = _build_depths(SIM_DEPTHS_M, _observed_depths_m(hourly))
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths)

    heights = _build_heights(canopy_height, reference_height, first.(height_series))

    canopy = canopy_model(site_name, resolved, heights)
    model = micro_model(site_name, site, soil_profile, depths, heights, canopy, rainfall_schedule_choice; litter_model)

    # Initial soil state from the site's own first observed Ts/Sws (falls
    # back to the mean air temperature / a mid-range moisture guess for any
    # depth with no observation to anchor to).
    initial_soil_temperature = initial_temperature_profile(_hourly_depth_series(hourly, "Ts"), depths, mean_air_temperature)
    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2, soil_profile; scale=sws_obs_scale(site_name, soil_profile))
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)  # bottom node always saturated (no drainage below the modelled soil column)

    inputs = MicroInputs(;
        environment_minmax = nothing, environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    problem = MicroProblem(model, inputs; days = days_of_year(hourly.DateTime),
        time_mode = ConsecutiveDayMode(; spinup_first_day = false))

    return (; site_name, years, resolved, hourly, t_model = hourly.DateTime, ndays, n, heights, depths, height_series,
              forcing_used, filled_mask, problem)
end

# Hourly rows from the first midnight at or after the first hour with every
# forcing variable observed, to the last 23:00 at or before the last such hour;
# pressure isn't needed (site_pressure).
function _trim_unforced_days(hourly, forcing)
    cols = [forcing.air_temperature, forcing.humidity, forcing.wind_speed,
            forcing.downwelling_shortwave, forcing.downwelling_longwave]
    forced = [all(c -> !ismissing(hourly[i, c]), cols) for i in 1:nrow(hourly)]
    i0, i1 = findfirst(forced), findlast(forced)
    i0 === nothing && error("no hour with every forcing variable observed")
    first_row = findfirst(i -> i >= i0 && Dates.hour(hourly.DateTime[i]) == 0, 1:nrow(hourly))
    last_row = findlast(i -> i <= i1 && Dates.hour(hourly.DateTime[i]) == 23, 1:nrow(hourly))
    return hourly[first_row:last_row, :]
end

# Day of year of each day in an hourly series starting at midnight -- the
# solar geometry's `days`, so a run that doesn't start on 1 January still
# sees the right sun.
function days_of_year(t)
    Dates.hour(first(t)) == 0 || error("hourly series must start at 00:00, starts at $(first(t))")
    return Dates.dayofyear.(t[1:24:end])
end

# MultilayerCanopy for woody sites (NoCanopy for grass, crop and pasture), from the config.jl
# model choices and the site's resolved traits. Shared by every prepare_site* so the forcing
# modes run the same canopy.
function canopy_model(site_name, resolved, heights)
    is_woody(site_name) || return NoCanopy()
    (; canopy_height, leaf_area_index, leaf_reflectance, leaf_transmittance) = resolved
    canopy_heights = heights[heights .<= canopy_height]
    plant_area_index = plant_area_index_profile(pai_shape(site_name), canopy_heights, canopy_height, leaf_area_index)
    return MultilayerCanopy(; canopy_height, plant_area_index,
        shortwave_model = TwoStreamRadiation(; leaf_reflectance, leaf_transmittance),
        longwave_model = longwave_model_choice, wind_model = canopy_wind_model_choice,
        air_profile_model = canopy_air_profile_model_choice, interception_model = interception_model_choice,
        leaf_parameters = leaf_parameters(site_name), stomatal_model = MoistureResponsiveStomatalConductance(),
        leaf_convection_model = leaf_convection_model_choice, convergence_model = canopy_convergence_model_choice,
        ground_coupling = ground_coupling_choice)
end

# MicroModel from the config.jl model choices; `rainfall_schedule` differs by
# forcing source (tower Precip is per-timestep, SILO a daily total).
function micro_model(site_name, site, soil_profile, depths, heights, canopy, rainfall_schedule;
    litter_model=site_litter(site_name))
    soil = SoilModel(; depths, profile = soil_profile,
        properties_model = soil_properties_model_choice, hydraulic_model = soil_hydraulic_model(site_name),
        moisture_strategy = soil_moisture_strategy_choice, energy_model = soil_energy_model_choice)
    atmosphere = AtmosphereModel(; heights, boundary_layer_model = boundary_layer_model_choice)
    config = MicroConfig(; convergence = convergence_choice, rainfall_schedule,
        canopy_soil_coupling = canopy_soil_coupling_choice)
    surface = SurfaceModel(; litter_model)
    return MicroModel(; hours = 0:1:23, site, soil, surface, atmosphere, canopy, config,
        vapour_pressure_equation = vapour_pressure_equation_choice)
end

# Depth series (matching discover_depth_series's (; top_m, bottom_m, values)) from the hourly
# table's `<prefix>_<depth_name>` columns, with the column name
function _hourly_depth_series(hourly, prefix)
    pattern = Regex("^" * prefix * raw"_([\d.]+)(?:-([\d.]+))?m$")
    result = NamedTuple{(:top_m, :bottom_m, :column, :values),Tuple{Float64,Float64,String,Vector{Union{Float64,Missing}}}}[]
    for name in names(hourly)
        m = match(pattern, name)
        m === nothing && continue
        top_m = parse(Float64, m[1])
        push!(result, (; top_m, bottom_m=m[2] === nothing ? top_m : parse(Float64, m[2]), column=name, values=hourly[!, name]))
    end
    return sort(result; by=r -> (r.top_m, r.bottom_m))
end

# First available value of each series in `depth_series`, at the middle of its depth range,
# sorted by depth.
function _known_depth_points(depth_series)
    known_d = Float64[]
    known_v = Float64[]
    for r in depth_series
        idx = findfirst(!ismissing, r.values)
        idx === nothing && continue
        push!(known_d, (r.top_m + r.bottom_m) / 2)
        push!(known_v, r.values[idx])
    end
    order = sortperm(known_d)
    return known_d[order], known_v[order]
end

# Linear interpolation at `dm`, flat-extrapolated beyond the shallowest/deepest point.
function _interp_flat(known_d, known_v, dm)
    dm <= known_d[1] && return known_v[1]
    dm >= known_d[end] && return known_v[end]
    j = searchsortedlast(known_d, dm)
    t = (dm - known_d[j]) / (known_d[j+1] - known_d[j])
    return known_v[j] + t * (known_v[j+1] - known_v[j])
end

# Kelvin, not °C: initial_soil_temperature feeds the soil ODE state directly,
# and °C (an affine unit) has no multiplicative identity -- the solver errors.
function initial_temperature_profile(ts_depths, depths, fallback_temperature)
    known_d, known_v_C = _known_depth_points(ts_depths)
    isempty(known_d) && return fill(u"K"(fallback_temperature), length(depths))
    return [u"K"(_interp_flat(known_d, known_v_C, ustrip(u"m", d)) * u"°C") for d in depths]
end

# Clamped to each depth's saturation, 1 - bulk density / solid density
function initial_moisture_profile(sws_depths, depths, fallback_fraction, soil_profile; scale=1.0)
    saturation = 1 .- soil_profile.bulk_density ./ soil_profile.mineral_density
    known_d, known_v = _known_depth_points(sws_depths)
    isempty(known_d) && return min.(fallback_fraction, saturation)
    known_v = known_v .* scale  # before clamping -- see sws_obs_scale (utils.jl)
    return [clamp(_interp_flat(known_d, known_v, ustrip(u"m", d)), 0.02, s) for (d, s) in zip(depths, saturation)]
end

function announce_and_prepare_site(site_name, years; max_days=nothing, litter_model=site_litter(site_name))
    println("\n" * "="^72)
    println("Site $site_name  years=$years  litter=$(nameof(typeof(litter_model)))")
    println("="^72)
    print_corrections(site_name)
    return prepare_site(site_name, years; max_days, litter_model)
end

# Thread-safe: prepare_site does the file and network reads
function solve_site(prep)
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  $(prep.site_name): $(prep.n) hours solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output, solve_time)
end

run_site(site_name, years; kw...) = solve_site(announce_and_prepare_site(site_name, years; kw...))
