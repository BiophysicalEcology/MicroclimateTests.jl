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
    return clamp(ustrip(u"Pa", p_v) / ustrip(u"Pa", p_sat) * 100.0, 0.0, 100.0)
end

_site_year_paths(site_name, years) =
    [joinpath(OZFLUX_DATA_DIR, "$(site_name)_$(y)_L3.nc") for y in years]

function resolve_site(site_name, years)
    haskey(SITE_LEAF_AREA_INDEX, site_name) ||
        error("No leaf_area_index configured for site \"$site_name\" -- add it to SITE_LEAF_AREA_INDEX in config.jl")
    paths = _site_year_paths(site_name, years)
    all(isfile, paths) || error("Missing data file(s): " * join(filter(!isfile, paths), ", "))

    global_attrib = Dict(NCDatasets.NCDataset(paths[1]).attrib)  # materialized -- not tied to the (about to close) dataset handle
    raw_canopy_height = _parse_length_attr(global_attrib["canopy_height"])
    tower_height  = _parse_length_attr(global_attrib["tower_height"])
    latitude  = parse(Float64, string(global_attrib["latitude"])) * u"°"
    longitude = parse(Float64, string(global_attrib["longitude"])) * u"°"

    # reference_height: the height the forcing Ta/RH/Ws are actually valid at
    # (see config.jl's SITE_REFERENCE_HEIGHT_M comment) -- NOT always
    # tower_height (Wallaby's only wind/temp sensor sits at 5m, well below
    # its 12m tower). forcing_vars: which tower variable supplies each,
    # for sites where the bare merged composite's real height doesn't match.
    reference_height = get(SITE_REFERENCE_HEIGHT_M, site_name, ustrip(u"m", tower_height)) * u"m"
    forcing_vars = get(SITE_FORCING_VARS, site_name, (ta="Ta", ws="Ws", ah=nothing))

    # canopy_height can never legitimately exceed reference_height --
    # Microclimate.jl always sets reference_height=last(heights), and
    # _build_heights forces canopy_height onto that same grid as an interior
    # point, so a taller canopy_height would silently push reference_height
    # off the top of the grid. Clamp + warn rather than let that happen
    # quietly (Wallaby: file says canopy_height=8-10m, but its only wind/temp
    # sensor genuinely sits at just 5m -- no above-canopy forcing exists for
    # this stand, so the simulated canopy is capped at what the forcing can
    # actually support).
    canopy_height = min(raw_canopy_height, reference_height)
    if canopy_height < raw_canopy_height
        @warn "canopy_height clamped to reference_height for $site_name" raw_canopy_height reference_height
    end

    return (; site_name, years, paths, canopy_height, tower_height, reference_height, forcing_vars,
              latitude, longitude, global_attrib,
              leaf_area_index = SITE_LEAF_AREA_INDEX[site_name], leaf_reflectance = SITE_LEAF_REFLECTANCE[site_name],
              leaf_transmittance = SITE_LEAF_TRANSMITTANCE[site_name],
              height_series = get(SITE_HEIGHT_SERIES, site_name, Tuple{Float64,String}[]))
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
# grouped, which is correct regardless of starting phase. This must land on
# the same on-the-hour convention SILO's t_model uses, or every
# DateTime-based alignment/gap-fill match against a SILO donor silently fails.
_hour_group_start(t) = Dates.minute(t) == 30 ? floor(t, Dates.Hour) : floor(t, Dates.Hour) - Dates.Hour(1)

function _aggregate_hourly(df)
    n_per_hour = 60 ÷ OZFLUX_NATIVE_STEP_MINUTES
    groups = Dict{DateTime,Vector{Int}}()
    for (i, t) in enumerate(df.DateTime)
        push!(get!(groups, _hour_group_start(t), Int[]), i)
    end
    # Drop any partial hour at the edges (fewer than n_per_hour sub-steps) --
    # an incomplete average would silently bias that hour rather than
    # reflecting a real QC gap; only affects the very first/last hour of the
    # requested span, or a genuine gap in the raw time index.
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
                (ismissing(v) || (has_qc && df[r, qc_name] != OZFLUX_GOOD_QC)) && continue
                push!(vals, v)
            end
            isempty(vals) || (hourly[k] = mean(vals))
        end
        out[!, name] = hourly
    end
    return out
end

# Fills `missing` entries in `tower_vals` (at `tower_t` DateTimes) from
# `donor_vals` (at `donor_t` DateTimes) where a matching timestamp exists;
# leaves them `missing` otherwise (fill_gaps handles the remainder).
function _gap_fill_from_donor(tower_t, tower_vals, donor_t, donor_vals)
    idx = Dict(t => i for (i, t) in enumerate(donor_t))
    out = Vector{Union{Float64,Missing}}(missing, length(tower_vals))
    for i in eachindex(tower_vals)
        if !ismissing(tower_vals[i])
            out[i] = tower_vals[i]
        elseif haskey(idx, tower_t[i])
            out[i] = donor_vals[idx[tower_t[i]]]
        end
    end
    return out
end

# SILO has no hourly rain, only a daily total -- unlike the other donor
# variables this can't fill individual missing hours. Only applied to tower
# days missing ALL 24 hours (a partially-observed day's known hours would
# make the true within-day distribution ambiguous either way); the whole
# day's SILO total is dumped at midnight, matching how daily totals are
# conventionally applied (e.g. oznet's DailyRainfall schedule).
function _gap_fill_precip_daily(tower_t, tower_precip, donor_dates, donor_daily_mm)
    out = collect(Union{Float64,Missing}, tower_precip)
    donor_idx = Dict(d => i for (i, d) in enumerate(donor_dates))
    for day_start in 1:24:(length(out) - 23)
        day_range = day_start:(day_start + 23)
        all(ismissing, @view out[day_range]) || continue
        d = Date(tower_t[day_start])
        haskey(donor_idx, d) || continue
        out[day_range] .= 0.0
        out[day_start] = donor_daily_mm[donor_idx[d]]
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

# Depth-m values from _hourly_depth_series's (depth_m, values) pairs, for
# both Ts and Sws -- the set of depths _build_depths should land nodes on.
_observed_depths_m(hourly) = vcat([r.depth_m for r in _hourly_depth_series(hourly, "Ts")],
                                   [r.depth_m for r in _hourly_depth_series(hourly, "Sws")])

# Live SLGA fetch + pedotransfer (same as oznet/pipeline.jl) -- SLGA has no
# point-query API, so a uniform profile is derived from a small bounding box
# around the site coordinates. Organic litter/crust cap on the top ~2.5 cm
# (nodes 1:3) matches oznet's/NicheMapR's default cap=1 behaviour. Prints a
# sanity check against the site's own metadata BulkDensity attribute (present
# on at least some sites, e.g. Whroo) when available -- SLGA vs a real
# point measurement, not a substitute for it (SLGA is a 90 m grid mean).
#
# soil_source can also be :slga_uniform -- the real per-depth SLGA profile
# (with the same organic_cap/BulkDensity corrections as :slga), collapsed to
# one depth-weighted-mean slab via flatten_soil_profile (utils.jl). Real
# per-depth SLGA data previously made micropoint's RunModelFull hang
# indefinitely (see micropoint/ozflux/README), and even flattening it
# per-model let the two models silently average it differently --
# :slga_uniform lets both Microclimate.jl and micropoint run on the exact
# same flat numbers for a genuine like-for-like comparison, while plain
# :slga (full depth resolution) stays available for Microclimate.jl-only runs.
#
# Cached per (site_name, depths, source) -- prepare_site/prepare_site_silo
# both fetch the same site's same depth grid, no need to hit SLGA twice; the
# source is part of the key since :slga and :slga_uniform must not share a
# cached result for the same (site_name, depths).
const _SOIL_CACHE = Dict{Tuple{String,Vector{Float64},Symbol},Any}()

function _fetch_soil_profile(site_name, latitude, longitude, depths, global_attrib)
    source = soil_source(site_name)
    cache_key = (site_name, ustrip.(u"m", depths), source)
    haskey(_SOIL_CACHE, cache_key) && return _SOIL_CACHE[cache_key]
    if source !== :slga && source !== :slga_uniform
        soil_profile = soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[source], depths)
        if organic_cap(site_name)
            soil_profile.mineral_conductivity[1:3] .= 0.05u"W/m/K"
            soil_profile.mineral_heat_capacity[1:3] .= 1920.0u"J/kg/K"
        end
        _SOIL_CACHE[cache_key] = soil_profile
        return soil_profile
    end
    b = soil_area_buffer_deg
    lon_dd, lat_dd = ustrip(u"°", longitude), ustrip(u"°", latitude)
    area = Extent(X = (lon_dd - b, lon_dd + b), Y = (lat_dd - b, lat_dd + b))
    soil_result = build_soil_profile(SLGA, area; depths, pedotransfer_model = pedotransfer_model_choice)
    soil_profile = soil_result.soil_profile
    if organic_cap(site_name)
        soil_profile.mineral_conductivity[1:3] .= 0.05u"W/m/K"
        soil_profile.mineral_heat_capacity[1:3] .= 1920.0u"J/kg/K"
    end
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
    if source === :slga_uniform
        soil_profile = flatten_soil_profile(soil_profile, depths)
    end
    _SOIL_CACHE[cache_key] = soil_profile
    return soil_profile
end

# PAI -> shade fraction via Beer-Lambert (Campbell & Norman), used only by
# canopy_mode=:legacy -- MultilayerCanopy (canopy_mode=:full) computes its
# own radiative transfer and doesn't use `shade` at all.
_legacy_shade(leaf_area_index, extinction_coefficient) = 1.0 - exp(-extinction_coefficient * leaf_area_index)

# Reads + aggregates all requested years, builds the MicroProblem. Doesn't
# run anything -- that's run_site. `max_days` truncates the run for fast
# debug-loop iteration (unlike oznet, ozflux is keyed on whole-year files,
# not a continuous date range, so there's no sim_start/sim_end to shorten).
# `gap_fill_donor` (see silo_gapfill_donor): substitutes real SILO-forced
# hourly values for missing tower Ta/RH/Ws/Fsd before falling back to plain
# interpolation (fill_gaps) for whatever's still missing.
# `canopy_mode`: :full (MultilayerCanopy) or :legacy (NoCanopy + a PAI-derived
# shade fraction + a wind knockdown + a large horizon angle -- how forest
# sites were approximated before MultilayerCanopy existed; run both to see
# what the full canopy model buys over that approximation).
function prepare_site(site_name, years; max_days=nothing, gap_fill_donor=nothing, canopy_mode=:full)
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, forcing_vars, latitude, longitude, leaf_area_index,
        leaf_reflectance, leaf_transmittance, height_series) = resolved

    println("  Reading ", length(resolved.paths), " file(s) for $site_name $years...")
    raw_parts = [read_ozflux_nc(path; height_series) for path in resolved.paths]
    raw = sort(vcat(raw_parts...; cols=:union), :DateTime)
    hourly = _aggregate_hourly(raw)
    ndays = nrow(hourly) ÷ 24
    max_days !== nothing && (ndays = min(ndays, max_days))
    hourly = hourly[1:(ndays * 24), :]
    n = ndays * 24
    println("  $ndays days ($n hours) after 30-min -> hourly aggregation")

    donor(col, name) = gap_fill_donor === nothing ? col :
        _gap_fill_from_donor(hourly.DateTime, col, gap_fill_donor.t_model, getproperty(gap_fill_donor, name))

    # ── Forcing (gap-filled; must be complete for the solver) ────────────────
    # forcing_vars.ta/.ws: the tower variable actually at reference_height
    # (bare "Ta"/"Ws" for most sites; a dedicated tower-top sensor where the
    # bare composite's real height doesn't match -- see config.jl).
    ta_col = hourly[!, forcing_vars.ta]
    ws_col = hourly[!, forcing_vars.ws]
    reference_temperature = fill_gaps(donor(ta_col, :reference_temperature_C)) .* u"°C"

    # RH: sites with a forcing_vars.ah override have no bare "RH" of verified
    # height either, so always derive from AH rather than trusting it. Other
    # sites keep bare RH (gaps preferentially filled from the tower's own Ah,
    # different instrument, same site/hour) before falling back to the SILO
    # donor -- also a direct sanity check on the RH scale fix (utils.jl's
    # _fix_rh_scale).
    ta_C = ustrip.(u"°C", reference_temperature)
    rh_raw = if forcing_vars.ah !== nothing
        ah_col = hourly[!, forcing_vars.ah]
        [ismissing(ah) ? missing : _rh_from_ah(ah, t) for (ah, t) in zip(ah_col, ta_C)]
    else
        rr = hourly.RH
        if "Ah" in names(hourly)
            ah_derived = [ismissing(ah) ? missing : _rh_from_ah(ah, t) for (ah, t) in zip(hourly.Ah, ta_C)]
            both = [(Float64(rh), ad) for (rh, ad) in zip(hourly.RH, ah_derived) if !ismissing(rh) && !ismissing(ad)]
            if length(both) > 1
                obs = first.(both); der = last.(both)
                @printf("  RH cross-check vs Ah-derived: r=%+.3f  bias=%+.2f%%  (n=%d)\n", cor(obs, der), mean(der .- obs), length(both))
            end
            rr = [ismissing(rh) ? ad : rh for (rh, ad) in zip(hourly.RH, ah_derived)]
        end
        rr
    end
    reference_humidity = clamp.(fill_gaps(donor(rh_raw, :reference_humidity_pct)) ./ 100.0, 0.0, 1.0)
    lp = legacy_params(site_name)
    wind_multiplier = canopy_mode == :legacy ? lp.wind_multiplier : 1.0
    reference_wind_speed  = max.(fill_gaps(donor(ws_col, :reference_wind_speed_ms)) .* wind_multiplier, 0.1) .* u"m/s"
    global_radiation      = max.(fill_gaps(donor(hourly.Fsd, :global_radiation_Wm2)), 0.0) .* u"W/m^2"
    longwave_radiation    = max.(fill_gaps(donor(hourly.Fld, :longwave_radiation_Wm2)), 0.0) .* u"W/m^2"  # donor backs Fld out of the SILO run's sky_temperature (see silo_gapfill_donor)
    pressure              = fill_gaps(hourly.ps) .* u"kPa" .|> p -> uconvert(u"Pa", p)
    precip = gap_fill_donor === nothing || !haskey(gap_fill_donor, :daily_rainfall_mm) ? hourly.Precip :
        _gap_fill_precip_daily(hourly.DateTime, hourly.Precip, gap_fill_donor.daily_dates, gap_fill_donor.daily_rainfall_mm)
    rainfall_hourly       = coalesce.(precip, 0.0) .* u"kg/m^2"  # any remaining missing (partial-coverage days) assumed zero, not extrapolated

    mean_air_temperature = mean(ustrip.(u"°C", reference_temperature)) * u"°C"

    # What actually went into the model (post gap-fill) + which hours were
    # originally missing from the tower file (donor- or interpolation-filled,
    # either way) -- report.jl's forcing-sanity plots color these separately.
    forcing_used = (;
        Ta = ustrip.(u"°C", reference_temperature), RH = reference_humidity .* 100.0,
        Ws = ustrip.(u"m/s", reference_wind_speed), Fsd = ustrip.(u"W/m^2", global_radiation),
        Fld = ustrip.(u"W/m^2", longwave_radiation), Precip = ustrip.(u"kg/m^2", rainfall_hourly),
    )
    filled_mask = (;
        Ta = ismissing.(ta_col), RH = ismissing.(rh_raw), Ws = ismissing.(ws_col),
        Fsd = ismissing.(hourly.Fsd), Fld = ismissing.(hourly.Fld), Precip = ismissing.(precip),
    )

    environment_hourly = HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity, reference_wind_speed,
        global_radiation, longwave_radiation,
        cloud_cover = fill(0.5, n),  # unused once global_radiation/longwave_radiation are both supplied directly
        rainfall = rainfall_hourly,
        zenith_angle = nothing,
    )
    daily_rainfall = [sum(ustrip.(u"kg/m^2", rainfall_hourly[(d-1)*24+1 : d*24])) for d in 1:ndays] .* u"kg/m^2"
    shade_level = canopy_mode == :legacy ? _legacy_shade(leaf_area_index, lp.extinction_coefficient) : 0.0
    environment_daily = DailyTimeseries(;
        shade = fill(shade_level, ndays),
        soil_wetness = fill(0.0, ndays),  # unused -- DynamicSoilMoisture simulates it instead
        surface_emissivity = fill(emissivity, ndays),
        cloud_emissivity = fill(emissivity, ndays),
        rainfall = daily_rainfall,
        deep_soil_temperature = fill(mean_air_temperature, ndays),  # no borehole climatology available; site-year mean air temp as a stand-in
        leaf_area_index = fill(leaf_area_index, ndays),  # only read by NoCanopy's ET partitioning; MultilayerCanopy uses its own plant_area_index
    )

    horizon = canopy_mode == :legacy ? fill(lp.horizon_angle, 24) : fill(0.0u"°", 24)
    roughness_height = canopy_mode == :legacy ? lp.roughness_height : 0.004u"m"
    site = Site(;
        latitude, longitude, elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = horizon, sky_view_fraction = 1.0,
        albedo = site_albedo(site_name), roughness_height,
        atmospheric_pressure = mean(pressure),
    )

    depths = _build_depths(SIM_DEPTHS_M, _observed_depths_m(hourly))
    println("  Fetching SLGA soil texture...")
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths, resolved.global_attrib)

    extra_heights_m = [h for (h, _) in height_series]
    heights = _build_heights(canopy_height, reference_height, extra_heights_m)

    canopy_model = if canopy_mode == :legacy
        NoCanopy()
    else
        pai_shape = get(SITE_PAI_SHAPE, site_name, :uniform)
        canopy_heights = heights[heights .<= canopy_height]
        plant_area_index = plant_area_index_profile(pai_shape, canopy_heights, canopy_height, leaf_area_index)
        MultilayerCanopy(; canopy_height, plant_area_index,
            shortwave_model=TwoStreamRadiation(leaf_reflectance, leaf_transmittance),
            leaf_convection_model = leaf_convection_model_choice, interception_model = interception_model_choice,
            leaf_parameters = leaf_parameters(site_name), longwave_model = longwave_model_choice,
            stomatal_model = MoistureResponsiveStomatalConductance(),
            convergence_model = canopy_convergence_model_choice, air_profile_model = canopy_air_profile_model_choice)
    end

    config = MicroConfig(;
        convergence = convergence_choice,
        rainfall_schedule = rainfall_schedule_choice,
        soil_moisture_strategy = soil_moisture_strategy_choice,
        canopy_soil_convergence = canopy_soil_convergence_choice,
    )
    model = MicroModel(;
        hours = 0:1:23, depths, heights,
        soil_properties_model = soil_properties_model_choice,
        soil_hydraulic_model = soil_hydraulic_model(site_name),
        canopy_model, config,
    )

    # Initial soil state from the site's own first observed Ts/Sws (falls
    # back to the mean air temperature / a mid-range moisture guess for any
    # depth with no observation to anchor to).
    initial_soil_temperature = initial_temperature_profile(_hourly_depth_series(hourly, "Ts"), depths, mean_air_temperature)
    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2)
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)  # bottom node always saturated (no drainage below the modelled soil column)

    inputs = MicroInputs(;
        site, soil_profile, environment_minmax = nothing, environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    problem = MicroProblem(model, inputs; days = collect(1:ndays), time_mode = ConsecutiveDayMode(; spinup_first_day = false))

    return (; site_name, years, resolved, hourly, t_model = hourly.DateTime, ndays, n, heights, depths, height_series,
              canopy_mode, forcing_used, filled_mask, problem)
end

# Depth-series NamedTuples (matching discover_depth_series's format) from the
# hourly-aggregated frame's renamed `<prefix>_<depth>m` columns.
function _hourly_depth_series(hourly, prefix)
    pattern = Regex("^" * prefix * raw"_([\d.]+)m$")
    result = NamedTuple{(:depth_m, :values),Tuple{Float64,Vector{Union{Float64,Missing}}}}[]
    for name in names(hourly)
        m = match(pattern, name)
        m === nothing && continue
        push!(result, (depth_m=parse(Float64, m.captures[1]), values=hourly[!, name]))
    end
    return sort(result; by=r -> r.depth_m)
end

# First available (depth_m, value) pair from each series in `depth_series`,
# sorted by depth.
function _known_depth_points(depth_series)
    known_d = Float64[]
    known_v = Float64[]
    for (depth_m, values) in depth_series
        idx = findfirst(!ismissing, values)
        idx === nothing && continue
        push!(known_d, depth_m)
        push!(known_v, values[idx])
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

function initial_moisture_profile(sws_depths, depths, fallback_fraction)
    known_d, known_v = _known_depth_points(sws_depths)
    isempty(known_d) && return fill(fallback_fraction, length(depths))
    return [clamp(_interp_flat(known_d, known_v, ustrip(u"m", d)), 0.02, 0.5) for d in depths]
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
function prepare_site_silo(site_name, years; max_days=nothing)
    resolved = resolve_site(site_name, years)
    (; canopy_height, reference_height, latitude, longitude, leaf_area_index,
        leaf_reflectance, leaf_transmittance, height_series) = resolved
    sim_start, sim_end = _year_span(years)

    # Tower obs read first -- needed for observed depths (_build_depths) and
    # initial soil state before the SILO fetch/soil-texture fetch below.
    raw_parts = [read_ozflux_nc(path; height_series) for path in resolved.paths]
    raw = sort(vcat(raw_parts...; cols=:union), :DateTime)
    hourly = _aggregate_hourly(raw)
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
    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2)

    site = Site(;
        latitude, longitude, elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), sky_view_fraction = 1.0,
        albedo = site_albedo(site_name), roughness_height = 0.004u"m",
        atmospheric_pressure = fetched.site.atmospheric_pressure,
    )
    println("  Fetching SLGA soil texture...")
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths, resolved.global_attrib)
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)  # bottom node always saturated (no drainage below the modelled soil column)

    config = MicroConfig(; convergence = convergence_choice, rainfall_schedule = DailyRainfall(),  # SILO rainfall is a daily total, not per-timestep
        soil_moisture_strategy = soil_moisture_strategy_choice, canopy_soil_convergence = canopy_soil_convergence_choice)
    model = MicroModel(; hours = 0:1:23, depths, heights,
        soil_properties_model = soil_properties_model_choice, soil_hydraulic_model = soil_hydraulic_model(site_name),
        #canopy_model,
        config)

    inputs = MicroInputs(;
        site, soil_profile,
        environment_minmax = fetched.environment_minmax, environment_daily = fetched.environment_daily,
        environment_hourly = fetched.environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    problem = MicroProblem(model, inputs; days = collect(1:ndays), time_mode = ConsecutiveDayMode(; spinup_first_day = false))

    # Model spans the full sim_start:sim_end calendar range (SILO forcing);
    # `hourly` (tower obs, for initial conditions + comparison) may cover a
    # different/shorter span -- report.jl aligns the two by DateTime, not position.
    t_model = [DateTime(sim_start) + Day(d) + Hour(h) for d in 0:(ndays - 1) for h in 0:23]
    # canopy_mode=:legacy (see prepare_site) isn't wired up for the SILO path.
    # sim_start/environment_daily kept (full sim_start:sim_end span, not
    # truncated by max_days) for silo_gapfill_donor's daily rainfall donor.
    return (; site_name, years, resolved, hourly, t_model, ndays, n, heights, depths, height_series, canopy_mode=:full,
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
    silo = run_site_silo(site_name, years; max_days)
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

function run_site(site_name, years; max_days=nothing, canopy_mode=:full)
    println("\n" * "="^72)
    println("Site $site_name  years=$years  canopy_mode=$canopy_mode")
    println("="^72)
    prep = prepare_site(site_name, years; max_days, canopy_mode)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output)
end
