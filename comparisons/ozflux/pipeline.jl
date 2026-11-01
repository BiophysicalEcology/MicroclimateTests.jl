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
    haskey(SITE_LEAF_AREA_INDEX, site_name) ||
        error("No leaf_area_index configured for site \"$site_name\" -- add it to SITE_LEAF_AREA_INDEX in config.jl")
    path = ozflux_site_path(OZFLUX_DATA_DIR, site_name)
    available = ozflux_years(path)
    missing_years = setdiff(years, available)
    isempty(missing_years) || error("$path has no data for year(s) $missing_years (has $available)")

    global_attrib = NCDatasets.NCDataset(ds -> Dict(ds.attrib), path)  # materialized -- not tied to the dataset handle
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
    forcing_backups = get(SITE_FORCING_BACKUPS, site_name, (;))

    # canopy_height can never legitimately exceed reference_height --
    # Microclimate.jl always sets reference_height=last(heights), and
    # _build_heights forces canopy_height onto that same grid as an interior
    # point, so a taller canopy_height would silently push reference_height
    # off the top of the grid. Clamp + warn rather than let that happen
    # quietly (Wallaby: file says canopy_height=8-10m, but its only wind/temp
    # sensor sits at just 5m -- no above-canopy forcing exists for
    # this stand, so the simulated canopy is capped at what the forcing can
    # actually support).
    canopy_height = min(raw_canopy_height, reference_height)
    if canopy_height < raw_canopy_height
        @warn "canopy_height clamped to reference_height for $site_name" raw_canopy_height reference_height
    end

    return (; site_name, years, path, canopy_height, tower_height, reference_height, forcing_vars, forcing_backups,
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

# Soil profile from soil_source(site_name): a fixed CAMPBELL_NORMAN_TEXTURES
# texture, or a live SLGA fetch (:slga, _slga_soil_profile in mapper.jl,
# needs MicroclimateMapper). Organic litter/crust cap on the top ~2.5 cm
# (nodes 1:3) matches oznet's/NicheMapR's default cap=1 behaviour.
#
# soil_source can also be :slga_uniform -- the real per-depth SLGA profile
# (with the same organic_cap/BulkDensity corrections as :slga), collapsed to
# one depth-weighted-mean slab via flatten_soil_profile (utils.jl). Real
# per-depth SLGA data previously made micropoint's RunModelFull hang
# indefinitely (see micropoint/ozflux/README), and even flattening it
# per-model let the two models silently average it differently --
# :slga_uniform lets both Microclimate.jl and micropoint run on the exact
# same flat numbers for a like-for-like comparison, while plain
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
    soil_profile = if source === :slga || source === :slga_uniform
        _slga_soil_profile(latitude, longitude, depths, global_attrib)  # mapper.jl (needs MicroclimateMapper)
    else
        soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[source], depths)
    end
    if organic_cap(site_name)
        soil_profile.mineral_conductivity[1:3] .= 0.05u"W/m/K"
        soil_profile.mineral_heat_capacity[1:3] .= 1920.0u"J/kg/K"
    end
    macropores(site_name) && apply_macropore_boost!(soil_profile, depths;
        max_boost=MACROPORE_MAX_BOOST, taper_depth=MACROPORE_TAPER_DEPTH)
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

# Forcing columns with gaps filled from config.jl's SITE_FORCING_BACKUPS, before
# any donor or interpolation, so the gap masks downstream mark only hours no
# tower instrument covers. Backup variables are summed (missing if any is).
function _fill_forcing_from_backups(hourly, forcing_vars, forcing_backups)
    primary_name = (; ta=forcing_vars.ta, ah=forcing_vars.ah, fsd="Fsd")
    for (key, backup_names) in pairs(forcing_backups)
        name = primary_name[key]
        backup = sum(hourly[!, b] for b in backup_names)
        hourly[!, name] = _fill_from_backup(hourly[!, name], backup)
    end
    return hourly
end

# Real, gap-filled hourly forcing: tower obs, SILO-donor-filled (see
# silo_gapfill_donor) wherever the tower itself is missing, then fill_gaps
# for any remainder. Shared by prepare_site (consumed hourly) and
# prepare_site_daily (degraded to a daily min/max/mean/total envelope
# afterward) so both start from the exact same real+gap-filled values --
# forcing_vars.ta/.ws: the tower variable actually at reference_height (bare
# "Ta"/"Ws" for most sites; a dedicated tower-top sensor where the bare
# composite's real height doesn't match -- see config.jl).
function _gapfilled_hourly_forcing(site_name, hourly, forcing_vars, gap_fill_donor, canopy_mode)
    donor(col, name) = gap_fill_donor === nothing ? col :
        _gap_fill_from_donor(hourly.DateTime, col, gap_fill_donor.t_model, getproperty(gap_fill_donor, name))

    ta_col = hourly[!, forcing_vars.ta]
    #ta_col = fill(missing, nrow(hourly))
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
                #@printf("  RH cross-check vs Ah-derived: r=%+.3f  bias=%+.2f%%  (n=%d)\n", cor(obs, der), mean(der .- obs), length(both))
            end
            rr = [ismissing(rh) ? ad : rh for (rh, ad) in zip(hourly.RH, ah_derived)]
        end
        rr
    end
    reference_humidity = fill_gaps(donor(rh_raw, :reference_humidity_pct)) ./ 100.0
    lp = legacy_params(site_name)
    wind_multiplier = canopy_mode == :legacy ? lp.wind_multiplier : 1.0
    reference_wind_speed  = fill_gaps(donor(ws_col, :reference_wind_speed_ms)) .* wind_multiplier .* u"m/s"
    global_radiation      = max.(fill_gaps(donor(hourly.Fsd, :global_radiation_Wm2)), 0.0) .* u"W/m^2"
    longwave_radiation    = max.(fill_gaps(donor(hourly.Fld, :longwave_radiation_Wm2)), 0.0) .* u"W/m^2"  # donor backs Fld out of the SILO run's sky_temperature (see silo_gapfill_donor)
    precip = gap_fill_donor === nothing || !haskey(gap_fill_donor, :daily_rainfall_mm) ? hourly.Precip :
        _gap_fill_precip_daily(hourly.DateTime, hourly.Precip, gap_fill_donor.daily_dates, gap_fill_donor.daily_rainfall_mm)
    rainfall_hourly       = coalesce.(precip, 0.0) .* u"kg/m^2"  # any remaining missing (partial-coverage days) assumed zero, not extrapolated

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

    return (; reference_temperature, reference_humidity, reference_wind_speed, global_radiation,
              longwave_radiation, rainfall_hourly, forcing_used, filled_mask)
end

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
    (; canopy_height, reference_height, forcing_vars, latitude, longitude, leaf_area_index, height_series) = resolved

    println("  Reading $(resolved.path) for $site_name $years...")
    extra_vars = collect(String, Iterators.flatten(values(resolved.forcing_backups)))
    hourly = _fill_forcing_from_backups(_aggregate_hourly(read_ozflux_nc(resolved.path, years; height_series, extra_vars)),
        forcing_vars, resolved.forcing_backups)
    # Without a donor, fill_gaps would hold a leading/trailing forcing gap
    # flat at the nearest observation for days (sun-up Fsd all night, fixed
    # Ta) -- start and end on whole observed days instead.
    gap_fill_donor === nothing && (hourly = _trim_unforced_days(hourly, forcing_vars))
    ndays = nrow(hourly) ÷ 24
    max_days !== nothing && (ndays = min(ndays, max_days))
    hourly = hourly[1:(ndays * 24), :]
    n = ndays * 24
    println("  $ndays days ($n hours) from $(first(hourly.DateTime)) after 30-min -> hourly aggregation")

    (; reference_temperature, reference_humidity, reference_wind_speed, global_radiation,
       longwave_radiation, rainfall_hourly, forcing_used, filled_mask) =
        _gapfilled_hourly_forcing(site_name, hourly, forcing_vars, gap_fill_donor, canopy_mode)
    pressure = fill_gaps(hourly.ps) .* u"kPa" .|> p -> uconvert(u"Pa", p)
    mean_air_temperature = mean(ustrip.(u"°C", reference_temperature)) * u"°C"
    lp = legacy_params(site_name)

    environment_hourly = HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity, reference_wind_speed,
        global_radiation, longwave_radiation,
        cloud_cover = fill(0.0, n),  # unused once global_radiation/longwave_radiation are both supplied directly
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
    soil_profile = _fetch_soil_profile(site_name, latitude, longitude, depths, resolved.global_attrib)

    extra_heights_m = [h for (h, _) in height_series]
    heights = _build_heights(canopy_height, reference_height, extra_heights_m)

    canopy = canopy_model(site_name, resolved, heights, canopy_mode)
    model = micro_model(site_name, site, soil_profile, depths, heights, canopy, rainfall_schedule_choice)

    # Initial soil state from the site's own first observed Ts/Sws (falls
    # back to the mean air temperature / a mid-range moisture guess for any
    # depth with no observation to anchor to).
    initial_soil_temperature = initial_temperature_profile(_hourly_depth_series(hourly, "Ts"), depths, mean_air_temperature)
    initial_soil_moisture = initial_moisture_profile(_hourly_depth_series(hourly, "Sws"), depths, 0.2; scale=sws_obs_scale(site_name, soil_profile))
    initial_soil_moisture[length(initial_soil_moisture)] = 1.0 - last(soil_profile.bulk_density) / last(soil_profile.mineral_density)  # bottom node always saturated (no drainage below the modelled soil column)

    inputs = MicroInputs(;
        environment_minmax = nothing, environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    problem = MicroProblem(model, inputs; days = days_of_year(hourly.DateTime),
        time_mode = ConsecutiveDayMode(; spinup_first_day = false))

    return (; site_name, years, resolved, hourly, t_model = hourly.DateTime, ndays, n, heights, depths, height_series,
              canopy_mode, forcing_used, filled_mask, problem)
end

# Hourly rows from the first midnight at or after the first hour with every
# forcing variable observed, to the last 23:00 at or before the last such hour.
function _trim_unforced_days(hourly, forcing_vars)
    cols = [forcing_vars.ta, forcing_vars.ws, something(forcing_vars.ah, "RH"), "Fsd", "Fld", "ps"]
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

# MultilayerCanopy for canopy_mode=:full (NoCanopy for :legacy), from the
# config.jl model choices and the site's resolved traits. Shared by every
# prepare_site* so the forcing modes run the same canopy.
function canopy_model(site_name, resolved, heights, canopy_mode)
    canopy_mode == :legacy && return NoCanopy()
    (; canopy_height, leaf_area_index, leaf_reflectance, leaf_transmittance) = resolved
    pai_shape = get(SITE_PAI_SHAPE, site_name, :uniform)
    canopy_heights = heights[heights .<= canopy_height]
    plant_area_index = plant_area_index_profile(pai_shape, canopy_heights, canopy_height, leaf_area_index)
    return MultilayerCanopy(; canopy_height, plant_area_index,
        shortwave_model = TwoStreamRadiation(; leaf_reflectance, leaf_transmittance),
        longwave_model = longwave_model_choice, wind_model = canopy_wind_model_choice,
        air_profile_model = canopy_air_profile_model_choice, interception_model = interception_model_choice,
        leaf_parameters = leaf_parameters(site_name), stomatal_model = MoistureResponsiveStomatalConductance(),
        leaf_convection_model = leaf_convection_model_choice, convergence_model = canopy_convergence_model_choice)
end

# MicroModel from the config.jl model choices; `rainfall_schedule` differs by
# forcing source (tower Precip is per-timestep, SILO a daily total).
function micro_model(site_name, site, soil_profile, depths, heights, canopy, rainfall_schedule)
    soil = SoilModel(; depths, profile = soil_profile,
        properties_model = soil_properties_model_choice, hydraulic_model = soil_hydraulic_model(site_name),
        moisture_strategy = soil_moisture_strategy_choice)
    atmosphere = AtmosphereModel(; heights, boundary_layer_model = boundary_layer_model_choice)
    config = MicroConfig(; convergence = convergence_choice, rainfall_schedule)
    return MicroModel(; hours = 0:1:23, site, soil, atmosphere, canopy, config,
        vapour_pressure_equation = vapour_pressure_equation_choice)
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

function initial_moisture_profile(sws_depths, depths, fallback_fraction; scale=1.0)
    known_d, known_v = _known_depth_points(sws_depths)
    isempty(known_d) && return fill(fallback_fraction, length(depths))
    known_v = known_v .* scale  # before clamping -- see sws_obs_scale (utils.jl)
    return [clamp(_interp_flat(known_d, known_v, ustrip(u"m", d)), 0.02, 0.5) for d in depths]
end

function run_site(site_name, years; max_days=nothing, canopy_mode=:full)
    println("\n" * "="^72)
    println("Site $site_name  years=$years  canopy_mode=$canopy_mode")
    println("="^72)
    prep = prepare_site(site_name, years; max_days, canopy_mode)
    println("  Solving MicroProblem ($(prep.n) hours)...")
    solve_time = @elapsed output = Microclimate.solve(prep.problem)
    println("  Solved in $(round(solve_time, digits=1)) s")
    return (; prep..., output, solve_time)
end
