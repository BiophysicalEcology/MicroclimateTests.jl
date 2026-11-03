# utils.jl — generic analysis helpers (shared with oznet/) plus the OzFlux
# NetCDF reader/aggregator, shared by comparison.jl / single_site.jl.

# ── Reused verbatim from oznet/utils.jl ──────────────────────────────────────

struct ModelStats
    r    :: Float64
    rmse :: Float64
    bias :: Float64   # mean(pred − obs)
    n    :: Int
end

function compute_stats(obs_vec, pred_vec)
    pairs = [(Float64(o), Float64(p))
             for (o, p) in zip(obs_vec, pred_vec)
             if !ismissing(o) && !ismissing(p) && isfinite(Float64(o)) && isfinite(Float64(p))]
    length(pairs) < 2 && return ModelStats(NaN, NaN, NaN, 0)
    o = first.(pairs);  p = last.(pairs)
    ModelStats(cor(o, p), sqrt(mean((p .- o).^2)), mean(p .- o), length(pairs))
end

fmt_stat(s::ModelStats) =
    s.n < 2 ? "         —         " :
    @sprintf("r=%+.3f  RMSE=%6.3f  bias=%+7.4f  (n=%d)", s.r, s.rmse, s.bias, s.n)

# Factor taking a site's observed Sws to volumetric water content: the top layer's porosity
# where data_corrections.jl's SaturationFraction says it is a fraction of saturation. Applied
# to both the obs comparison (report.jl) and the initial-condition read (pipeline.jl).
sws_obs_scale(site_name, soil_profile) = has_correction(site_corrections(site_name), SaturationFraction, "Sws") ?
    1 - soil_profile.bulk_density[1] / soil_profile.mineral_density[1] : 1.0

# Truncates a DailyMinMaxEnvironment's per-day forcing series to the first
# `ndays` days. Needed by prepare_site_silo's max_days path: fetch_silo_forcing
# fetches whole years, but the solve itself only iterates
# `ndays` days -- left untruncated, environment_minmax's forcing arrays stay
# longer than solar_radiation_out (sized to `days`), a BoundsError in
# evaluate!/solar_day once the solver reaches day ndays+1.
# minmax_forcings' derived _mean fields are lazy Broadcasted (not materialized
# Vectors, so they track later in-place mutation of the min/max arrays --
# see its own comment) -- materialize before slicing since Broadcasted
# doesn't support arbitrary indexing. Safe here: any such mutation is
# internal to fetch_silo_forcing, already done by the time it returns.
_truncate_days(v::AbstractVector, ndays) = v[1:ndays]
_truncate_days(v, ndays) = Base.Broadcast.materialize(v)[1:ndays]
_truncate_forcing(f::DielForcing, ndays) = DielForcing(f.curve, map(v -> _truncate_days(v, ndays), f.values))
_truncate_forcing(d::Derived, ndays) = d
truncate_environment_minmax(env::DailyMinMaxEnvironment, ndays) =
    DailyMinMaxEnvironment(; forcings = map(f -> _truncate_forcing(f, ndays), env.forcings))

# Same max_days issue, one level down: fetch_silo_forcing's environment_hourly
# (the SILO stub -- only pressure is ever populated, everything else nothing)
# is also fetched full-year and needs truncating to `n` hours to match the
# truncated solve. Generic over whichever fields are actually Vectors, so it
# doesn't assume pressure is the only one.
function truncate_environment_hourly(h::HourlyTimeseries, n)
    kwargs = NamedTuple(f => (v = getfield(h, f); v isa AbstractVector ? v[1:n] : v) for f in fieldnames(typeof(h)))
    return HourlyTimeseries(; kwargs...)
end

# ── Canopy/soil grid + profile construction (generic, not site-specific) ────

# Soil depths (m): NicheMapR's 10-node scheme (cm) plus a midpoint between
# each pair -- same 19-node default as Microclimate.DEFAULT_DEPTHS.
const NMR_DEP_CM = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]
const SIM_DEPTHS_M = let d = NMR_DEP_CM
    vcat([d[1]], reduce(vcat, [[(d[i] + d[i+1]) / 2, d[i+1]] for i in 1:(length(d) - 1)])) ./ 100.0
end

# base_depths_m plus extra site-specific points (e.g. observed sensor depths).
_build_depths(base_depths_m, extra_depths_m) = sort(unique(vcat(base_depths_m, extra_depths_m))) .* u"m"

# Canopy height grid: N_CANOPY_LAYERS evenly-spaced points from 0 to
# canopy_height -- the same scheme micropoint's own R model uses (see
# write_ozflux_micropoint_inputs.jl's micropoint_canopy_layers), matched at
# the same N=50 so the two models' PAI/layer structure is directly
# comparable. A single special-cased near-ground insertion (kept resolution
# to ~15 layers otherwise) was tried first but created a genuine artifact:
# splitting one ~1.87m bottom layer into a normal layer plus one very thin
# (0.1m) one made that thin layer's *total* PAI look anomalously low even
# though its density was correctly high -- an unequal-layer-thickness
# artifact, not a bug in the PAI math. Uniform N=50 spacing avoids it
# entirely (matches micropoint's own approach) at the cost of more layers.
const N_CANOPY_LAYERS = 20  # wind_attenuation_profile requires at least 10
const N_ABOVE_CANOPY_POINTS = 6  # graded points from canopy_height up to reference_height

# Full height grid: N_CANOPY_LAYERS evenly-spaced canopy points up to
# canopy_height, then N_ABOVE_CANOPY_POINTS graded up to reference_height
# (reference_height inserted exactly -- Microclimate.jl treats last(heights)
# as reference height).
const NEAR_ZERO_FLOOR_M = 0.01  # mm-scale spacing has crashed micropoint's solver before

# Forward-floors an ascending grid to a minimum step, starting from pts[1]
# itself (no implicit zero anchor). Self-limiting: only touches points where
# natural spacing is already below the floor. Skipped entirely (natural,
# possibly sub-floor spacing kept) when the grid's own span can't fit
# length(pts)-1 full floor-sized steps -- otherwise the cascade would push
# the last point past the domain's own top (e.g. a canopy shorter than
# length(pts)*floor_m); that would break the pts[end] == canopy_height
# invariant callers rely on, worse than the mm-scale spacing this guards
# against. Not triggered by any current site (shortest canopy is 0.5m,
# exactly the no-op boundary at N_CANOPY_LAYERS=50/NEAR_ZERO_FLOOR_M=0.01).
function _floor_ascending!(pts, floor_m)
    n = length(pts)
    (pts[n] - pts[1]) < floor_m * (n - 1) && return pts
    pts[1] = max(pts[1], floor_m)
    for k in 2:n
        pts[k] = max(pts[k], pts[k - 1] + floor_m)
    end
    return pts
end

function _build_heights(canopy_height, reference_height)
    canopy_m = ustrip(u"m", canopy_height)
    reference_m = ustrip(u"m", reference_height)
    canopy_pts = _floor_ascending!([(i / N_CANOPY_LAYERS) * canopy_m for i in 1:N_CANOPY_LAYERS], NEAR_ZERO_FLOOR_M)
    graded_above = canopy_m < reference_m ?
        collect(range(canopy_m, reference_m; length=N_ABOVE_CANOPY_POINTS + 1)[2:end]) : Float64[]
    above_pts = filter(>(canopy_m), unique(vcat(graded_above, reference_m)))
    return sort(vcat(canopy_pts, above_pts)) .* u"m"
end

# Vertical PAI density shapes (unnormalized, per-metre); plant_area_index_profile rescales to target_pai.
_pai_shape_top_heavy(layer_heights, canopy_height) = @. exp(3.0 * (layer_heights / canopy_height)) / u"m"
_pai_shape_bottom_heavy(layer_heights, canopy_height) = @. exp(-3.0 * (layer_heights / canopy_height)) / u"m"
_pai_shape_mid_crown(layer_heights, canopy_height) =
    @. exp(-0.5 * ((layer_heights - 0.7 * canopy_height) / (0.25 * canopy_height))^2) / u"m"
_pai_shape_uniform(layer_heights, canopy_height) = fill(1.0 / u"m", length(layer_heights))
const PAI_SHAPES = Dict(
    :top_heavy => _pai_shape_top_heavy, :bottom_heavy => _pai_shape_bottom_heavy,
    :mid_crown => _pai_shape_mid_crown, :uniform => _pai_shape_uniform,
)

# Per-layer PAI vector for MultilayerCanopy's plant_area_index, from a named shape rescaled to sum to target_pai.
function plant_area_index_profile(shape_kind, heights, canopy_height, target_pai)
    n_layers = count(h -> h <= canopy_height, heights)
    (; layer_heights) = Microclimate.canopy_layer_heights(heights, canopy_height, n_layers)
    density = PAI_SHAPES[shape_kind](layer_heights, canopy_height)
    raw = plant_area_index_from_density(density, heights, canopy_height)
    return raw .* (target_pai / sum(raw))
end

# example_campbell_hydraulic_profile()'s root_density is fixed to DEFAULT_DEPTHS's
# 19 nodes -- interpolate onto arbitrary `depths` instead of reusing it directly.
function _root_density_for_depths(depths)
    ref = Microclimate.example_campbell_hydraulic_profile()
    ref_depths_m = ustrip.(u"m", Microclimate.DEFAULT_DEPTHS)
    ref_root = ustrip.(u"m/m^3", ref.root_density)
    target_m = ustrip.(u"m", depths)
    n = length(ref_depths_m)
    interp1(x) = x <= ref_depths_m[1] ? ref_root[1] :
                 x >= ref_depths_m[n] ? ref_root[n] :
                 let i = clamp(searchsortedlast(ref_depths_m, x), 1, n - 1)
                     t = (x - ref_depths_m[i]) / (ref_depths_m[i+1] - ref_depths_m[i])
                     ref_root[i] + t * (ref_root[i+1] - ref_root[i])
                 end
    return interp1.(target_m) .* u"m/m^3"
end

# air_entry_water_potential stored positive (matches SLGA's sign convention).
function soil_profile_from_texture(texture::NamedTuple, depths;
    bulk_density=1.3u"Mg/m^3", mineral_density=2.560u"Mg/m^3",
    mineral_conductivity=2.5u"W/m/K", mineral_heat_capacity=870.0u"J/kg/K",
    root_density=_root_density_for_depths(depths),
)
    n = length(depths)
    SoilProfile(;
        bulk_density=fill(bulk_density, n), mineral_density=fill(mineral_density, n),
        mineral_conductivity=fill(mineral_conductivity, n), mineral_heat_capacity=fill(mineral_heat_capacity, n),
        hydraulics=CampbellHydraulicProfile(;
            air_entry_water_potential=fill(texture.air_entry, n),
            saturated_hydraulic_conductivity=fill(texture.Ksat, n),
            campbell_b_parameter=fill(texture.b, n), root_density,
        ),
    )
end

# SoilProfile from measured soil horizons (tables/sites.toml's soil_horizons):
# each node takes the bulk density and texture of the horizon it sits in, with
# Campbell hydraulics from `pedotransfer_model`. Nodes below the deepest horizon
# take the deepest horizon's values.
function soil_profile_from_horizons(horizons, depths, pedotransfer_model;
    mineral_density=2.560u"Mg/m^3", mineral_conductivity=2.5u"W/m/K", mineral_heat_capacity=870.0u"J/kg/K",
    root_density=_root_density_for_depths(depths),
)
    n = length(depths)
    horizon = [min(searchsortedfirst(horizons.bottom, d), length(horizons.bottom)) for d in depths]
    bulk_density = horizons.bulk_density[horizon]
    SoilProfile(;
        bulk_density, mineral_density=fill(mineral_density, n),
        mineral_conductivity=fill(mineral_conductivity, n), mineral_heat_capacity=fill(mineral_heat_capacity, n),
        hydraulics=CampbellHydraulicProfile(pedotransfer_model; clay_percent=horizons.clay[horizon],
            silt_percent=horizons.silt[horizon], sand_percent=horizons.sand[horizon], bulk_density, root_density),
    )
end

# Layer thickness per node (m), from midpoint boundaries between neighbouring
# depths, clipped to [0, totaldepth_m] at the ends -- used to depth-weight an
# uneven depths grid (SIM_DEPTHS_M packs many nodes near the surface) rather
# than let shallow nodes dominate an unweighted per-node mean.
function _depth_weights(depths_m, totaldepth_m)
    n = length(depths_m)
    edges = n == 1 ? [0.0, totaldepth_m] :
        [0.0; [(depths_m[i] + depths_m[i + 1]) / 2 for i in 1:(n - 1)]; totaldepth_m]
    return diff(edges)
end

_wmean(vals::AbstractVector{<:Quantity}, weights) =
    (sum(weights .* ustrip.(unit(first(vals)), vals)) / sum(weights)) * unit(first(vals))
_wmean(vals::AbstractVector{<:Real}, weights) = sum(weights .* vals) / sum(weights)
_wgeomean(vals::AbstractVector{<:Quantity}, weights) =
    exp(sum(weights .* log.(ustrip.(unit(first(vals)), vals))) / sum(weights)) * unit(first(vals))

# Depth-weighted mean of `soil_profile` over 0..active_depth, returned as a
# new SoilProfile with every field replaced by its own uniform mean but on
# the SAME depths grid (same array length) -- a drop-in replacement wherever
# the depth-varying original was used (Microclimate.jl's solver or
# micropoint's single-slab input alike), for a like-for-like
# comparison. active_depth defaults to 0.3m -- the near-surface layer that
# dominates diurnal heat/moisture exchange -- not micropoint's own
# totalDepth=2m column parameter (just where its deep boundary condition
# sits). Averaging over the full column instead pulls in deeper, denser,
# more compacted subsoil and understates near-surface porosity. Saturated_
# hydraulic_conductivity uses a weighted geometric mean (spans orders of
# magnitude across a real profile); everything else is weighted-arithmetic.
# root_density is left depth-varying -- a rooting profile, not a soil
# property, still meaningful under an otherwise-uniform soil.
function flatten_soil_profile(soil_profile, depths; active_depth=0.3u"m")
    depths_m = ustrip.(u"m", depths)
    td = ustrip(u"m", active_depth)
    nodes = findall(<=(td), depths_m)
    weights = _depth_weights(depths_m[nodes], td)
    n = length(depths)
    fillrep(x) = fill(x, n)
    return SoilProfile(;
        bulk_density = fillrep(_wmean(soil_profile.bulk_density[nodes], weights)),
        mineral_density = fillrep(_wmean(soil_profile.mineral_density[nodes], weights)),
        mineral_conductivity = fillrep(_wmean(soil_profile.mineral_conductivity[nodes], weights)),
        mineral_heat_capacity = fillrep(_wmean(soil_profile.mineral_heat_capacity[nodes], weights)),
        hydraulics = CampbellHydraulicProfile(;
            air_entry_water_potential = fillrep(_wmean(soil_profile.hydraulics.air_entry_water_potential[nodes], weights)),
            saturated_hydraulic_conductivity = fillrep(_wgeomean(soil_profile.hydraulics.saturated_hydraulic_conductivity[nodes], weights)),
            campbell_b_parameter = fillrep(_wmean(soil_profile.hydraulics.campbell_b_parameter[nodes], weights)),
            root_density = soil_profile.hydraulics.root_density,
        ),
    )
end

# ── OzFlux NetCDF reading ─────────────────────────────────────────────────────

# One OzFlux L3 NetCDF per site holding its whole record (TERN OzFlux
# releases), named `<Site>[<n>]_L3.nc` (e.g. Whroo1_L3.nc). Site names carry
# no digits, so a trailing number is the file's own suffix, not the site.
const OZFLUX_FILE_PATTERN = r"^([A-Za-z]+)\d*_L3\.nc$"

function ozflux_site_path(data_dir, site_name)
    matches = filter(f -> (m = match(OZFLUX_FILE_PATTERN, f); m !== nothing && m.captures[1] == site_name), readdir(data_dir))
    length(matches) == 1 || error("expected one OzFlux L3 file for $site_name in $data_dir, found $(length(matches)): $matches")
    return joinpath(data_dir, only(matches))
end

# Calendar years with at least one record in an OzFlux file.
ozflux_years(path) = NCDataset(ds -> sort(unique(Dates.year.(ds["time"][:]))), path)

# Scans `data_dir` for OzFlux files, returns Dict{String,Vector{Int}} of site
# name => sorted available years (from each file's own time axis).
# Data-driven rather than a hand-maintained site/year list, since files get
# added over time.
function discover_site_years(data_dir)
    sites = Dict{String,Vector{Int}}()
    for f in readdir(data_dir)
        m = match(OZFLUX_FILE_PATTERN, f)
        m === nothing && continue
        sites[m.captures[1]] = ozflux_years(joinpath(data_dir, f))
    end
    return sites
end

const OZFLUX_MISSING = (-9999, -9999.0)

# One OzFlux variable over the record rows `rows`, missing-sentinel-replaced.
function _read_var(ds, name::AbstractString, rows)
    haskey(ds, name) || return nothing
    replace(Vector{Union{Missing,Float64}}(ds[name][rows]), OZFLUX_MISSING[1] => missing, OZFLUX_MISSING[2] => missing)
end

# _read_var with data_corrections.jl's `corrections` applied, so every series read (forcing,
# targets, height and depth series, Ah) gets the same corrections
function _read_clean_var(ds, name::AbstractString, rows, corrections)
    vals = _read_var(ds, name, rows)
    vals === nothing && return nothing
    return correct_values(corrections, name, vals, _time_step(ds))
end

# The file's record interval (30 min at most sites, 60 at Otway and Tumbarumba)
_time_step(ds) = Minute(parse(Int, string(ds.attrib["time_step"])))

# Heights (m) in an OzFlux variable's `height` attribute ("36m", "8.8 m", "30cm", "-0.08m",
# CumberlandPlain's "-<.08m", Yanco's "8m,4m"), as distances from the ground surface
function attribute_heights(ds, name)
    haskey(ds, name) || return Float64[]
    text = string(get(ds[name].attrib, "height", ""))
    return [parse(Float64, m[1]) * (m[2] == "cm" ? 0.01 : 1.0) for m in eachmatch(r"(\d+(?:\.\d*)?|\.\d+)\s*(cm|m)?", text)]
end

# The single height (m) of an OzFlux variable; nothing when its attribute gives none or several
function variable_height(ds, name)
    heights = attribute_heights(ds, name)
    return length(heights) == 1 ? only(heights) : nothing
end

# Depth (m) a soil sensor measures at, given the depth it is named or labelled for. A TCAV
# averaging thermocouple is named for the heat flux plate depth (8 cm) but its junctions sit
# at 2 and 6 cm above the plates (Campbell Scientific HFT3/TCAV installation), so it measures
# the 2-6 cm mean, ~4 cm.
_sensor_depth(ds, name, depth_m) = occursin("TCAV", get(ds[name].attrib, "instrument", "")) ? 0.04 : depth_m

# Depth-suffixed series for `prefix` (e.g. "Sws", "Ts"): variable names
# `<prefix>_<depth>(cm|m)`, optionally followed by a single replicate letter
# (a/b/c, ...). A bare (no-letter) variable is the site's own official
# per-depth average and is used directly; when only replicates exist (e.g.
# Calperum's Ts has no bare Ts_10cm, only Ts_10cma/b/c) they're averaged
# here, element-wise, skipping missing. Files with no depth-suffixed series
# (most current TERN releases) give the bare `<prefix>` composite at its
# `height` attribute's depth.
function discover_depth_series(ds, prefix::AbstractString, rows, corrections)
    varnames = keys(ds)
    pattern = Regex("^" * prefix * raw"_(\d+(?:\.\d+)?)(cm|m)([a-z]?)$")
    bare = Dict{Float64,String}()
    replicates = Dict{Float64,Vector{String}}()
    for name in varnames
        m = match(pattern, name)
        m === nothing && continue
        depth_m = _sensor_depth(ds, name, parse(Float64, m.captures[1]) * (m.captures[2] == "cm" ? 0.01 : 1.0))
        if m.captures[3] == ""
            bare[depth_m] = name
        else
            push!(get!(replicates, depth_m, String[]), name)
        end
    end
    if isempty(bare) && isempty(replicates)
        depth_m = variable_height(ds, prefix)
        depth_m === nothing || (bare[_sensor_depth(ds, prefix, depth_m)] = prefix)
    end
    depths_m = sort(collect(union(keys(bare), keys(replicates))))
    result = NamedTuple{(:depth_m, :values),Tuple{Float64,Vector{Union{Float64,Missing}}}}[]
    for d in depths_m
        values = if haskey(bare, d)
            _read_clean_var(ds, bare[d], rows, corrections)
        else
            cols = [_read_clean_var(ds, name, rows, corrections) for name in replicates[d]]
            _rowwise_mean(cols)
        end
        values === nothing || push!(result, (depth_m=d, values=values))
    end
    return result
end

# Element-wise mean across columns, skipping missing (NaN if all missing at a row).
function _rowwise_mean(cols::Vector)
    n = length(cols[1])
    out = Vector{Union{Float64,Missing}}(missing, n)
    @inbounds for i in 1:n
        vals = Float64[]
        for col in cols
            v = col[i]
            ismissing(v) || push!(vals, v)
        end
        isempty(vals) || (out[i] = mean(vals))
    end
    return out
end

# Reads the records of `years` from one OzFlux L3 file's validation
# variables, per-depth Ts/Sws series, and any requested multi-height series
# into one DataFrame (one row per native timestep, missing-replaced,
# `<name>_QCFlag` columns carried alongside each variable so
# _aggregate_hourly's QC filtering applies to all of them uniformly). A
# record belongs to the year of the hour it averages into
# (_hour_group_start, pipeline.jl), so a year's last hour keeps its
# end-labelled 00:00 record from the next year. `height_series`:
# (height_m, variable_name) pairs from sites.jl's discover_height_series;
# height_m itself isn't needed here, only the variable name. `extra_vars`:
# further variable names to read: the forcing sensors forcing_sources picked and the site's
# forcing_backups (tables/sites.toml). Every series gets
# the site's data_corrections.jl `corrections`.
function read_ozflux_nc(path, years, corrections; height_series=Tuple{Float64,String}[], extra_vars=String[])
    NCDataset(path) do ds
        all_time = ds["time"][:]
        rows = findall(t -> Dates.year(_hour_group_start(t)) in years, all_time)
        isempty(rows) && error("$path has no records in $years")
        rows == first(rows):last(rows) || error("$path time axis is not sorted")
        rows = first(rows):last(rows)
        df = DataFrame(DateTime = all_time[rows])
        DataFrames.metadata!(df, "time_step", _time_step(ds); style=:note)

        target_vars = ("Fsu", "Flu", "Fn", "Fh", "Fe", "Fg", "ustar")
        for v in unique((target_vars..., last.(height_series)..., extra_vars...))
            vals = _read_clean_var(ds, v, rows, corrections)
            vals === nothing && continue
            df[!, v] = vals
            qc = _read_var(ds, v * "_QCFlag", rows)
            qc === nothing || (df[!, v * "_QCFlag"] = qc)
        end

        # Depth-suffixed soil series get renamed (Ts_10cma/b/c -> Ts_0.1m) since
        # discover_depth_series already resolves replicates/units; no per-replicate
        # QC filtering here (sentinel-missing only) -- see README caveats.
        for prefix in ("Ts", "Sws")
            for (depth_m, values) in discover_depth_series(ds, prefix, rows, corrections)
                df[!, "$(prefix)_$(depth_m)m"] = values
            end
        end

        return df
    end
end

# Parses e.g. "25m" / "3.5 m" / "25" -> 25.0u"m". Range attributes (Wallaby's
# canopy_height = "8.0 - 10.0m") take the first (lower-bound) number.
function _parse_length_attr(s)
    m = match(r"([\d.]+)", string(s))
    m === nothing && throw(ArgumentError("could not parse a length from $(repr(s))"))
    parse(Float64, m.captures[1]) * u"m"
end

# Read a CSV if it exists; return nothing and print a notice if absent.
function try_read(path, label)
    if isfile(path)
        df = DataFrame(CSV.File(path, normalizenames=true))
        println("  Found $label ($(nrow(df)) rows)")
        return df
    else
        println("  Not found: $label")
        return nothing
    end
end