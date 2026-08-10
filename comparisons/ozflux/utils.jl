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

# ── Canopy/soil grid + profile construction (generic, not site-specific) ────

# Soil depths (m): NicheMapR's 10-node scheme (cm) plus a midpoint between
# each pair -- same 19-node default as Microclimate.DEFAULT_DEPTHS.
const NMR_DEP_CM = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]
const SIM_DEPTHS_M = let d = NMR_DEP_CM
    vcat([d[1]], reduce(vcat, [[(d[i] + d[i+1]) / 2, d[i+1]] for i in 1:(length(d) - 1)])) ./ 100.0
end

# base_depths_m plus extra site-specific points (e.g. observed sensor depths).
_build_depths(base_depths_m, extra_depths_m) = sort(unique(vcat(base_depths_m, extra_depths_m))) .* u"m"

# Canopy height grid: fixed absolute spacing, graded near the ground (steepest
# log-law gradient) up to 2m, then uniform coarse steps to canopy_height.
const CANOPY_NEAR_GROUND_M = [1.00, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.00]#[0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50, 2.00]
const CANOPY_COARSE_STEP_M = 1.5
const MIN_CANOPY_LAYERS = 10  # wind_attenuation_profile requires at least this many
const N_ABOVE_CANOPY_POINTS = 6  # graded points from canopy_height up to reference_height

# Pads `pts` to `min_n` points by bisecting the largest gap (anchored at 0).
function _pad_to_min_count(pts, min_n)
    pts = sort(unique(pts))
    while length(pts) < min_n
        anchored = vcat(0.0, pts)
        gaps = diff(anchored)
        i = argmax(gaps)
        push!(pts, (anchored[i] + anchored[i + 1]) / 2)
        sort!(pts); unique!(pts)
    end
    return pts
end

# Full height grid: graded canopy points up to canopy_height, then
# N_ABOVE_CANOPY_POINTS graded up to reference_height (plus any extras and
# reference_height itself, inserted exactly -- reference_height must land on
# a grid point since Microclimate.jl treats last(heights) as reference height).
function _build_heights(canopy_height, reference_height, extra_heights_m=Float64[])
    canopy_m = ustrip(u"m", canopy_height)
    reference_m = ustrip(u"m", reference_height)
    near_ground = filter(<=(canopy_m), CANOPY_NEAR_GROUND_M)
    coarse_top = isempty(near_ground) ? 0.0 : near_ground[end]
    coarse_pts = coarse_top < canopy_m ? collect((coarse_top + CANOPY_COARSE_STEP_M):CANOPY_COARSE_STEP_M:canopy_m) : Float64[]
    canopy_pts = unique(vcat(near_ground, coarse_pts, canopy_m))
    length(canopy_pts) < MIN_CANOPY_LAYERS && (canopy_pts = _pad_to_min_count(canopy_pts, MIN_CANOPY_LAYERS))
    graded_above = canopy_m < reference_m ?
        collect(range(canopy_m, reference_m; length=N_ABOVE_CANOPY_POINTS + 1)[2:end]) : Float64[]
    above_pts = filter(>(canopy_m), unique(vcat(extra_heights_m, graded_above, reference_m)))
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
    mineral_conductivity=1.25u"W/m/K", mineral_heat_capacity=870.0u"J/kg/K",
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
# micropoint's single-slab input alike), for a genuine like-for-like
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

# Scans `data_dir` for `<Site>_<Year>_L3.nc`, returns Dict{String,Vector{Int}}
# of site name => sorted available years. Data-driven rather than a
# hand-maintained site/year list, since files get added over time.
function discover_site_years(data_dir)
    sites = Dict{String,Vector{Int}}()
    for f in readdir(data_dir)
        m = match(r"^(.+)_(\d{4})_L3\.nc$", f)
        m === nothing && continue
        push!(get!(sites, m.captures[1], Int[]), parse(Int, m.captures[2]))
    end
    for years in values(sites)
        sort!(years)
    end
    return sites
end

const OZFLUX_MISSING = (-9999, -9999.0)

# One OzFlux variable, missing-sentinel-replaced.
function _read_var(ds, name::AbstractString)
    sym = Symbol(name)
    haskey(ds, sym) || return nothing
    replace(vec(collect(ds[sym])), OZFLUX_MISSING[1] => missing, OZFLUX_MISSING[2] => missing)
end

# CapeTribulation's RH is a fraction (0-1) for some stretches of the file and
# a percentage (0-100) for others -- a units bug in the source, carrying
# QCFlag=0 ("good") throughout, confirmed by the low stretch clustering
# tightly in 0-0.95 (mean 0.73). Any value <= 1.5 is unambiguously the
# fractional form (percentage RH essentially never falls that low at these
# sites) and is rescaled -- real data recovered, not discarded.
_fix_rh_scale(vals) = [ismissing(v) || v > 1.5 ? v : v * 100.0 for v in vals]

# Physically implausible values that still pass QC (PLAUSIBLE_RANGE,
# config.jl) -- replaced with `missing`, same as the -9999 sentinel, so
# gap-filling (not the raw sensor fault) determines what feeds the model.
function _clip_implausible(name, vals)
    haskey(PLAUSIBLE_RANGE, name) || return vals
    lo, hi = PLAUSIBLE_RANGE[name]
    return [ismissing(v) || (lo <= v <= hi) ? v : missing for v in vals]
end

# Stuck-sensor detection: FLATLINE_MIN_RUN+ identical consecutive raw
# readings is a fault signature real sensor noise essentially never produces
# -- not caught by a range check (a plausible-looking value can still be
# stuck). Exact zero is excluded: calm wind, night Fsd, and dry spells are
# all legitimately flat at zero for hours at a time.
const FLATLINE_MIN_RUN = 8  # raw 30-min steps = 4 hours

function _clip_flatline(vals; min_run=FLATLINE_MIN_RUN)
    out = collect(vals)
    n = length(out)
    i = 1
    while i <= n
        j = i
        while j < n && !ismissing(out[i]) && !ismissing(out[j+1]) && out[j+1] == out[i]
            j += 1
        end
        if !ismissing(out[i]) && out[i] != 0.0 && (j - i + 1) >= min_run
            out[i:j] .= missing
        end
        i = j + 1
    end
    return out
end

# _read_var + RH scale fix + PLAUSIBLE_RANGE clip + flatline clip, in one
# place so every read site (forcing/target vars, height series, Ah) gets the
# same cleaning rather than only the main forcing loop.
function _read_clean_var(ds, name::AbstractString)
    vals = _read_var(ds, name)
    vals === nothing && return nothing
    name == "RH" && (vals = _fix_rh_scale(vals))
    return _clip_flatline(_clip_implausible(name, vals))
end

# Depth-suffixed series for `prefix` (e.g. "Sws", "Ts"): variable names
# `<prefix>_<depth>(cm|m)`, optionally followed by a single replicate letter
# (a/b/c, ...). A bare (no-letter) variable is the site's own official
# per-depth average and is used directly; when only replicates exist (e.g.
# Calperum's Ts has no bare Ts_10cm, only Ts_10cma/b/c) they're averaged
# here, element-wise, skipping missing.
function discover_depth_series(ds, prefix::AbstractString)
    varnames = String.(collect(names(ds)))
    pattern = Regex("^" * prefix * raw"_(\d+(?:\.\d+)?)(cm|m)([a-z]?)$")
    bare = Dict{Float64,String}()
    replicates = Dict{Float64,Vector{String}}()
    for name in varnames
        m = match(pattern, name)
        m === nothing && continue
        depth_m = parse(Float64, m.captures[1]) * (m.captures[2] == "cm" ? 0.01 : 1.0)
        if m.captures[3] == ""
            bare[depth_m] = name
        else
            push!(get!(replicates, depth_m, String[]), name)
        end
    end
    depths_m = sort(collect(union(keys(bare), keys(replicates))))
    result = NamedTuple{(:depth_m, :values),Tuple{Float64,Vector{Union{Float64,Missing}}}}[]
    for d in depths_m
        values = if haskey(bare, d)
            _read_var(ds, bare[d])
        else
            cols = [_read_var(ds, name) for name in replicates[d]]
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

# Similarly discovers height-suffixed variants of a base variable (e.g. Ws ->
# Ws_RMY2m_Av at 2m, Ws_RMY10m_Av at 10m, bare Ws / Ws_SONIC_Av at tower
# height) for the multi-height profile comparison. Site-specific naming (no
# single convention across sites), so this takes an explicit list of
# (height_m, variable_name) pairs rather than trying to auto-detect them.
function read_height_series(ds, entries)
    [(height_m=h, values=_read_clean_var(ds, name)) for (h, name) in entries if _read_clean_var(ds, name) !== nothing]
end

# Reads one OzFlux L3 file's forcing + validation variables, per-depth Ts/Sws
# series, and any requested multi-height series into one DataFrame (one row
# per native timestep, missing-replaced, `<name>_QCFlag` columns carried
# alongside each variable so _aggregate_hourly's QC filtering applies to all
# of them uniformly). `height_series`: (height_m, variable_name) pairs from
# config.jl's SITE_HEIGHT_SERIES; height_m itself isn't needed here, only
# the variable name.
function read_ozflux_nc(path; height_series=Tuple{Float64,String}[])
    ds = RasterStack(path)
    time = DateTime.(lookup(ds, Ti))
    df = DataFrame(DateTime = time)

    forcing_vars = ("Ta", "RH", "Ws", "Fsd", "Fld", "Precip", "ps")
    target_vars = ("Fsu", "Flu", "Fn", "Fh", "Fe", "Fg")
    for v in (forcing_vars..., target_vars..., last.(height_series)...)
        vals = _read_clean_var(ds, v)
        vals === nothing && continue
        df[!, v] = vals
        qc = _read_var(ds, v * "_QCFlag")
        qc === nothing || (df[!, v * "_QCFlag"] = qc)
    end

    # minimum wind speed capped at 0.1 m/s (sensor limit)
    wind_cols = filter(names(df)) do col
        (col == "Ws" || startswith(col, "Ws_")) &&
            !endswith(col, "_QCFlag") # reject QCFlag columns
    end

    for col in wind_cols
        df[!, col] = map(x -> ismissing(x) ? missing : max(x, 0.1), df[!, col])
    end 

    # Absolute humidity, cased inconsistently across sites (Ah vs AH) --
    # normalized to "Ah" (g/m^3). Used in pipeline.jl as an RH cross-check
    # and secondary gap-fill source (same tower, different instrument).
    ah_name = haskey(ds, :Ah) ? "Ah" : haskey(ds, :AH) ? "AH" : nothing
    if ah_name !== nothing
        df[!, "Ah"] = _read_clean_var(ds, ah_name)
        ah_qc = _read_var(ds, ah_name * "_QCFlag")
        ah_qc === nothing || (df[!, "Ah_QCFlag"] = ah_qc)
    end

    # Depth-suffixed soil series get renamed (Ts_10cma/b/c -> Ts_0.1m) since
    # discover_depth_series already resolves replicates/units; no per-replicate
    # QC filtering here (sentinel-missing only) -- see README caveats.
    for prefix in ("Ts", "Sws")
        for (depth_m, values) in discover_depth_series(ds, prefix)
            df[!, "$(prefix)_$(depth_m)m"] = values
        end
    end

    return df
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
