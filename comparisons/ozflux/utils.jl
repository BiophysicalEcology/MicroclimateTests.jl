# utils.jl — generic analysis helpers (shared with oznet/) plus the OzFlux
# NetCDF reader/aggregator, shared by comparison.jl / debug_site.jl.

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
