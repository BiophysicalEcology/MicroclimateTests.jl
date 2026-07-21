# utils.jl — generic analysis helpers for comparison.jl / debug_site.jl

# Raw per-day series for one input of a bound DielForcing (e.g. min series).
function _minmax_series(env, name::Symbol, want::Symbol)
    spec = getproperty(Microclimate.MINMAX_FORCING_MODEL, name)
    idx = findfirst(==(want), Microclimate.forcing_inputs(spec))
    return getproperty(env.forcings, name).values[idx]
end

# Per-day (min, max) from a sub-daily series, for sources like BARRA that are
# natively sub-daily -- Microclimate.jl solves straight off environment_hourly
# for these and leaves environment_minmax = nothing, but NicheMapR (and our
# own day-1 initial-condition estimate) still need daily min/max.
function _daily_minmax(series::AbstractVector, ndays::Int)
    per_day = length(series) ÷ ndays
    mn = similar(series, ndays)
    mx = similar(series, ndays)
    @inbounds for d in 1:ndays
        lo, hi = extrema(view(series, (d - 1) * per_day + 1 : d * per_day))
        mn[d] = lo
        mx[d] = hi
    end
    return mn, mx
end

# Linear interpolation through (xs, ys) at x; clamps to endpoints outside range.
# Works for any numeric ys, including Unitful quantities.
function _linterp(xs::Vector{Float64}, ys::Vector, x::Float64)
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i] + t * (ys[i+1] - ys[i])
end

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

# Match hourly model output to irregularly-sampled obs by rounding obs to nearest hour.
function align_model_to_obs(obs_dt::Vector{DateTime},
                             t_model::Vector{DateTime},
                             model_cols::Dict{Symbol, Vector{Float64}})
    idx_map = Dict(t => i for (i, t) in enumerate(t_model))
    result  = Dict(sym => Vector{Union{Float64,Missing}}(missing, length(obs_dt))
                   for sym in keys(model_cols))
    for (row, dt) in enumerate(obs_dt)
        idx = get(idx_map, round(dt, Hour(1)), 0)
        idx == 0 && continue
        for (sym, vec) in model_cols
            idx <= length(vec) && (result[sym][row] = vec[idx])
        end
    end
    result
end

# Remove spikes from an obs vector (Union{Float64,Missing}).
#
# Algorithm:
#   1. Optional: set negative values to 0 (clip_negative=true).
#   2. Optional: remove values above the (1-tail_quantile) quantile (tail_quantile > 0).
#   3. Flag any value whose absolute difference from ANY non-missing neighbour
#      within ±`halfwin` steps exceeds `thresh`; flag both the value and its
#      offending neighbour.  Flagging is done simultaneously (not sequentially)
#      so the result is order-independent.
#   4. Replace flagged values with missing.
#
# Returns a copy; does not modify in place.
function remove_spikes(y::AbstractVector;
                       thresh::Real,
                       halfwin::Int       = 2,
                       clip_negative::Bool = false,
                       tail_quantile::Real = 0.0)
    v = Vector{Union{Float64,Missing}}(y)
    n = length(v)

    if clip_negative
        for i in 1:n
            !ismissing(v[i]) && v[i] < 0 && (v[i] = 0.0)
        end
    end

    if tail_quantile > 0
        vals = Float64[x for x in v if !ismissing(x)]
        if length(vals) >= 10
            upper = quantile(vals, 1.0 - tail_quantile)
            for i in 1:n
                !ismissing(v[i]) && v[i] > upper && (v[i] = missing)
            end
        end
    end

    flagged = falses(n)
    for i in 1:n
        ismissing(v[i]) && continue
        for Δ in -halfwin:halfwin
            Δ == 0 && continue
            j = i + Δ
            (j < 1 || j > n) && continue
            ismissing(v[j]) && continue
            if abs(v[i] - v[j]) > thresh
                flagged[i] = true
                flagged[j] = true
            end
        end
    end
    for i in 1:n
        flagged[i] && (v[i] = missing)
    end
    v
end

# Apply remove_spikes to named columns of a DataFrame in-place.
# col_specs: iterable of (colname, thresh; kwargs...) or (colname => (thresh=…, …)).
# Simpler: pass a vector of (Symbol, NamedTuple) pairs.
function despiked!(df::DataFrame, specs::Vector{<:Pair})
    for (col, kw) in specs
        col in propertynames(df) || continue
        n_before = count(!ismissing, df[!, col])
        df[!, col] = remove_spikes(df[!, col]; kw...)
        n_after  = count(!ismissing, df[!, col])
        n_before - n_after > 0 &&
            println("    despiked $col: removed $(n_before - n_after) values")
    end
end

# Read a CSV if it exists; return nothing and print a notice if absent.
function try_read(path, label)
    if isfile(path)
        df = DataFrame(CSV.File(path, normalizenames=true))
        println("  Found $label ($(nrow(df)) rows)")
        return df
    else
        println("  Not found: $label — NMR comparison skipped for this variable")
        return nothing
    end
end
