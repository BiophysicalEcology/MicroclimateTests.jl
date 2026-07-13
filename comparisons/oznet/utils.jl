# utils.jl — generic analysis helpers + the OzNet observation-file parser,
# shared by comparison.jl / debug_site.jl.

# Raw per-day series for one input of a bound DielForcing (e.g. min series).
function _minmax_series(env, name::Symbol, want::Symbol)
    spec = getproperty(Microclimate.MINMAX_FORCING_MODEL, name)
    idx = findfirst(==(want), Microclimate.forcing_inputs(spec))
    return getproperty(env.forcings, name).values[idx]
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

# ── OzNet observation file parser ────────────────────────────────────────────
#
# File layout: line 1 = site display name (ignored), line 2 = comma-separated
# header (ignored -- column *count*, plus the site's own `temp.depths.cm`/
# `soil.moisture.depths.cm` from oznetsiteinfo2.csv, are what actually
# determine column meaning, not this header text. The header's depth labels
# don't always match the true depths for sites sharing a column layout --
# e.g. site m5's true deepest temperature depth is 62.5 cm even though its
# file uses the same 15-column layout as sites whose header says 75 cm.)
#
# Data rows are whitespace-delimited:
#   "DD/MM/YYYY HH:MM <temp cols> <moisture cols> [<suction cols>] <rainfall>"
# Across every OzNet site there are exactly two layouts (verified against all
# 38 files): 9 data columns (2 temp, 4 moisture, no suction) or 15 data
# columns (4 temp, 4 moisture, 4 suction -- suction is read but dropped, not
# used by either model here). `temp_depths_cm` labels the first
# `length(temp_depths_cm)` temperature columns; any temperature column beyond
# that (only site y3, whose file has 4 temp columns but whose metadata lists
# only 3 depths) is left unlabelled and silently dropped, matching the
# original NicheMapR-based study's behaviour (it only ever loops over the
# listed depths). `moist_depths_cm` always has exactly 4 entries.
#
# Missing/invalid readings are coded -99 in the raw files. Following the
# original study's thresholds: temperature < -20 °C and moisture < 0 %vol are
# both physically impossible at these sites and are converted to `missing`.
function read_oznet_obs(path, temp_depths_cm::Vector{Float64}, moist_depths_cm::Vector{Float64})
    lines = readlines(path)
    data_lines = @view lines[3:end]
    isempty(data_lines) && throw(ArgumentError("$path: no data rows"))

    n_fields = length(split(first(data_lines)))
    n_data_cols = n_fields - 2   # minus DATE, TIME
    has_suction = n_data_cols == 15 - 2
    n_temp_cols = has_suction ? 4 : 2
    n_moist_cols = 4
    n_suction_cols = has_suction ? 4 : 0
    expected = n_temp_cols + n_moist_cols + n_suction_cols + 1
    n_data_cols == expected || throw(ArgumentError(
        "$path: $n_data_cols data column(s) (unrecognised layout; expected 9 or 15 total incl. date+time)"))
    length(temp_depths_cm) <= n_temp_cols || throw(ArgumentError(
        "$path: $(length(temp_depths_cm)) labelled temp depth(s) but only $n_temp_cols temp column(s) in file"))

    ndates = length(data_lines)
    dt = Vector{DateTime}(undef, ndates)
    temp_vals  = [Vector{Union{Float64,Missing}}(undef, ndates) for _ in temp_depths_cm]
    moist_vals = [Vector{Union{Float64,Missing}}(undef, ndates) for _ in moist_depths_cm]
    rain_vals  = Vector{Union{Float64,Missing}}(undef, ndates)

    for (i, line) in enumerate(data_lines)
        tok = split(line)
        dt[i] = DateTime(tok[1] * " " * tok[2], dateformat"d/m/y H:M")
        for (j, _) in enumerate(temp_depths_cm)
            v = parse(Float64, tok[2 + j])
            temp_vals[j][i] = v < -20.0 ? missing : v
        end
        for (j, _) in enumerate(moist_depths_cm)
            v = parse(Float64, tok[2 + n_temp_cols + j])
            moist_vals[j][i] = v < 0.0 ? missing : v
        end
        rv = parse(Float64, tok[2 + n_temp_cols + n_moist_cols + n_suction_cols + 1])
        rain_vals[i] = rv <= -99.0 ? missing : rv
    end

    obs = DataFrame(DateTime = dt)
    for (j, d) in enumerate(temp_depths_cm)
        obs[!, Symbol("Temp_$(d)cm")] = temp_vals[j]
    end
    for (j, d) in enumerate(moist_depths_cm)
        obs[!, Symbol("SM_$(d)cm")] = moist_vals[j]
    end
    obs[!, :Rainfall] = rain_vals
    return obs
end
