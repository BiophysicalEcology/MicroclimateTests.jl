# Pre-download and rechunk SILO annual files for grid-scale point-extraction
# runs (e.g. egg_model/demos/points_australia.jl). SILO's on-disk chunking is
# (nx, ny, 1) -- one chunk per day, covering the entire national grid --
# efficient for reading a spatial map at one time, terrible for reading a
# full time series at many scattered points (each day's chunk gets
# re-decompressed once per point instead of once total). Rechunks in place
# using the same time-major approach already validated for BARRA
# (rechunk_barra.jl's `rechunk_time_major`, reused here unchanged).
#
# Usage: julia --project=<env> prefetch_silo.jl <year> [<year> ...]

using RasterDataSources
using NCDatasets
using Dates

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", "z:/")

include(joinpath(@__DIR__, "rechunk_barra.jl"))  # reuses rechunk_time_major

const SILO_LAYERS_TO_PREFETCH = (:max_temp, :min_temp, :daily_rain, :rh_tmax, :rh_tmin, :radiation)

# Heuristic: original SILO files chunk time as 1 (one day per chunk);
# rechunk_time_major always writes the *entire* time axis as one chunk, so
# a last-dim chunk size >1 means this file has already been fixed.
function _is_time_major_chunked(path, varname)
    NCDataset(path) do ds
        chunking = NCDatasets.chunking(ds[varname].var)
        chunking isa Tuple && last(chunking[2]) > 1
    end
end

function prefetch_and_rechunk_silo(years; layers = SILO_LAYERS_TO_PREFETCH)
    for layer in layers, yr in years
        println("Fetching SILO $layer $yr...")
        path = getraster(SILO, layer; date = Date(yr, 1, 1))
        if _is_time_major_chunked(path, String(layer))
            println("  already time-major chunked, skipping rechunk.")
            continue
        end
        tmp_path = path * ".rechunk_tmp"
        println("  rechunking $path...")
        @time rechunk_time_major(path, tmp_path, String(layer))
        mv(tmp_path, path; force = true)
        println("  done.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    years = parse.(Int, ARGS)
    isempty(years) && error("Usage: julia prefetch_silo.jl <year> [<year> ...]")
    prefetch_and_rechunk_silo(years)
end
