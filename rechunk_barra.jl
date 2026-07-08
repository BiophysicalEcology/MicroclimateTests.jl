# Rechunk a BARRA NetCDF file so a full time series at one point can be read
# in a single chunk fetch, instead of ~60 small chunk reads per file.
#
# The source files are chunked (168, 130, 12) — 12 timesteps per chunk, so a
# full month (~700+ hours) is split across ~60 chunks. That's efficient for
# "read a spatial map at a few times" but terrible for "read a full time
# series at one point", which is our access pattern. This rewrites the main
# data variable with chunks of (168, 130, ntime) — same spatial tiling, but
# one chunk covers the *entire* time axis per tile. Processes one spatial
# tile at a time (~130 MB), so memory stays bounded regardless of file size.
#
# Works on raw (`.var`) variables throughout — not the CF-translated view —
# since CF-decoded types (e.g. `time` as `DateTime`) aren't valid NetCDF
# on-disk types for `defVar`. Attributes (units, calendar, scale_factor, …)
# are copied through unchanged, so the rechunked file decodes identically.
#
# Usage: julia rechunk_barra.jl <input.nc> <output.nc> <varname>

using NCDatasets

function rechunk_time_major(src_path::String, dst_path::String, varname::String;
                            xy_chunk::Tuple{Int,Int} = (168, 130))
    NCDataset(src_path) do src
        srcvar = src[varname].var
        dims = NCDatasets.dimnames(srcvar)
        sz = size(srcvar)
        length(dims) == 3 || error("expected a 3-D variable, got dims $dims")
        nx, ny, nt = sz
        xchunk, ychunk = min(xy_chunk[1], nx), min(xy_chunk[2], ny)

        NCDataset(dst_path, "c") do dst
            for (dname, dlen) in src.dim
                defDim(dst, dname, dlen)
            end
            for name in keys(src)
                v = src[name].var
                vdims = NCDatasets.dimnames(v)
                fillvalue = get(v.attrib, "_FillValue", nothing)
                # Scalar (dimensionless) variables like `crs` can't be
                # chunked/compressed — only arrays with dims support that.
                dst_v = if isempty(vdims)
                    defVar(dst, name, eltype(v), vdims; fillvalue)
                else
                    chunksizes = name == varname ? [xchunk, ychunk, nt] : nothing
                    defVar(dst, name, eltype(v), vdims;
                        chunksizes, deflatelevel = 1, shuffle = true, fillvalue)
                end
                for (k, val) in v.attrib
                    k == "_FillValue" && continue
                    dst_v.attrib[k] = val
                end
                if name != varname
                    dst_v[:] = v[:]
                end
            end
            for (k, val) in src.attrib
                dst.attrib[k] = val
            end

            dst_v = dst[varname].var
            for xi in 1:xchunk:nx, yi in 1:ychunk:ny
                xr = xi:min(xi + xchunk - 1, nx)
                yr = yi:min(yi + ychunk - 1, ny)
                dst_v[xr, yr, :] = srcvar[xr, yr, :]
            end
        end
    end
    return dst_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    src, dst, varname = ARGS
    rechunk_time_major(src, dst, varname)
    println("wrote ", dst)
end
