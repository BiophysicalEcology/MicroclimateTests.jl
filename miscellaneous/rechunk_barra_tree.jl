# Walks a BARRA archive and rechunks every time-dimensioned NetCDF file so
# each spatial tile holds one chunk spanning the *entire* time axis, instead
# of the native (168, 130, 12) shape (~60 small chunks per file for an hourly
# month). Static files with no time dimension (e.g. fx/orog) are left alone.
#
# Run this ON THE SERVER, against local disk paths, with plenty of RAM
# available — each file is read/written as one bulk in-memory array (see
# below), not tiled, so needs ~1-2 GB free per file in flight.
#
# Overwrites files in place: each file is rechunked to a `.rechunk_tmp`
# sibling (never opening the same path for both read and write — `"c"` mode
# truncates on open, so writing directly over `src_path` would destroy it
# before the read completes), then atomically renamed over the original.
# Only one extra file's worth of disk space is needed at a time, not a full
# second copy of the archive.
#
# Resumable: a file already in the target chunk shape (chunk spans the full
# time axis) is detected from its on-disk chunking and skipped, so an
# interrupted run can just be restarted.
#
# Variable name is inferred from each file's parent directory, matching the
# archive's own layout (.../1hr/psl/*.nc -> "psl", .../day/mrsol/*.nc ->
# "mrsol", ...).
#
# Processes files sequentially, one at a time. Concurrent *reads* of
# different NetCDF files from multiple threads are safe (verified directly
# against BARRA files), but concurrent *writes*/creates are not — HDF5's
# internal state isn't thread-safe across simultaneous file-creation calls,
# even to distinct files, and corrupts other threads' in-flight tmp files.
# That's fine here: this runs against local disk, not the network share, so
# there's no latency to hide with concurrency in the first place.
#
# Usage: julia rechunk_barra_tree.jl <root>

using NCDatasets

function rechunk_time_major(src_path::String, dst_path::String, varname::String;
                            xy_chunk::Tuple{Int,Int} = (168, 130))
    NCDataset(src_path) do src
        srcvar = src[varname].var
        vdims = NCDatasets.dimnames(srcvar)
        sz = size(srcvar)
        # Time is always the last dimension; x/y (lon/lat) the first two.
        # Anything in between (e.g. mrsol/tsl's 4-layer `depth`) is a small
        # extra axis that isn't itself chunked — one chunk just covers it
        # in full, same as it already does for the full time axis.
        length(vdims) >= 3 || error("expected at least a 3-D variable, got dims $vdims")
        nx, ny, nt = sz[1], sz[2], sz[end]
        xchunk, ychunk = min(xy_chunk[1], nx), min(xy_chunk[2], ny)
        chunksizes_full = (xchunk, ychunk, sz[3:end-1]..., nt)
        idx = ntuple(_ -> Colon(), length(sz))

        # Bulk single-shot read/write, not a manual per-tile loop: letting
        # the library iterate its own native chunks internally is ~3x
        # faster than issuing hundreds of small Julia-level chunk reads
        # (verified: 139s vs 444s for a full BARRA grid file). Needs the
        # whole array in memory at once (~1-2 GB decompressed per file).
        data = srcvar[idx...]

        NCDataset(dst_path, "c") do dst
            for (dname, dlen) in src.dim
                defDim(dst, dname, dlen)
            end
            for name in keys(src)
                v = src[name].var
                vd = NCDatasets.dimnames(v)
                fillvalue = get(v.attrib, "_FillValue", nothing)
                dst_v = if isempty(vd)
                    defVar(dst, name, eltype(v), vd; fillvalue)
                else
                    chunksizes = name == varname ? collect(chunksizes_full) : nothing
                    defVar(dst, name, eltype(v), vd;
                        chunksizes, deflatelevel = 1, shuffle = true, fillvalue)
                end
                for (k, val) in v.attrib
                    k == "_FillValue" && continue
                    dst_v.attrib[k] = val
                end
                if name == varname
                    dst_v[idx...] = data
                else
                    dst_v[:] = v[:]
                end
            end
            for (k, val) in src.attrib
                dst.attrib[k] = val
            end
        end
    end
    return dst_path
end

function has_time_dim(path::String, varname::String)
    NCDataset(path) do ds
        "time" in NCDatasets.dimnames(ds[varname].var)
    end
end

function already_rechunked(path::String, varname::String)
    NCDataset(path) do ds
        v = ds[varname].var
        nt = size(v)[end]
        storage, chunksizes = NCDatasets.chunking(v)
        storage == :chunked && chunksizes[end] == nt
    end
end

function rechunk_tree(root::String)
    files = String[]
    for (dir, _, fnames) in walkdir(root)
        for f in fnames
            endswith(f, ".nc") && push!(files, joinpath(dir, f))
        end
    end
    @info "found $(length(files)) NetCDF files"

    for path in files
        varname = basename(dirname(path))
        if !has_time_dim(path, varname)
            continue
        end

        if already_rechunked(path, varname)
            @info "skip (already rechunked): $path"
            continue
        end

        @info "rechunking: $path (var=$varname)..."
        t0 = time()
        tmp_path = path * ".rechunk_tmp"
        rechunk_time_major(path, tmp_path, varname)
        mv(tmp_path, path; force = true)
        @info "  done: $path ($(round(time() - t0, digits = 1))s)"
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    rechunk_tree(ARGS[1])
    println("done.")
end
