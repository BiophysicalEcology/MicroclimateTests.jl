# utils.jl — I/O + interpolation helpers for the micropoint (R) vs
# Microclimate.jl vs real OzFlux obs 3-way comparison. No dependencies on
# comparisons/ozflux -- only report.jl and ozflux_comparison.jl use these.

read_veg(outdir, name) = Matrix{Float64}(coalesce.(Matrix(CSV.read(joinpath(outdir, "veg_$(name).csv"), DataFrame; missingstring="NA")), NaN))

# Reads back every veg_*.csv micropoint wrote (write_ozflux_micropoint_inputs.jl
# -> run_micropoint_ozflux_vegetated.R -> here) plus its own z/depth grids,
# bundled into one NamedTuple for report_site_results3.
function read_micropoint_output(outdir)
    tair, tleaf, relhum, wind =
        read_veg(outdir, "tair"), read_veg(outdir, "tleaf"), read_veg(outdir, "relhum"), read_veg(outdir, "windspeed")
    Rdirdown, Rdifdown, Rswup, Rlwdown, Rlwup =
        read_veg(outdir, "Rdirdown"), read_veg(outdir, "Rdifdown"), read_veg(outdir, "Rswup"),
        read_veg(outdir, "Rlwdown"), read_veg(outdir, "Rlwup")
    tsoil, theta = read_veg(outdir, "tsoil"), read_veg(outdir, "theta")
    layers = CSV.read(joinpath(outdir, "canopy_layers.csv"), DataFrame)
    soil_depths = CSV.read(joinpath(outdir, "soil_depths.csv"), DataFrame)
    return (; tair, tleaf, relhum, wind, Rdirdown, Rdifdown, Rswup, Rlwdown, Rlwup, tsoil, theta,
              z_m = layers.z_m, paii = layers.paii, depth_m = soil_depths.depth_m)
end

# micropoint's own flat grid onto an arbitrary target height/depth -- same
# linear-interpolation style as comparisons/ozflux/report.jl's _profile_series,
# but on a plain evenly-spaced grid rather than Julia's uneven, reversed-column
# one.
function _interp_col(mat, grid_m, target_m)
    j = searchsortedlast(grid_m, target_m)
    if j <= 0
        return mat[:, 1]
    elseif j >= length(grid_m)
        return mat[:, end]
    else
        lo, hi = grid_m[j], grid_m[j + 1]
        t = (target_m - lo) / (hi - lo)
        return mat[:, j] .* (1 - t) .+ mat[:, j + 1] .* t
    end
end
