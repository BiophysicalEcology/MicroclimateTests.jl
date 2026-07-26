# Stage 8, HPC split (1/2): solves and caches ONE leg's microclimate --
# historical (member index 0) or one ACCESS-S2 forecast member (1..n_ensembles)
# -- and nothing else. No egg model here; that's points_australia_forecast_eggmodel.jl,
# run afterward once every microclimate job below has completed.
#
# Usage: julia points_australia_forecast_microclimate.jl <member>
#   member = 0        -> solve the historical (SILO) leg
#   member = 1..99     -> solve that ACCESS-S2 forecast member's leg
#
# Intended as a SLURM job array over 0:n_ensembles, with the forecast members
# (1:n_ensembles) depending on member 0 finishing first (they need its
# final soil state as their `init=`).

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", "c:/Spatial_Data/")

include(joinpath(@__DIR__, "points_australia_forecast_setup.jl"))

member = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 0

if member == 0
    # sanity check before committing to the full (expensive) solve below --
    # only the historical (member 0) job runs this, so concurrent array
    # tasks don't race to write the same domain-check png.
    using Plots
    let
        domain = Extent(X=(minimum(lon_range) - 1, maximum(lon_range) + 1),
                         Y=(minimum(lat_range) - 1, maximum(lat_range) + 1))
        elv_domain = crop(CRUCL2_ELV; to=domain, touches=true)
        p = plot(elv_domain; title="CRUCL2 elevation + grid points (green=kept, red=rejected)")
        rejected = setdiff(all_grid_points, points)
        scatter!(p, first.(points), last.(points); markersize=2, markerstrokewidth=0, color=:green, label="kept")
        scatter!(p, first.(rejected), last.(rejected); markersize=2, markerstrokewidth=0, color=:red, label="rejected")
        savefig(p, joinpath(output_dir, "history_forecast_splice_domain_check.png"))
    end

    println("[member 0] Solving historical SILO microclimate for $n points: $oviposition_date to $issue_date...")
    solve_batched(build_historical_model(), historical_label(), points, historical_dates,
        (; soil_moisture=fill(0.2, length(depths))))
    println("[member 0] Historical leg done.")
else
    1 <= member <= n_ensembles || error("member=$member out of range 1:$n_ensembles")

    # Needs the historical leg's final soil state to seed this member's init=
    # -- expected to already be cached (this job should be scheduled with a
    # dependency on member 0), so this is just a fast deserialize, not a re-solve.
    println("[member $member] Loading historical leg for splice continuity...")
    historical_raw = solve_batched(build_historical_model(), historical_label(), points, historical_dates,
        (; soil_moisture=fill(0.2, length(depths))))
    # MicroVectorProblem takes one shared init.soil_moisture/soil_temperature
    # vector for all points, not a point x depth matrix -- average the
    # per-point historical endpoints across points (a known approximation,
    # not per-point continuity) to seed the forecast leg.
    now_soil_temperature = reduce(+, historical_raw.final_soil_temperature) ./ n
    now_soil_moisture = reduce(+, historical_raw.final_soil_moisture) ./ n

    println("[member $member] Solving ACCESS-S2 forecast microclimate for $n points...")
    solve_batched(build_forecast_model(member), forecast_label(member), points, forecast_dates,
        (; soil_moisture=now_soil_moisture, soil_temperature=now_soil_temperature))
    println("[member $member] Forecast leg done.")
end
