# Stage 8, HPC split (2/2): loads the microclimate results already solved and
# cached by points_australia_forecast_microclimate.jl (one job per leg/member
# -- run that first, for member=0 then 1:n_ensembles) and runs the (cheap,
# fast) egg model over all of them, producing the stats CSV and maps.
#
# Safe to also run standalone/serially: solve_batched falls back to solving
# anything not already cached, it just won't get the job-array parallelism.

ENV["RASTERDATASOURCES_PATH"] = get(ENV, "RASTERDATASOURCES_PATH", "c:/Spatial_Data/")

using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Statistics
using NaturalEarth, GeoInterface

include(joinpath(@__DIR__, "points_australia_forecast_setup.jl"))

# no trajectories needed -- outputs are maps + a stats CSV, not per-point
# development/mass/temperature curves.
save_trajectory = false

# ── egg model, identical config to point_silo_deterministic.jl ──

geometry = Ellipsoid(initial_egg_mass, egg_density, axis_ratio, axis_ratio)
arrest = ProportionWindowArrest(;
    cold_temperature, diapause_window, quiescence_windows,
    cold_hour_threshold, diapause_hour_threshold,
    desiccation_tolerance,
)
dm = arrhenius_development_model(;
    T_A, T_AL, T_AH, T_L, T_H, T_ref, rate_at_reference, rate_unit,
)
stage = SteppedHydricStage(;
    conductance_threshold, wetness_threshold, dormant_conductance, active_conductance,
    dormant_wetness, active_wetness,
)
pars = EggParameters(;
    hydraulic_conductance, specific_hydration, conduction_fraction, skin_wetness,
    initial_egg_mass, minimum_egg_mass,
)
survival_model = CombinedSurvival(
    HardTemperatureLimit(; lower_lethal_temperature=u"K"(-5.0u"°C"), upper_lethal_temperature=u"K"(52.0u"°C")),
    DesiccationLimit(; dry_mass=0.1 * initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

# ── load historical microclimate (cached by the member=0 microclimate job) ──

println("Loading historical SILO microclimate for $n points...")
historical_raw = solve_batched(build_historical_model(), historical_label(), points, historical_dates,
    (; soil_moisture=fill(0.2, length(depths))))

historical_day_range = 1:size(historical_raw.per_point[1].soil_temperature, 1)
historical_tspan = (0.0u"hr", length(historical_day_range) * 1.0u"hr")
historical_forcings = [egg_nest_forcing(historical_raw.per_point[i], historical_day_range, nest_node, environment_pars)
                        for i in 1:n]

# ── egg-model cache pool: one integrator per thread, reinit! per work item ──

nworkers = min(Threads.nthreads(), n)
build_cache() = init_egg_cache(egg_model, pars, initial_state, soil_hydraulics,
    historical_forcings[1], historical_tspan; save_trajectory)
cache_pool = Channel{typeof(build_cache())}(nworkers)
for _ in 1:nworkers
    put!(cache_pool, build_cache())
end

function run_threaded!(results, indices, f)
    work = Channel{Int}(length(indices))
    for i in indices
        put!(work, i)
    end
    close(work)
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            cache = take!(cache_pool)
            for i in work
                results[i] = f(cache, i)
            end
            put!(cache_pool, cache)
        end
    end
end

println("Running historical egg model at $n points, lay date $oviposition_date...")
historical_egg_results = Vector{Any}(undef, n)
@time run_threaded!(historical_egg_results, 1:n, (cache, i) ->
    simulate_egg!(cache, initial_state, soil_hydraulics, historical_forcings[i], historical_tspan))

n_hatched_historically = count(r -> r.hatched, historical_egg_results)
n_died_historically = count(r -> r.died, historical_egg_results)
println("$n_hatched_historically/$n hatched and $n_died_historically/$n died before $issue_date; " *
        "$(n - n_hatched_historically - n_died_historically)/$n continue into the forecast.")
forecast_point_indices = [i for i in 1:n
                           if !historical_egg_results[i].hatched && !historical_egg_results[i].died]

# ── forecast leg: load each ACCESS-S2 member's cached microclimate (solved
# by that member's own microclimate job), egg model threaded across points ──

forecast_outcomes = Matrix{Any}(nothing, n, length(members))

if !isempty(forecast_point_indices)
    for (member_index, member) in pairs(members)
        println("\nACCESS-S2 member $member/$(length(members)): loading microclimate for $n points...")
        forecast_raw = solve_batched(build_forecast_model(member), forecast_label(member), points, forecast_dates,
            (; soil_moisture=fill(0.2, length(depths))))   # init unused on a cache hit

        forecast_day_range = 1:size(forecast_raw.per_point[1].soil_temperature, 1)
        forecast_forcings = Dict(i => egg_nest_forcing(forecast_raw.per_point[i], forecast_day_range, nest_node, environment_pars)
                                  for i in forecast_point_indices)

        member_outcomes = Vector{Any}(undef, n)
        run_threaded!(member_outcomes, forecast_point_indices, (cache, i) ->
            simulate_egg!(cache, historical_egg_results[i].final_state, soil_hydraulics,
                forecast_forcings[i], forecast_tspan))
        forecast_outcomes[forecast_point_indices, member_index] .= member_outcomes[forecast_point_indices]
    end
end

# ── per-point summary: median/earliest/latest hatching member, death-cause
# counts, saved to one combined CSV ──

function summarize_point(historical_egg_result, point_outcomes, members)
    if historical_egg_result.hatched
        return (; status=:hatched_historically, counts=Dict(:hatched => 1),
            earliest=nothing, median=nothing, latest=nothing)
    elseif historical_egg_result.died
        return (; status=Symbol("died_historically_$(historical_egg_result.death_cause)"),
            counts=Dict(historical_egg_result.death_cause => 1),
            earliest=nothing, median=nothing, latest=nothing)
    end
    counts = Dict{Symbol,Int}()
    for r in point_outcomes
        code = r.hatched ? :hatched : r.died ? r.death_cause : :timeout
        counts[code] = get(counts, code, 0) + 1
    end
    hatched = [(member, r) for (member, r) in zip(members, point_outcomes) if r.hatched]
    if isempty(hatched)
        return (; status=:no_hatch, counts, earliest=nothing, median=nothing, latest=nothing)
    end
    sorted = sort(hatched; by=x -> x[2].hatch_time)
    (; status=:forecast, counts,
        earliest=sorted[1], median=sorted[cld(length(sorted), 2)], latest=sorted[end])
end

hatch_date_of(member_result, issue_date) =
    issue_date + Day(round(Int, ustrip(u"d", member_result[2].hatch_time)))
mass_mg_of(member_result) = ustrip(u"mg", member_result[2].final_state.egg_mass)

point_summaries = [summarize_point(historical_egg_results[i],
    [forecast_outcomes[i, m] for m in eachindex(members) if forecast_outcomes[i, m] !== nothing], members)
    for i in 1:n]

stats_path = joinpath(output_dir, "history_forecast_splice_stats.csv")
open(stats_path, "w") do io
    println(io, "lon,lat,status,n_members,n_hatched,n_died_cold,n_died_heat,n_died_desiccation,n_timeout," *
                 "earliest_hatch_date,earliest_hatch_mass_mg,median_hatch_date,median_hatch_mass_mg," *
                 "latest_hatch_date,latest_hatch_mass_mg")
    for i in 1:n
        s = point_summaries[i]
        n_run = s.earliest === nothing && s.status != :no_hatch ? 0 : length(members)
        c(sym) = get(s.counts, sym, 0)
        fields = if s.earliest === nothing
            ("", "", "", "", "", "")
        else
            (string(hatch_date_of(s.earliest, issue_date)), string(round(mass_mg_of(s.earliest); digits=2)),
             string(hatch_date_of(s.median, issue_date)), string(round(mass_mg_of(s.median); digits=2)),
             string(hatch_date_of(s.latest, issue_date)), string(round(mass_mg_of(s.latest); digits=2)))
        end
        println(io, join((points[i][1], points[i][2], s.status, n_run, c(:hatched), c(:cold), c(:heat),
            c(:desiccation), c(:timeout), fields...), ","))
    end
end
println("Saved $stats_path")

# ── maps: heatmap over the regular grid + basemap (points_australia.jl) ──

using Plots

point_to_index = Dict(p => i for (i, p) in enumerate(all_grid_points))
to_heatmap_z(vals) = permutedims(reshape(vals, length(lon_range), length(lat_range)))

const GI = GeoInterface
function _collect_boundary!(xs, ys, geom)
    trait = GI.geomtrait(geom)
    if trait isa GI.AbstractPointTrait
        push!(xs, GI.x(geom)); push!(ys, GI.y(geom))
    elseif trait isa GI.LinearRingTrait || trait isa GI.LineStringTrait
        for pt in GI.getpoint(geom)
            push!(xs, GI.x(pt)); push!(ys, GI.y(pt))
        end
        push!(xs, NaN); push!(ys, NaN)
    else
        for sub in GI.getgeom(geom)
            _collect_boundary!(xs, ys, sub)
        end
    end
end
au_states = naturalearth("admin_1_states_provinces", 10)
state_xs, state_ys = Float64[], Float64[]
for i in 1:length(au_states)
    au_states[i].iso_a2 == "AU" || continue
    _collect_boundary!(state_xs, state_ys, au_states[i].geometry)
end

map_towns = ["Birdsville, Queensland", "Roma, Queensland", "Charleville, Queensland",
             "Dubbo, New South Wales", "Broken Hill, New South Wales", "Bourke, New South Wales",
             "Mildura, Victoria", "Canberra, Australia"]
towns_cache_path = joinpath(output_dir, "points_australia_towns.jls")
towns = if isfile(towns_cache_path)
    deserialize(towns_cache_path)
else
    t = map(name -> geocode(name), map_towns)
    serialize(towns_cache_path, t)
    t
end

function add_basemap!(p)
    plot!(p, state_xs, state_ys; color=:black, linewidth=0.75, label=nothing)
    scatter!(p, [t.lon for t in towns], [t.lat for t in towns]; color=:black, markersize=2, label=nothing)
    for t in towns
        annotate!(p, t.lon, t.lat, text(first(split(t.display_name, ",")), 6, :left, :bottom))
    end
    xlims!(p, extrema(lon_range)...)
    ylims!(p, extrema(lat_range)...)
    p
end

historical_hatch_date(r) = oviposition_date + Day(round(Int, ustrip(u"d", r.hatch_time)))

median_dates = fill(NaN, length(all_grid_points))
earliest_dates = fill(NaN, length(all_grid_points))
latest_dates = fill(NaN, length(all_grid_points))
median_mass = fill(NaN, length(all_grid_points))
earliest_mass = fill(NaN, length(all_grid_points))
latest_mass = fill(NaN, length(all_grid_points))
frac_hatched = fill(NaN, length(all_grid_points))
frac_cold = fill(NaN, length(all_grid_points))
frac_heat = fill(NaN, length(all_grid_points))
frac_desiccation = fill(NaN, length(all_grid_points))
for i in 1:n
    gi = point_to_index[points[i]]
    s = point_summaries[i]
    if s.status == :hatched_historically
        d = Dates.value(historical_hatch_date(historical_egg_results[i]))
        median_dates[gi] = earliest_dates[gi] = latest_dates[gi] = d
        m = ustrip(u"mg", historical_egg_results[i].final_state.egg_mass)
        median_mass[gi] = earliest_mass[gi] = latest_mass[gi] = m
    elseif s.earliest !== nothing
        median_dates[gi] = Dates.value(hatch_date_of(s.median, issue_date))
        earliest_dates[gi] = Dates.value(hatch_date_of(s.earliest, issue_date))
        latest_dates[gi] = Dates.value(hatch_date_of(s.latest, issue_date))
        median_mass[gi] = mass_mg_of(s.median)
        earliest_mass[gi] = mass_mg_of(s.earliest)
        latest_mass[gi] = mass_mg_of(s.latest)
    end
    # historically-resolved points are a single outcome (n=1), not an ensemble
    total = s.status == :hatched_historically || startswith(string(s.status), "died_historically") ? 1 : length(members)
    frac_hatched[gi] = get(s.counts, :hatched, 0) / total
    frac_cold[gi] = get(s.counts, :cold, 0) / total
    frac_heat[gi] = get(s.counts, :heat, 0) / total
    frac_desiccation[gi] = get(s.counts, :desiccation, 0) / total
end

# Colourbars can't show date strings (Plots.jl GR limitation) -- build a
# series legend from dummy points instead, same trick as points_australia.jl.
function date_heatmap(title, ordinals)
    valid = filter(!isnan, ordinals)
    isempty(valid) && return add_basemap!(heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude"))
    lo, hi = extrema(valid)
    gradient = cgrad(:plasma)
    p = heatmap(lon_range, lat_range, to_heatmap_z(ordinals);
        title, xlabel="Longitude", ylabel="Latitude", color=gradient, colorbar=false)
    tick_vals = round.(Int, range(lo, hi; length=min(6, length(unique(round.(Int, valid))))))
    for v in unique(tick_vals)
        scatter!(p, [NaN], [NaN]; color=gradient[(v - lo) / max(hi - lo, 1)], markersize=6,
            markerstrokewidth=0, label=string(Date(Dates.UTD(v))))
    end
    add_basemap!(p)
end

function mass_heatmap(title, vals)
    p = heatmap(lon_range, lat_range, to_heatmap_z(vals);
        title, xlabel="Longitude", ylabel="Latitude", color=cgrad(:viridis))
    add_basemap!(p)
end
function frac_heatmap(title, vals)
    p = heatmap(lon_range, lat_range, to_heatmap_z(vals);
        title, xlabel="Longitude", ylabel="Latitude", color=cgrad(:viridis), clims=(0, 1))
    add_basemap!(p)
end

plot_title = "Lay date $oviposition_date, issue date $issue_date, diapause=$diapause"

hatch_date_panel = plot(
    date_heatmap("Median hatch date", median_dates),
    date_heatmap("Earliest hatch date", earliest_dates),
    date_heatmap("Latest hatch date", latest_dates);
    layout=(2, 2), size=(1400, 1050), legendfontsize=6, plot_title,
)
hatch_date_path = joinpath(output_dir, "history_forecast_splice_hatch_dates.png")
savefig(hatch_date_panel, hatch_date_path)
println("Saved $hatch_date_path")

mass_panel = plot(
    mass_heatmap("Median egg mass at hatch (mg)", median_mass),
    mass_heatmap("Earliest egg mass at hatch (mg)", earliest_mass),
    mass_heatmap("Latest egg mass at hatch (mg)", latest_mass);
    layout=(2, 2), size=(1400, 1050), plot_title,
)
mass_path = joinpath(output_dir, "history_forecast_splice_mass.png")
savefig(mass_panel, mass_path)
println("Saved $mass_path")

mortality_panel = plot(
    frac_heatmap("Fraction hatched", frac_hatched),
    frac_heatmap("Fraction died of heat", frac_heat),
    frac_heatmap("Fraction died of cold", frac_cold),
    frac_heatmap("Fraction died of desiccation", frac_desiccation);
    layout=(2, 2), size=(1400, 1050), plot_title,
)
mortality_path = joinpath(output_dir, "history_forecast_splice_mortality.png")
savefig(mortality_panel, mortality_path)
println("Saved $mortality_path")
