# Stage 6 Pass A: broad first-look egg-development sweep over a regular grid
# of points across south-eastern Australia, via MicroVectorProblem (SILO),
# assuming one uniform/optimal Campbell & Norman soil texture rather than
# fetching real per-point texture (that's Pass B, raster-mode, once this
# broad pass identifies a region worth zooming into).
#
# Reuses point_silo_deterministic.jl's exact model/egg-model setup, just with
# many points instead of one, and reuses egg_model/src/phases.jl's
# init_egg_cache/simulate_egg! cache-reuse API (one integrator cache per
# thread, reinit! per point) instead of building a fresh integrator per point
# -- the same Channel-based worker-pool pattern MicroclimateMapper.jl's own
# grid solve uses (src/common.jl).

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_properties_model, example_soil_hydraulic_model
using ThermalPhysiology
using BiophysicalGeometry
using HeatExchange
using BiophysicalBehaviour: EnvironmentForcing
using Rasters, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using DataInterpolations
using Dates, Unitful
using Serialization

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))

# SILO weather: either the live DataDrill point-query API (network, one HTTP
# call per point per variable -- ~12000 calls for a 2000-point sweep, and
# fails hard on the whole batch if even one point/response is bad, as seen
# with the ocean points above) or SILO's regular grid-mode loader, reading
# pre-downloaded annual NetCDF files (one file per variable per year, cropped
# and per-point extracted locally -- far fewer, more robust network calls).
# The user has pre-downloaded SILO (and CRUCL2) under z:/; missing years are
# fetched on demand from the same z:/ path.
use_local_silo = true
if use_local_silo
    ENV["RASTERDATASOURCES_PATH"] = "z:/"
else
    ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"
    ENV["SILO_EMAIL"] = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")
    MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()
end

# ── microclimate: SILO points run, same setup as point_silo_deterministic.jl ──

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
nest_depth = 5.0u"cm"

output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_heat_capacity, :soil),
    LayerSpec(:soil_bulk_density, :soil),
    LayerSpec(:soil_humidity, :soil),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:diffuse_fraction, :scalar),
    LayerSpec(:reference_temperature, :scalar),
    LayerSpec(:pressure, :scalar),
    LayerSpec(:zenith_angle, :solar),
)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
        config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
    ),
    # CRUCL2's :elv band (~18 km, monthly-climate grid) instead of SRTM: SRTM at
    # native resolution across the whole south-eastern Australia bounding box
    # (135-153E, 39-24S) is too large to materialise in memory (OOM: "required
    # memory ... greater than system memory"), and compute_terrain=false means
    # only a flat elevation reference is needed here anyway, not real terrain.
    dem_source              = CRUCL2,
    weather_source          = SILO,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
    output_layers,
)

# regular grid over south-eastern Australia (broad plague-locust range: NSW,
# VIC, SA, southern QLD), thinned to a manageable point count. n_lon*n_lat
# below (not a fixed total) so the grid stays reshape-able into a Raster
# directly -- no need for Rasters.rasterize since this is already regular.
lon_range = range(135.0, 153.0; length=50)
lat_range = range(-39.0, -24.0; length=40)
all_grid_points = vec([(lon, lat) for lon in lon_range, lat in lat_range])   # 2000 points, column-major (lon fastest)

# Reject points that fall in the ocean before handing anything to
# MicroVectorProblem/SILO: MicroVectorProblem fetches one shared bounding-box
# tile set for the whole batch of points, and SILO's DataDrill point API only
# covers the Australian mainland, so even one oceanic point crashes the DEM
# fetch or the SILO point query for *every* point in the run, not just that
# one -- flagged as something to fix upstream in MicroclimateMapper.jl
# (per-point masking, matching how MicroRasterProblem already skips
# masked/oceanic grid cells), worked around here for now since that's out of
# scope to touch right now.
#
# An earlier version of this check used RasterDataSources' 5°x5° SRTM tile
# inventory (`HAS_SRTM_TILE`), but that's too coarse -- a whole tile counts as
# "land" even when a specific point in it is offshore, and it let 0 of the
# real ocean points here get filtered (confirmed: kept 2000/2000). Using
# CRUCL2's `:elv` layer instead (missing over ocean, ~18 km resolution --
# already being fetched anyway now that dem_source=CRUCL2, see above) is a
# much finer land mask: it correctly flags 355/2000 of this grid's points as
# ocean, which is what was actually crashing SILO's point query.
const CRUCL2_ELV = read(Raster(RasterDataSources.getraster(CRUCL2); name=:elv, lazy=true))
has_crucl2_land(lon, lat) = !ismissing(CRUCL2_ELV[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_crucl2_land(p...), all_grid_points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the CRUCL2 land-mask pre-check.")

# Sanity-check plot: CRUCL2 elevation over the domain, with kept (land) points
# in green and rejected (ocean) points in red, so a bad mask is visible before
# committing to a multi-hour SILO sweep.
using Plots
let
    domain = Extent(X = (minimum(lon_range) - 1, maximum(lon_range) + 1),
                     Y = (minimum(lat_range) - 1, maximum(lat_range) + 1))
    elv_domain = crop(CRUCL2_ELV; to = domain, touches = true)
    p = plot(elv_domain; title = "CRUCL2 elevation + grid points (green=kept, red=rejected)")
    rejected = setdiff(all_grid_points, points)
    scatter!(p, first.(points), last.(points); markersize = 2, markerstrokewidth = 0, color = :green, label = "kept")
    scatter!(p, first.(rejected), last.(rejected); markersize = 2, markerstrokewidth = 0, color = :red, label = "rejected")
    savefig(p, joinpath(@__DIR__, "points_australia_domain_check.png"))
    display(p)
end

# cap how many of the (already land-filtered) points actually get
# queried/solved -- start small to verify correctness cheaply (SILO point
# queries are network-bound, not CPU-bound, so a couple thousand of them is a
# genuinely long first real run). Set test_mode=false for the real Pass A
# sweep over the full (filtered) grid once this checks out.
test_mode = false
run_points = if test_mode
    # a small 3x3 cluster around Bendigo (definitely on land, already
    # validated at point scale) -- a safe, quick correctness check, unlike
    # slicing the first few points of the full grid (tried that: it landed in
    # the Southern Ocean south of SA, where SRTM has no tiles at all).
    bendigo_lon, bendigo_lat = 144.2826718, -36.7590183
    vec([(bendigo_lon + dlon, bendigo_lat + dlat) for dlon in -0.2:0.2:0.2, dlat in -0.2:0.2:0.2])
else
    points
end
n_points_to_run = length(run_points)

# whole calendar year(s) only, per current SILO/point-query constraints.
dates = Date(2010, 1, 1):Day(1):Date(2011, 12, 31)
use_cache = true

# ── uniform soil texture (Pass A assumption -- see module docstring above) ──

const CAMPBELL_NORMAN_TEXTURES = (
    sand             = (air_entry=0.7u"J/kg", b=1.7, Ksat=5.8e-3u"kg*s/m^3", field_capacity=0.09, wilting_point=0.03),
    loamy_sand       = (air_entry=0.9u"J/kg", b=2.1, Ksat=1.7e-3u"kg*s/m^3", field_capacity=0.13, wilting_point=0.06),
    sandy_loam       = (air_entry=1.5u"J/kg", b=3.1, Ksat=7.2e-4u"kg*s/m^3", field_capacity=0.21, wilting_point=0.10),
    loam             = (air_entry=1.1u"J/kg", b=4.5, Ksat=3.7e-4u"kg*s/m^3", field_capacity=0.27, wilting_point=0.12),
    silt_loam        = (air_entry=2.1u"J/kg", b=4.7, Ksat=1.9e-4u"kg*s/m^3", field_capacity=0.33, wilting_point=0.13),
    sandy_clay_loam  = (air_entry=2.8u"J/kg", b=4.0, Ksat=1.2e-3u"kg*s/m^3", field_capacity=0.26, wilting_point=0.15),
    clay_loam        = (air_entry=2.6u"J/kg", b=5.2, Ksat=6.4e-5u"kg*s/m^3", field_capacity=0.32, wilting_point=0.20),
    silty_clay_loam  = (air_entry=3.3u"J/kg", b=6.6, Ksat=4.2e-5u"kg*s/m^3", field_capacity=0.37, wilting_point=0.32),
    sandy_clay       = (air_entry=2.9u"J/kg", b=6.0, Ksat=3.3e-5u"kg*s/m^3", field_capacity=0.34, wilting_point=0.24),
    silty_clay       = (air_entry=3.4u"J/kg", b=7.9, Ksat=2.5e-5u"kg*s/m^3", field_capacity=0.39, wilting_point=0.25),
    clay             = (air_entry=3.7u"J/kg", b=7.6, Ksat=1.7e-5u"kg*s/m^3", field_capacity=0.40, wilting_point=0.27),
)

function soil_profile_from_texture(texture::NamedTuple, depths;
    bulk_density=1.3u"Mg/m^3", mineral_density=2.560u"Mg/m^3",
    mineral_conductivity=1.25u"W/m/K", mineral_heat_capacity=870.0u"J/kg/K",
    root_density=Microclimate.example_campbell_hydraulic_profile(depths).root_density,
)
    n = length(depths)
    Microclimate.SoilProfile(;
        bulk_density=fill(bulk_density, n), mineral_density=fill(mineral_density, n),
        mineral_conductivity=fill(mineral_conductivity, n), mineral_heat_capacity=fill(mineral_heat_capacity, n),
        hydraulics=Microclimate.CampbellHydraulicProfile(;
            air_entry_water_potential=fill(-texture.air_entry, n),
            saturated_hydraulic_conductivity=fill(texture.Ksat, n),
            campbell_b_parameter=fill(texture.b, n),
            root_density,
        ),
    )
end

soil_source = :sandy_loam
soil_profile = soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[soil_source], depths)

problem = MicroVectorProblem(;
    model, points=run_points, dates, soil_profile,
    init = (; soil_moisture = fill(0.2, length(depths))),
)

cache_path = joinpath(@__DIR__, "points_australia_cache_n$(n_points_to_run).jls")
if isfile(cache_path) && use_cache
    println("Loading cached microclimate result from $cache_path...")
    output = deserialize(cache_path)
else
    println("Solving SILO microclimate at $(length(run_points)) points...")
    @time output = solve(problem)
    serialize(cache_path, output)
end

nest_node = nearest_node(nest_depth, depths)
environment_pars = example_environment_pars()

# result arrays are hourly (length(dates)*24 rows), not daily -- day_range
# must span the full hourly series or the forcing silently extrapolates a
# constant past the first length(dates) *hours*, not the full run (see
# point_silo_deterministic.jl's identical comment). Computed from one
# point-extracted probe (point=1), NOT the raw (still point-dimensioned)
# `output` directly -- output.soil_temperature has an extra `point` axis, so
# size(output.soil_temperature, 1) measures the wrong dimension entirely.
# Same for every point since they all share the same `dates`.
day_range = 1:size(output.soil_temperature[point=1], 1)

# one forcing closure per point, built once up front (cheap -- interpolators
# only, no heavy computation) and reused across every lay date if more than
# one is tried.
println("Building per-point forcing...")
forcings = map(1:length(run_points)) do i
    result_i = (;
        soil_temperature          = collect(output.soil_temperature[point=i]),
        soil_moisture             = collect(output.soil_moisture[point=i]),
        soil_water_potential      = collect(output.soil_water_potential[point=i]),
        soil_thermal_conductivity = collect(output.soil_thermal_conductivity[point=i]),
        soil_heat_capacity        = collect(output.soil_heat_capacity[point=i]),
        soil_bulk_density         = collect(output.soil_bulk_density[point=i]),
        soil_humidity             = collect(output.soil_humidity[point=i]),
        global_radiation          = collect(output.global_radiation[point=i]),
        sky_temperature           = collect(output.sky_temperature[point=i]),
        diffuse_fraction          = collect(output.diffuse_fraction[point=i]),
        reference_temperature     = collect(output.reference_temperature[point=i]),
        pressure                  = collect(output.pressure[point=i]),
        solar_radiation           = (; zenith_angle = collect(output.zenith_angle[point=i])),
    )
    egg_nest_forcing(result_i, day_range, nest_node, environment_pars)
end

# soil_hydraulics is the SAME for every point (Pass A's uniform-soil
# assumption) -- unlike forcing, no per-point array needed.
soil_hydraulics = (;
    air_entry_potential    = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b             = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)

# ── egg model, identical config to point_silo_deterministic.jl ──

shape_b = 0.69 / 1.82
geometry = Ellipsoid(0.0036u"g", 1000.0u"kg/m^3", 1 / shape_b, 1 / shape_b)
arrest = ProportionWindowArrest(;
    cold_temperature=u"K"(0.0u"°C"), diapause_window=(0.25, 0.30),
    quiescence_windows=((0.25, 0.30), (0.45, 0.50)),
    cold_hour_threshold=1000.0u"hr", diapause_hour_threshold=240.0u"hr",
    desiccation_tolerance=0.6,
)
dm = arrhenius_development_model(;
    T_A=6641.6175, T_AL=33600.0, T_AH=48000.0, T_L=289.15, T_H=314.65, T_ref=301.65,
    rate_at_reference=1 / 17.4, rate_unit=1.0u"d^-1",
)
base_K_e = 2.347802e-9 * u"kg/m^2/s/(J/kg)"
stage = SteppedHydricStage(;
    conductance_threshold=0.25, wetness_threshold=0.45,
    dormant_conductance=0.0u"kg/m^2/s/(J/kg)", active_conductance=base_K_e * 3,
    dormant_wetness=0.35 / 100, active_wetness=0.35,
)
pars = EggParameters(;
    hydraulic_conductance=base_K_e, specific_hydration=0.000304u"m^3/m^3/(J/kg)",
    conduction_fraction=0.5, skin_wetness=0.35 / 100,
    initial_egg_mass=0.0036u"g", minimum_egg_mass=0.0026u"g",
)
survival_model = CombinedSurvival(
    HardTemperatureLimit(; lower_lethal_temperature=u"K"(-5.0u"°C"), upper_lethal_temperature=u"K"(52.0u"°C")),
    DesiccationLimit(; dry_mass=0.1 * pars.initial_egg_mass, critical_water_ratio=0.6),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)

max_duration = 720.0u"d"
forcing_end_hr = length(day_range) * 1.0u"hr"
oviposition_date = Date(2020, 9, 1)
start_hr = oviposition_offset(oviposition_date, dates)
tspan = (start_hr, min(start_hr + max_duration, forcing_end_hr))
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)

# ── grid loop: one egg-model integrator cache per thread, reused across all
# points via reinit! (init_egg_cache/simulate_egg!, phases.jl) -- the same
# per-thread-not-per-cell pattern MicroclimateMapper.jl's own grid solve uses
# (src/common.jl's Channel-based worker pool).

println("Running egg model at $(length(run_points)) points, lay date $oviposition_date...")
n = length(run_points)
# in test_mode, cap workers well below n so cache reuse across *different*
# locations is actually exercised (with nworkers>=n every point would get its
# own dedicated cache, never reinit!'d for a second location) -- the real
# correctness property Stage 6 depends on.
nworkers = test_mode ? min(3, n) : min(Threads.nthreads(), n)

build_cache() = init_egg_cache(egg_model, pars, initial_state, soil_hydraulics, forcings[1], tspan)
cache_pool = Channel{typeof(build_cache())}(nworkers)
put!(cache_pool, build_cache())
for _ in 2:nworkers
    put!(cache_pool, build_cache())
end

work = Channel{Int}(n)
for i in 1:n
    put!(work, i)
end
close(work)

results = Vector{Any}(undef, n)
@time @sync for _ in 1:nworkers
    Threads.@spawn begin
        cache = take!(cache_pool)
        for i in work
            results[i] = simulate_egg!(cache, initial_state, soil_hydraulics, forcings[i], tspan)
        end
        put!(cache_pool, cache)
    end
end

for (i, r) in enumerate(results)
    lon, lat = run_points[i]
    outcome = if r.hatched
        hatch_date = first(dates) + Day(round(Int, ustrip(u"d", r.hatch_time)))
        "hatched on $hatch_date"
    elseif r.died
        "died of $(r.death_cause)"
    else
        "did not hatch in $(max_duration)"
    end
    println("  ($lon, $lat) -> $outcome")
end

# ── rasterize for plotting: the point grid is already regular (lon fastest,
# see `points` construction above), so a plain reshape into a Raster is all
# that's needed here -- no Rasters.rasterize (that's for genuinely irregular
# point sets, e.g. if this ever switches to a random scatter instead).
if n_points_to_run == length(points)
    using Plots
    hatch_days = map(results) do r
        r.hatched ? ustrip(u"d", r.hatch_time) : NaN
    end
    hatch_raster = Raster(reshape(hatch_days, length(lon_range), length(lat_range)), (X(lon_range), Y(lat_range)))
    hatch_plot = plot(hatch_raster; title="Hatch time (hours since $(first(dates))), lay date $oviposition_date")
    savefig(hatch_plot, joinpath(@__DIR__, "points_australia_hatch.png"))
    display(hatch_plot)
else
    println("\nSkipping the hatch-date raster plot -- only $(n_points_to_run)/$(length(points)) points were run " *
            "(set test_mode = false for the full Pass A sweep).")
end
