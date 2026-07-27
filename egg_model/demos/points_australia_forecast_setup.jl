# Shared setup for the split microclimate/egg-model HPC scripts
# (points_australia_forecast_microclimate.jl, points_australia_forecast_eggmodel.jl).
# Grid, land masks, soil texture, and microclimate batch-caching -- everything
# both scripts need to agree on to index into the same cache files.

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_properties_model, example_soil_hydraulic_model
using Rasters, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using Dates, Unitful
using Serialization

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))
include(joinpath(@__DIR__, "..", "src", "access_s2.jl"))

include(joinpath(@__DIR__, "..", "params", "chortoicetes.jl"))

output_dir = get(ENV, "LOCUST_FORECAST_OUTPUT_DIR", joinpath(@__DIR__, "output"))
mkpath(output_dir)

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

n_ensembles = 25    # how many of the up-to-99 ACCESS-S2 members to run
diapause = true
oviposition_date = Date(2026, 4, 25)
use_cache = true
batch_size = 100

if diapause
    nest_depth = 5.0u"cm"
    cold_hour_threshold = 30u"d"
    diapause_hour_threshold = 240.0u"d"
else
    nest_depth = 10.0u"cm"
    cold_hour_threshold = 0.0u"d"
    diapause_hour_threshold = 0.0u"d"
end

# raw caches only keep depths down to this value -- set past whatever
# nest_depth values you want to compare (deeper needs a fresh solve).
cache_max_depth = 15.0u"cm"
CACHED_DEPTH_RANGE = 1:nearest_node(cache_max_depth, depths)
nest_node = nearest_node(nest_depth, depths)
nest_node in CACHED_DEPTH_RANGE || error(
    "nest_depth=$nest_depth is deeper than cache_max_depth=$cache_max_depth -- raise cache_max_depth and re-solve")

# Just what egg_nest_forcing (forcing.jl) actually reads, plus soil_moisture
# (kept for a possible future moisture threshold).
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
)

issue_date = Date(2026, 7, 1)   # ACCESS-S2 issue date -- also used below for the ACCESS-S2 land-mask probe

# ── regular grid over south-eastern Australia, land-masked (points_australia.jl) ──

lon_range = range(135.0, 153.0; length=50)
lat_range = range(-39.0, -24.0; length=40)
all_grid_points = vec([(lon, lat) for lon in lon_range, lat in lat_range])   # column-major (lon fastest)

# TODO work out what DEM SILO uses and use that instead, probably only need
# that and not the CRUCL2_ELV filter.
const CRUCL2_ELV = read(Raster(RasterDataSources.getraster(CRUCL2); name=:elv, lazy=true))
has_crucl2_land(lon, lat) = !ismissing(CRUCL2_ELV[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_crucl2_land(p...), all_grid_points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the CRUCL2 land-mask pre-check.")

const SILO_MAXTEMP_PROBE = read(Raster(RasterDataSources.getraster(SILO, :max_temp; date=Date(2025, 1, 1)); name=:max_temp, lazy=true)[Ti(1)])
has_silo_land(lon, lat) = !ismissing(SILO_MAXTEMP_PROBE[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_silo_land(p...), points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the SILO land-mask pre-check.")

# ACCESS-S2's own land/sea representation doesn't always agree with
# CRUCL2/SILO's -- a point that passes both those masks can still be
# `missing` in ACCESS-S2 at its (coarser) resolution, which propagates into
# the soil ODE as degenerate forcing rather than a clean error.
const ACCESS_S2_TMAX_PROBE = read(Raster(joinpath(ENV["RASTERDATASOURCES_PATH"], "ACCESS-S2", "$(Dates.format(issue_date, "yyyymmdd"))_tmax.nc");
    name=:tmax, lazy=true)[Dim{:ensemble}(Rasters.At(1)), Ti(1)])
has_access_s2_land(lon, lat) = !ismissing(ACCESS_S2_TMAX_PROBE[X(Near(lon)), Y(Near(lat))])
points = filter(p -> has_access_s2_land(p...), points)
println("$(length(points))/$(length(all_grid_points)) grid points kept after the ACCESS-S2 land-mask pre-check.")

# ── uniform soil texture (Campbell & Norman) -- Pass A, no per-point texture ──

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
soil_source = :sandy_loam

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
            campbell_b_parameter=fill(texture.b, n), root_density,
        ),
    )
end
soil_profile = soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[soil_source], depths)
soil_hydraulics = (;
    air_entry_potential    = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b             = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)
environment_pars = example_environment_pars()

historical_dates = oviposition_date:Day(1):issue_date
members = 1:n_ensembles
n = length(points)

forecast_horizon_days = 214
forecast_dates = issue_date:Day(1):(issue_date + Day(forecast_horizon_days - 1))
forecast_tspan = (0.0u"hr", length(forecast_dates) * 24.0u"hr")

# ── microclimate model configs ──

build_historical_model() = MicroMapModel(;
    micro_model=MicroModel(;
        depths, heights,
        soil_properties_model=example_soil_properties_model(),
        soil_hydraulic_model=example_soil_hydraulic_model(),
        snow_model=NoSnow(),
        config=MicroConfig(soil_moisture_strategy=DynamicSoilMoisture()),
    ),
    dem_source=CRUCL2, weather_source=SILO,
    surface_albedo_source=0.15, roughness_height_source=0.004u"m",
    compute_terrain=false, output_layers,
)

build_forecast_model(member) = MicroMapModel(;
    micro_model=MicroModel(;
        depths, heights,
        soil_properties_model=example_soil_properties_model(),
        soil_hydraulic_model=example_soil_hydraulic_model(),
        snow_model=NoSnow(),
        config=MicroConfig(
            soil_moisture_strategy=DynamicSoilMoisture(),
            convergence=FixedSoilTemperatureIterations(1),
        ),
    ),
    dem_source=CRUCL2, weather_source=AccessS2(issue_date, member),
    surface_albedo_source=0.15, roughness_height_source=0.004u"m",
    compute_terrain=false, output_layers,
)

historical_label() = "splice_historical_lay$(oviposition_date)_issue$(issue_date)"
forecast_label(member) = "splice_forecast_member$(member)_issue$(issue_date)"

_raw_layers(output, i) = (;
    soil_temperature          = Float32.(collect(output.soil_temperature[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_moisture             = Float32.(collect(output.soil_moisture[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_water_potential      = Float32.(collect(output.soil_water_potential[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_thermal_conductivity = Float32.(collect(output.soil_thermal_conductivity[point=i, depth=CACHED_DEPTH_RANGE])),
    soil_humidity             = Float32.(collect(output.soil_humidity[point=i, depth=CACHED_DEPTH_RANGE])),
)

# Solves `model` over `points` in batches, caching each batch to disk under
# `label`. On a cache hit (the expected case for every call except the one
# job that actually solves a given leg/member), this just deserializes --
# `model`/`dates`/`init` are cheap to reconstruct regardless of whether
# they'll actually be used, so every script can safely call this the same
# way without needing to know whether the underlying solve already ran.
function solve_batched(model, label, points, dates, init)
    n = length(points)
    n_batches = cld(n, batch_size)
    per_point = Vector{Any}(undef, n)
    final_soil_temperature = Vector{Any}(undef, n)
    final_soil_moisture = Vector{Any}(undef, n)
    for b in 1:n_batches
        i_start, i_end = (b - 1) * batch_size + 1, min(b * batch_size, n)
        batch_points = points[i_start:i_end]
        batch_cache_path = joinpath(output_dir, "$(label)_batch$(b)of$(n_batches)_n$(n).jls")
        batch = if isfile(batch_cache_path) && use_cache
            println("Loading cached $label batch $b/$n_batches ($(length(batch_points)) points)...")
            deserialize(batch_cache_path)
        else
            println("Solving $label batch $b/$n_batches ($(length(batch_points)) points)...")
            batch_problem = MicroVectorProblem(; model, points=batch_points, dates, soil_profile, init)
            @time batch_output = solve(batch_problem)
            # batch_output's raw axis order is (point, Ti, depth) -- size(...,1)
            # silently gives the point count, not the hour count; size(...,Ti)
            # looks the dimension up by name instead of position.
            last_hour = size(batch_output.soil_temperature, Ti)
            result = (;
                per_point = [_raw_layers(batch_output, i) for i in 1:length(batch_points)],
                # full (untrimmed) depth profile at the last hour, for splice continuity --
                # kept Float64 (unlike per_point's Float32 trim): this seeds
                # MicroVectorProblem's init=, which expects Float64 internally.
                final_soil_temperature = [collect(batch_output.soil_temperature[point=i, Ti=last_hour]) for i in 1:length(batch_points)],
                final_soil_moisture    = [collect(batch_output.soil_moisture[point=i, Ti=last_hour]) for i in 1:length(batch_points)],
            )
            serialize(batch_cache_path, result)
            batch_output = nothing
            batch_problem = nothing
            GC.gc()
            result
        end
        per_point[i_start:i_end] .= batch.per_point
        final_soil_temperature[i_start:i_end] .= batch.final_soil_temperature
        final_soil_moisture[i_start:i_end] .= batch.final_soil_moisture
    end
    (; per_point, final_soil_temperature, final_soil_moisture)
end
