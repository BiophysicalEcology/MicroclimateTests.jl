using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_hydraulic_model
using Rasters, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using Dates, Statistics, Unitful, Plots

ENV["RASTERDATASOURCES_PATH"] = "Z:"

depths = ([0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0] ./ 100.0) .* u"m"
heights = [0.01, 1.2]u"m"

bulk_density          = 1.3u"Mg/m^3"
mineral_density       = 2.56u"Mg/m^3"
saturation_moisture   = 1 - bulk_density/mineral_density # for reference; not a settable field
mineral_conductivity  = [0.2, 0.2, 0.2, 1.35, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5,
                          2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]u"W/m/K"
mineral_heat_capacity = [1920.0, 1920.0, 1920.0, 1395.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0,
                          870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0]u"J/kg/K"
air_entry_potential   = 1.1u"J/kg"
sat_hydraulic_cond    = 0.0037u"kg*s/m^3"
campbell_b            = 4.5
root_density          = [0.0, 0.0, 82000.0, 80000.0, 78000.0, 74000.0, 71000.0, 64000.0, 58000.0, 48000.0,
                          40000.0, 18000.0, 9000.0, 6000.0, 8000.0, 4000.0, 4000.0, 0.0, 0.0]u"m/m^3"

# Alice Springs — central Australia, hot arid climate.
points = [geocode("Alice Springs, Australia")]

# SILO is a real historical daily record, not a climatology — pick a real year.
dates = Date(2020, 1, 1):Day(1):Date(2020, 12, 31)

# Set true to fetch SILO via direct OPeNDAP point queries (PointDataSources.jl)
# instead of downloading/cropping the whole SILO grid -- an optional
# alternative to the default grid-crop path, not required. Wind still comes
# from CRUCL2's monthly climatology (no point-query source for it), fetched
# via a small area crop and expanded onto the daily Ti axis automatically.
# Requires ENV["SILO_EMAIL"] -- SILO's API requires an email address.
use_point_query = true
if use_point_query
    ENV["SILO_EMAIL"] = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")
    MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()
end

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights,
        soil_properties_model = CampbelldeVriesSoilProperties(;
            de_vries_shape_factor = 0.1,
            recirculation_power = 4.0,
            return_flow_threshold = 0.162,
        ),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
        config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
    ),
    dem_source              = SRTM,
    weather_source          = SILO,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
)

soil_profile = SoilProfile(;
    bulk_density = fill(bulk_density, length(depths)),
    mineral_density = fill(mineral_density, length(depths)),
    mineral_conductivity,
    mineral_heat_capacity,
    hydraulics = CampbellHydraulicProfile(;
        air_entry_water_potential = fill(air_entry_potential, length(depths)),
        saturated_hydraulic_conductivity = fill(sat_hydraulic_cond, length(depths)),
        campbell_b_parameter = fill(campbell_b, length(depths)),
        root_density,
    ),
)

problem = MicroVectorProblem(;
    model,
    points,
    dates,
    soil_profile,
    init = (; soil_moisture = fill(0.3, length(depths))),
)

@time output = solve(problem);

soil_T = output.soil_temperature[point=1]           # (Ti × depth)
air_T  = output.air_temperature[point=1, height=2]  # 1.2 m
rad    = output.global_radiation[point=1]
rh     = output.relative_humidity[point=1, height=2]
ws     = output.wind_speed[point=1, height=2]
soil_M = output.soil_moisture[point=1]

depth_labels = reshape(string.(depths), 1, :)

p1 = plot(collect(uconvert.(u"°C", soil_T)); label = depth_labels,
          title = "Soil temperature — Alice Springs (SILO)", legend = false)
p2 = plot(collect(uconvert.(u"°C", air_T));  title = "Air temperature (1.2 m)",   legend = false)
p3 = plot(collect(rad);    title = "Global radiation",           legend = false)
p4 = plot(collect(rh);     title = "Relative humidity (1.2 m)", legend = false)
p5 = plot(collect(ws);     title = "Wind speed (1.2 m) — CRUCL2 fallback", legend = false)
p6 = plot(collect(soil_M); label = depth_labels,
          title = "Soil moisture — Alice Springs (BARRA-C2)", legend = false)

display(plot(p1, p2, p3, p4, p5, p6; layout = (6, 1), size = (900, 1300), link = :x))
