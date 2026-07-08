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

# Set true to fetch SILO's own primary layers via direct point queries
# (PointDataSources.jl's DataDrill API) instead of downloading/cropping the
# whole SILO grid. Only correct for a single point (as here). Wind still
# comes from CRUCL2's raster climatology (no point-query source for it, and
# none needed -- it's a coarse monthly grid, not something worth avoiding a
# small crop for). Requires `ENV["SILO_EMAIL"]` (or edit `silo_username`
# below) -- SILO's API requires an email address.
use_point_query = true
silo_username = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")

if use_point_query
    lon0, lat0 = points[1].lon, points[1].lat
    date_range = (Date(year(first(dates)), 1, 1), Date(year(last(dates)), 12, 31))

    function _point_raster(T, name; kw...)
        @info "  loading $T $name (point query)..."
        nt = getpoint(T, name; lon = lon0, lat = lat0, kw...)
        Raster(reshape(nt.values, 1, 1, length(nt.values)), (X([lon0]), Y([lat0]), Ti(nt.times)))
    end

    # CRUCL2's monthly wind climatology has no point-query source, so it's
    # still fetched via the ordinary raster crop, then expanded onto the
    # point's daily Ti axis by hand -- `Rasters.resample` (used by the
    # generic `_monthly_to_daily`) requires a proper X/Y grid and rejects a
    # degenerate 1x1 point "grid", so the nearest CRUCL2 cell is selected
    # directly instead of resampling onto it.
    function _point_monthly_to_daily(layer::AbstractRaster, ref::AbstractRaster, years, lon, lat)
        cell = layer[X(Near(lon)), Y(Near(lat))]
        ti = dims(ref, Ti)
        out = zeros(eltype(cell), 1, 1, length(ti))
        years_v = collect(years)
        d = 0
        for (yi, y) in enumerate(years_v), m in 1:12
            val = cell[(yi - 1) * 12 + m]
            for _ in 1:Dates.daysinmonth(y, m)
                d += 1
                out[1, 1, d] = val
            end
        end
        Raster(out, (dims(ref, X), dims(ref, Y), ti))
    end

    function MicroclimateMapper._load_weather(T::Type{<:SILO}, area::Extent, years)
        primary_fields = (:max_temp, :min_temp, :daily_rain, :rh_tmax, :rh_tmin, :radiation)
        met = NamedTuple{primary_fields}(map(
            name -> _point_raster(T, name; date = date_range, username = silo_username),
            primary_fields))
        ref = first(values(met))

        @info "  loading CRUCL2 wnd..."
        wnd_stack = MicroclimateMapper._load_layers(
            MicroclimateMapper.SingleFileBands(), CRUCL2, (:wnd,), area, years)
        wnd = _point_monthly_to_daily(wnd_stack.wnd, ref, years, lon0, lat0)

        stack = RasterStack(merge(met, (; wnd)))
        return Rasters.replace_missing(stack, NaN)
    end
end

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights,
        soil_properties_model = CampbelldeVriesSoilProperties(;
            mineral_conductivity,
            mineral_heat_capacity,
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
