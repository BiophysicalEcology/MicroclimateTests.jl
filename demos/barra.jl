using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_hydraulic_model
using Rasters, RasterDataSources, PointDataSources
using NCDatasets
using Dates, Statistics, Unitful, Plots
using DataInterpolations: CubicSpline, ExtrapolationType

ENV["RASTERDATASOURCES_PATH"] = "Z:" # "c:/Spatial_Data/"

depths = ([0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0] ./ 100.0) .* u"m"
heights = [0.01, 1.2]u"m"

mineral_density       = 2.56u"Mg/m^3"
mineral_conductivity  = [0.2, 0.2, 0.2, 1.35, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5,
                          2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]u"W/m/K"
mineral_heat_capacity = [1920.0, 1920.0, 1920.0, 1395.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0,
                          870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0]u"J/kg/K"

points = [geocode("Daly Waters, NT, Australia")]

dates = Date(2025, 1, 1):Day(1):Date(2025, 12, 31)

# Set true to fetch weather via direct point queries (PointDataSources.jl)
use_point_query = false

if use_point_query
    lon0, lat0 = points[1].lon, points[1].lat
    date_range = (Date(year(first(dates)), 1, 1), Date(year(last(dates)), 12, 31))

    # One 1x1xTi point-Raster per field, matching the shape every other
    # WeatherLoader returns.
    function _point_raster(T, name; kw...)
        @info "  loading $T $name (point query)..."
        nt = getpoint(T, name; lon=lon0, lat=lat0, kw...)
        Raster(reshape(nt.values, 1, 1, length(nt.values)), (X([lon0]), Y([lat0]), Ti(nt.times)))
    end

    # mrsol/tsl are 4-D (lon, lat, depth, time) -- PointDataSources' generic
    # `getpoint` assumes a plain 3-D (lon, lat, time) variable (indexes
    # `var[loni, lati, i0:i1]`), so this reimplements the same OPeNDAP
    # point-lookup with the extra depth axis carried through. Returns
    # `values` as a (depth, time) matrix, one row per BARRA soil layer.
    function _point_barra_layered(T, layer::Symbol; date)
        @info "  loading $T $layer (point query, layered)..."
        start, finish = PointDataSources._point_daterange(date)
        dates_seq = PointDataSources._periods_spanned(RasterDataSources.date_step(T), Date(start), Date(finish))
        times = DateTime[]
        cols = Vector{Vector{Float64}}()
        for d in dates_seq
            PointDataSources._retry() do
                NCDataset(PointDataSources._barra_opendap_url(T, layer; date = d)) do ds
                    lonv, latv = Float64.(ds["lon"][:]), Float64.(ds["lat"][:])
                    loni = PointDataSources._nearest_index(lonv, lon0, PointDataSources._half_cell(lonv))
                    lati = PointDataSources._nearest_index(latv, lat0, PointDataSources._half_cell(latv))
                    var = ds[string(layer)]
                    t = DateTime.(ds["time"][:])
                    i0, i1 = findfirst(>=(start), t), findlast(<=(finish), t)
                    (i0 === nothing || i1 === nothing || i0 > i1) && return nothing
                    vals = Float64.(var[loni, lati, :, i0:i1])  # (depth, ntime_slice)
                    append!(times, t[i0:i1])
                    for k in axes(vals, 2)
                        push!(cols, vals[:, k])
                    end
                end
            end
        end
        return (; times, values = reduce(hcat, cols))
    end

    # Point-specific equivalents of MicroclimateMapper's `_daily_to_hourly_midnight`/
    # `_static_to_ti` -- those call `Rasters.resample` (GDAL), which requires
    # regular grid dimensions and rejects a degenerate 1x1 point "grid".
    # Resampling a point onto itself is a no-op anyway, so skip GDAL entirely.
    function _point_daily_to_hourly_midnight(daily::Raster, ref::Raster)
        ti = dims(ref, Ti)
        out = zeros(eltype(daily), 1, 1, length(ti))
        for d in 1:size(daily, 3)
            out[1, 1, (d - 1) * 24 + 1] = daily[1, 1, d]
        end
        Raster(out, (dims(ref, X), dims(ref, Y), ti))
    end

    function _point_static_to_ti(value::Real, ref::Raster)
        ti = dims(ref, Ti)
        Raster(fill(value, 1, 1, length(ti)), (dims(ref, X), dims(ref, Y), ti))
    end

    function MicroclimateMapper._load_weather(T::Type{<:BARRA{P, D}}, area::Rasters.Extents.Extent, years) where {P, D}
        hourly_fields = (:tas, :sfcWind, :hurs, :psl, :rsds, :rlds)
        met = NamedTuple{hourly_fields}(map(name -> _point_raster(T, name; date = date_range), hourly_fields))
        ref = first(values(met))

        pr_daily = _point_raster(BARRA{P, D, Day}, :pr; date = date_range)
        pr = _point_daily_to_hourly_midnight(pr_daily, ref)

        @info "  loading $T orog (point query)..."
        orog_pt = getpoint(T, :orog; lon = lon0, lat = lat0)
        orog = _point_static_to_ti(orog_pt.value, ref)

        stack = RasterStack(merge(met, (; pr, orog)))
        return Rasters.replace_missing(stack, NaN)
    end
end

# =============================================================================
# Fetch BARRA's own soil layers -- used both to seed initial conditions
# (below) and later to compare against the model's own predictions.
# =============================================================================
# BARRA (JULES land surface model) reports mrsol/tsl on 4 fixed layers —
# not available at Hour frequency, fetched here at Day frequency instead.
# Layer bounds: 0-10, 10-35, 35-100, 100-300 cm (BOM/JULES convention).
barra_layer_bounds_cm = ((0.0, 10.0), (10.0, 35.0), (35.0, 100.0), (100.0, 300.0))
barra_layer_mid_cm    = [(lo + hi) / 2 for (lo, hi) in barra_layer_bounds_cm]
barra_layer_thick_m   = [(hi - lo) / 100 for (lo, hi) in barra_layer_bounds_cm]

lon, lat = points[1].lon, points[1].lat
barra_area = Extent(X = (lon - 0.1, lon + 0.1), Y = (lat - 0.1, lat + 0.1))
year_months = Date(year(first(dates)), month(first(dates))):Month(1):Date(year(last(dates)), month(last(dates)))

if use_point_query
    mrsol_vals = _point_barra_layered(BARRA{BARRAC2, AUST04, Day}, :mrsol; date = date_range).values
    tsl_vals   = _point_barra_layered(BARRA{BARRAC2, AUST04, Day}, :tsl;   date = date_range).values
    mrsol_point = [mrsol_vals[i, :] for i in 1:4]
    tsl_point   = [tsl_vals[i, :]   for i in 1:4]
    mrsos_point = vec(collect(_point_raster(BARRA{BARRAC2, AUST04}, :mrsos; date = date_range)))
else
    function _fetch_barra_daily_point(layer, year_months, area, lon, lat)
        slices = map(year_months) do d
            path = getraster(BARRA{BARRAC2, AUST04, Day}, layer; date = d)
            read(crop(Raster(path; name = layer, lazy = true); to = area, touches = true))
        end
        full = cat(slices...; dims = Ti)
        layer_dim = first(Rasters.otherdims(full, (X, Y, Ti)))
        DimType = Rasters.DimensionalData.basetypeof(layer_dim)
        [full[X(Near(lon)), Y(Near(lat)), DimType(i)] for i in 1:length(layer_dim)]
    end

    # Plain 3-D (lon, lat, time) hourly variable -- no depth axis to peel off.
    function _fetch_barra_point(T, layer, year_months, area, lon, lat)
        slices = map(year_months) do d
            path = getraster(T, layer; date = d)
            read(crop(Raster(path; name = layer, lazy = true); to = area, touches = true))
        end
        full = cat(slices...; dims = Ti)
        vec(collect(full[X(Near(lon)), Y(Near(lat))]))
    end

    mrsol_point = _fetch_barra_daily_point(:mrsol, year_months, barra_area, lon, lat)
    tsl_point   = _fetch_barra_daily_point(:tsl,   year_months, barra_area, lon, lat)
    mrsos_point = _fetch_barra_point(BARRA{BARRAC2, AUST04}, :mrsos, year_months, barra_area, lon, lat)
end

# kg/m^2 over the layer -> volumetric water content (m^3/m^3).
barra_moisture = [mrsol_point[i] ./ (barra_layer_thick_m[i] * 1000.0) for i in 1:4]
barra_temp_C   = [tsl_point[i] .- 273.15 for i in 1:4]
# mrsos is the same 0-10 cm surface layer as mrsol's layer 1, but hourly.
mrsos_moisture = mrsos_point ./ (barra_layer_thick_m[1] * 1000.0)

depths_cm = ustrip.(u"cm", depths)

# =============================================================================
# Soil hydraulic profile: from SLGA soil texture, or specified directly
# =============================================================================
# Set true to derive the profile from SLGA soil texture via build_soil_profile;
# set false to specify soil hydraulic properties directly (uniform values),
# as before build_soil_profile existed.
use_slga_soil_profile = true

if use_slga_soil_profile
    (; soil_profile, campbell_b, air_entry_potential, saturated_conductivity, field_capacity, wilting_point) =
        build_soil_profile(SLGA, barra_area; depths, mineral_density)
    saturation_fraction = 1 .- ustrip.(soil_profile.bulk_density ./ soil_profile.mineral_density)
    @info "SLGA hydraulic properties" campbell_b air_entry_potential saturated_conductivity field_capacity wilting_point
else
    bulk_density              = 1.3u"Mg/m^3"
    air_entry_potential       = 1.1u"J/kg"
    saturated_conductivity    = 0.0037u"kg*s/m^3"
    campbell_b                = 4.5
    root_density              = Microclimate.example_campbell_hydraulic_profile().root_density

    saturation_fraction = 1 - ustrip(bulk_density / mineral_density)
    soil_profile = SoilProfile(;
        bulk_density    = fill(bulk_density, length(depths)),
        mineral_density = fill(mineral_density, length(depths)),
        mineral_conductivity,
        mineral_heat_capacity,
        hydraulics = CampbellHydraulicProfile(;
            air_entry_water_potential         = fill(air_entry_potential, length(depths)),
            saturated_hydraulic_conductivity  = fill(saturated_conductivity, length(depths)),
            campbell_b_parameter              = fill(campbell_b, length(depths)),
            root_density,
        ),
    )
end

# =============================================================================
# Initial soil moisture/temperature from BARRA's own first-day values
# =============================================================================
# BARRA only gives 4 fixed-layer values; the model has 19 depth nodes, so the
# 4 layer-midpoint values are splined across all of them. The BARRA layer
# midpoints (5, 22.5, 67.5, 200 cm) are irregularly spaced, and
# `extrapolation = ExtrapolationType.Extension` continues evaluating the
# boundary segment's cubic beyond the outermost knots (rather than erroring)
# -- the shallowest model depths (0-3.75 cm) sit above BARRA's shallowest
# layer midpoint (5 cm). Same pattern as
# MicroclimateMapper's `build_soil_profile` (`_interpolate_onto_depths` in
# soil_profile_builder.jl) uses for SLGA texture -- DataInterpolations.jl is
# already a MicroclimateMapper dependency, so no new package is added here.
moisture_day1_layers = [barra_moisture[i][1] for i in 1:4]
temp_day1_layers_K   = [tsl_point[i][1] for i in 1:4]  # tsl is native Kelvin

moisture_spline = CubicSpline(moisture_day1_layers, barra_layer_mid_cm; extrapolation = ExtrapolationType.Extension)
temp_spline     = CubicSpline(temp_day1_layers_K,   barra_layer_mid_cm; extrapolation = ExtrapolationType.Extension)

init_soil_moisture = clamp.(moisture_spline.(depths_cm), 0.01, saturation_fraction)
init_soil_temperature = temp_spline.(depths_cm) .* u"K"

@info "init: soil moisture (BARRA day 1, splined): $(round.(init_soil_moisture; digits=3))"
@info "init: soil temperature (BARRA day 1, splined): $(round.(ustrip.(u"°C", init_soil_temperature); digits=1)) °C"

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
    dem_source              = BARRA{BARRAC2, AUST04},
    weather_source          = BARRA{BARRAC2, AUST04},
    # weather_source          = SILO,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
)

problem = MicroVectorProblem(;
    model,
    points,
    dates,
    soil_profile,
    init = (; soil_moisture = init_soil_moisture, soil_temperature = init_soil_temperature),
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
          title = "Soil temperature — Alice Springs (BARRA-C2)", legend = false)
p2 = plot(collect(uconvert.(u"°C", air_T));  title = "Air temperature (1.2 m)",   legend = false)
p3 = plot(collect(rad);    title = "Global radiation",           legend = false)
p4 = plot(collect(rh);     title = "Relative humidity (1.2 m)", legend = false)
p5 = plot(collect(ws);     title = "Wind speed (1.2 m)",        legend = false)
p6 = plot(collect(soil_M); label = depth_labels,
          title = "Soil moisture — Alice Springs (BARRA-C2)", legend = false)

display(plot(p1, p2, p3, p4, p5, p6; layout = (6, 1), size = (900, 1300), link = :x))

# =============================================================================
# Compare modelled soil moisture/temperature against BARRA's own soil layers
# =============================================================================
# Nearest model depth node to each BARRA layer midpoint.
nearest_idx = [argmin(abs.(depths_cm .- mid)) for mid in barra_layer_mid_cm]

model_moisture = output.soil_moisture[point=1]            # (Ti × depth)
model_temp_C   = ustrip.(u"°C", output.soil_temperature[point=1])

# Explicit date axes per series -- model/mrsos are hourly, mrsol/tsl are
# daily, and plotting plain (`collect`-ed) arrays together previously left
# each series on its own implicit 1:N index, badly misaligning the daily
# series against the hourly ones. Building real DateTime x-axes fixes the
# alignment and reads as actual calendar dates instead of step indices.
hourly_dates = [DateTime(first(dates)) + Hour(h) for h in 0:(size(model_moisture, 1) - 1)]
daily_dates  = [DateTime(first(dates)) + Day(d) for d in 0:(length(barra_moisture[1]) - 1)]
mrsos_dates  = [DateTime(first(dates)) + Hour(h) for h in 0:(length(mrsos_moisture) - 1)]

p7 = plot(layout = (4, 1), size = (900, 900), link = :x)
p8 = plot(layout = (4, 1), size = (900, 900), link = :x)
for (i, idx) in enumerate(nearest_idx)
    lo, hi = barra_layer_bounds_cm[i]
    label = "$(lo)-$(hi) cm (model @ $(depth_labels[idx]))"
    plot!(p7[i], hourly_dates, collect(model_moisture[:, idx]);   label = "model", title = label)
    plot!(p7[i], daily_dates,  collect(barra_moisture[i]);        label = "BARRA mrsol", linestyle = :dash)
    i == 1 && plot!(p7[i], mrsos_dates, collect(mrsos_moisture);  label = "BARRA mrsos (hourly)", linestyle = :dot)
    plot!(p8[i], hourly_dates, collect(model_temp_C[:, idx]);     label = "model", title = label)
    plot!(p8[i], daily_dates,  collect(barra_temp_C[i]);          label = "BARRA tsl", linestyle = :dash)
end
display(plot(p7; plot_title = "Soil moisture (m³/m³) — model vs BARRA"))
display(plot(p8; plot_title = "Soil temperature (°C) — model vs BARRA"))
