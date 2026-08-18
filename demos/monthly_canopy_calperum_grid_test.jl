using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile, DEFAULT_DEPTHS,
    example_soil_properties_model, example_soil_hydraulic_model
using Rasters, RasterDataSources
using Rasters.Extents: Extent
using Dates, Statistics, Unitful, Plots

# Same evenly-spaced canopy grid now in comparisons/ozflux/utils.jl's
# _build_heights (N_CANOPY_LAYERS points from 0 to canopy_height, then
# N_ABOVE_CANOPY_POINTS graded up to reference_height) -- testing whether
# this grid + a Calperum-like short canopy (canopy_height=3m, LAI=0.3,
# :bottom_heavy) hangs even in this demo's fast 12-representative-day solve.
const N_CANOPY_LAYERS = 15
const N_ABOVE_CANOPY_POINTS = 6
function _evenly_spaced_heights(canopy_height, reference_height, extra_heights_m=Float64[])
    canopy_m = ustrip(u"m", canopy_height)
    reference_m = ustrip(u"m", reference_height)
    canopy_pts = [(i / N_CANOPY_LAYERS) * canopy_m for i in 1:N_CANOPY_LAYERS]
    graded_above = canopy_m < reference_m ?
        collect(range(canopy_m, reference_m; length=N_ABOVE_CANOPY_POINTS + 1)[2:end]) : Float64[]
    above_pts = filter(>(canopy_m), unique(vcat(extra_heights_m, graded_above, reference_m)))
    return sort(vcat(canopy_pts, above_pts)) .* u"m"
end

# Vertical PAI density shapes (unnormalized, per-metre); plant_area_index_profile rescales to target_pai.
_pai_shape_top_heavy(layer_heights, canopy_height) = @. exp(3.0 * (layer_heights / canopy_height)) / u"m"
_pai_shape_bottom_heavy(layer_heights, canopy_height) = @. exp(-3.0 * (layer_heights / canopy_height)) / u"m"
_pai_shape_mid_crown(layer_heights, canopy_height) =
    @. exp(-0.5 * ((layer_heights - 0.7 * canopy_height) / (0.25 * canopy_height))^2) / u"m"
_pai_shape_uniform(layer_heights, canopy_height) = fill(1.0 / u"m", length(layer_heights))
const PAI_SHAPES = Dict(
    :top_heavy => _pai_shape_top_heavy, :bottom_heavy => _pai_shape_bottom_heavy,
    :mid_crown => _pai_shape_mid_crown, :uniform => _pai_shape_uniform,
)
function plant_area_index_profile(shape_kind, heights, canopy_height, target_pai)
    n_layers = count(h -> h <= canopy_height, heights)
    (; layer_heights) = Microclimate.canopy_layer_heights(heights, canopy_height, n_layers)
    density = PAI_SHAPES[shape_kind](layer_heights, canopy_height)
    raw = plant_area_index_from_density(density, heights, canopy_height)
    return raw .* (target_pai / sum(raw))
end

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

depths = DEFAULT_DEPTHS
# Calperum-like: canopy_height/reference_height/LAI/shape from
# comparisons/ozflux/config.jl (canopy_height from file metadata, ~3m;
# reference_height=20m; leaf_area_index=0.3; pai_shape=:bottom_heavy).
canopy_height = 3.0u"m"
reference_height = 20.0u"m"

heights = _evenly_spaced_heights(canopy_height, reference_height)
canopy_heights = heights[heights .<= canopy_height]
println("n canopy layers: ", length(canopy_heights), "  canopy_heights: ", canopy_heights)

points = [geocode("Calperum, Australia")]  # location doesn't matter for this grid/convergence test
dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)

leaf_area_index = 0.3
pai_shape = :bottom_heavy
plant_area_index = plant_area_index_profile(pai_shape, canopy_heights, canopy_height, leaf_area_index)
println("plant_area_index: ", plant_area_index)
plot(plant_area_index, canopy_heights, ylabel="Height (cm)", xlabel="Plant Area Index", title="Canopy Structure")

canopy_convergence_model_choice = PicardCanopyConvergence(;
    convergence=IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=20), relaxation=0.7)
leaf_convection_model_choice = ElaborateLeafConvection()
interception_model_choice = LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")
canopy_air_profile_model_choice = RaupachLTheoryAirProfile(; far_field_mode=Val(:exact))
canopy_soil_convergence_choice = IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=20)
canopy_model = MultilayerCanopy(; canopy_height, plant_area_index,
            shortwave_model=TwoStreamRadiation(),
            leaf_convection_model = leaf_convection_model_choice, interception_model = interception_model_choice,
            convergence_model = canopy_convergence_model_choice, air_profile_model = canopy_air_profile_model_choice)

output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
    LayerSpec(:air_temperature, :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:diffuse_fraction, :scalar),
    LayerSpec(:reference_temperature, :scalar),
    LayerSpec(:pressure, :scalar),
    LayerSpec(:zenith_angle, :solar),
    LayerSpec(:leaf_temperature, :canopy),
    LayerSpec(:canopy_air_temperature, :canopy, :air_temperature),
    LayerSpec(:canopy_wind_speed, :canopy, :wind_speed),
    LayerSpec(:canopy_relative_humidity, :canopy, :relative_humidity),
)


model = MicroMapModel(;
    micro_model = MicroModel(;
        hours                 = 0:1:23,
        depths,
        heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
        canopy_model,
        config                = MicroConfig(;
            convergence             = IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=20),
            soil_moisture_strategy  = PrescribedSoilMoisture(),
            canopy_soil_convergence = canopy_soil_convergence_choice,
        ),
    ),
    dem_source              = CRUCL2,
    weather_source          = CRUCL2,
    soil_moisture_source    = CPCSoil,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.01u"m",  # Calperum's roughness_height (config.jl SITE_LEGACY_PARAMS)
    compute_terrain         = false,
    output_layers,
)

# model = MicroMapModel(;
#     micro_model = MicroModel(;
#         hours                 = 0:1:23,
#         depths,
#         heights,
#         soil_properties_model = example_soil_properties_model(),
#         soil_hydraulic_model  = example_soil_hydraulic_model(),
#         #snow_model            = NoSnow(),
#         canopy_model,
#     ),
#     dem_source              = CRUCL2,
#     weather_source          = CRUCL2,
#     soil_moisture_source    = CPCSoil,
#     surface_albedo_source   = 0.15,
#     roughness_height_source = 0.01u"m",  # Calperum's roughness_height (config.jl SITE_LEGACY_PARAMS)
#     #compute_terrain         = false,
#     output_layers,
# )

problem = MicroVectorProblem(;
    model,
    points,
    dates,
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths))),
)

@time output = solve(problem)

site = 1
rad    = output.global_radiation[point=site][:]
soil_T = output.soil_temperature[point=site][:, :]           # (Ti × depth)
air_T  = output.air_temperature[point=site][:, 1:length(heights)]  # (Ti × height)
rh     = output.relative_humidity[point=site][:, 1:length(heights)]
ws     = output.wind_speed[point=site][:, 1:length(heights)]
canopy_air_T  = output.canopy_air_temperature[point=site][:, :]
leaf_T  = output.leaf_temperature[point=site][:, :]        
canopy_rh = output.canopy_relative_humidity[point=site][:, :]
canopy_ws = output.canopy_wind_speed[point=site][:, :]  # 1.2 m

rad    = output.global_radiation[point=site][1:24]
soil_T = output.soil_temperature[point=site][1:24, :]           # (Ti × depth)
air_T  = output.air_temperature[point=site][1:24, 1:length(heights)]  # (Ti × height)
rh     = output.relative_humidity[point=site][1:24, 1:length(heights)]
ws     = output.wind_speed[point=site][1:24, 1:length(heights)]
canopy_air_T  = output.canopy_air_temperature[point=site][1:24, :]
leaf_T  = output.leaf_temperature[point=site][1:24, :]        
canopy_rh = output.canopy_relative_humidity[point=site][1:24, :]
canopy_ws = output.canopy_wind_speed[point=site][1:24, :]  # 1.2 m

ti           = lookup(soil_T, Ti)
depth_labels = reshape(string.(depths), 1, :)

p1 = plot(collect(uconvert.(u"°C", soil_T)); label = depth_labels,
          title = "Soil temperature", legend = false)
p2 = plot(collect(uconvert.(u"°C", air_T));  title = "Air temperature (1.2 m)",   legend = false)
p3 = plot(collect(rad);    title = "Global radiation",           legend = false)
p4 = plot(collect(rh);     title = "Relative humidity (1.2 m)", legend = false)
p5 = plot(collect(ws);     title = "Wind speed (1.2 m)",        legend = false)
p6 = plot(collect(uconvert.(u"°C", canopy_air_T)); title = "Canopy air temperature", legend = false)
p7 = plot(collect(uconvert.(u"°C", leaf_T)); title = "Leaf temperature", legend = false)
p8 = plot(collect(canopy_rh); title = "Canopy relative humidity", legend = false) 
p9 = plot(collect(canopy_ws); title = "Canopy wind speed", legend = false)

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, p9; layout = (3, 3), size = (1450, 900), font_size = 12, link = :x))

p1 = plot(canopy_ws[6, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 04:00", legend=false)
p1_1 = plot!(ws[6, :], heights)
p2 = plot(canopy_ws[7, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 05:00", legend=false)
p2_1 = plot!(ws[7, :], heights)
p3 = plot(canopy_ws[8, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 06:00", legend=false)
p3_1 = plot!(ws[8, :], heights)
p4 = plot(canopy_ws[9, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 07:00", legend=false)
p4_1 = plot!(ws[9, :], heights)
p5 = plot(canopy_ws[10, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 08:00", legend=false)
p5_1 = plot!(ws[10, :], heights)
p6 = plot(canopy_ws[11, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 09:00", legend=false)
p6_1 = plot!(ws[11, :], heights)
p7 = plot(canopy_ws[12, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 10:00", legend=false)
p7_1 = plot!(ws[12, :], heights)
p8 = plot(canopy_ws[13, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 11:00", legend=false)
p8_1 = plot!(ws[13, :], heights)
p9 = plot(canopy_ws[14, :], reverse(canopy_heights), ylabel="Height (m)", xlabel="Wind speed (m/s)", title="Canopy wind speed profile, 12:00", legend=false)
p9_1 = plot!(ws[14, :], heights)

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, p9; layout = (3, 3), size = (1450, 900), font_size = 12, link = :y))

p1 = plot(collect(uconvert.(u"°C", canopy_air_T[6, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 04:00", legend=false)
p1_1 = plot!(collect(uconvert.(u"°C", air_T[6, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 04:00", legend=false)
p2 = plot(collect(uconvert.(u"°C", canopy_air_T[7, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 05:00", legend=false)
p2_1 = plot!(collect(uconvert.(u"°C", air_T[7, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 05:00", legend=false)
p3 = plot(collect(uconvert.(u"°C", canopy_air_T[8, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 06:00", legend=false)
p3_1 = plot!(collect(uconvert.(u"°C", air_T[8, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 06:00", legend=false)
p4 = plot(collect(uconvert.(u"°C", canopy_air_T[9, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 07:00", legend=false)
p4_1 = plot!(collect(uconvert.(u"°C", air_T[9, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 07:00", legend=false)
p5 = plot(collect(uconvert.(u"°C", canopy_air_T[10, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 08:00", legend=false)
p5_1 = plot!(collect(uconvert.(u"°C", air_T[10, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 08:00", legend=false)
p6 = plot(collect(uconvert.(u"°C", canopy_air_T[11, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 09:00", legend=false)
p6_1 = plot!(collect(uconvert.(u"°C", air_T[11, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 09:00", legend=false)
p7 = plot(collect(uconvert.(u"°C", canopy_air_T[12, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 10:00", legend=false)
p7_1 = plot!(collect(uconvert.(u"°C", air_T[12, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 10:00", legend=false)
p8 = plot(collect(uconvert.(u"°C", canopy_air_T[13, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 11:00", legend=false)
p8_1 = plot!(collect(uconvert.(u"°C", air_T[13, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 11:00", legend=false)
p9 = plot(collect(uconvert.(u"°C", canopy_air_T[14, :])), reverse(canopy_heights), ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 12:00", legend=false)
p9_1 = plot!(collect(uconvert.(u"°C", air_T[14, :])), heights, ylabel="Height (m)", xlabel="Air temperature (°C)", title="Canopy air temperature profile, 12:00", legend=false)

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, p9; layout = (3, 3), size = (1450, 900), font_size = 12, link = :y))
