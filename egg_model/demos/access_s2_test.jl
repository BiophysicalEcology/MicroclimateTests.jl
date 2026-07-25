# Stage 7 smoke test: solve a MicroVectorProblem driven by the local
# ACCESS-S2 reader (egg_model/src/access_s2.jl) at a single point, and sanity
# check the resulting soil forcing against SILO for an overlapping period
# (ACCESS-S2's two locally-available issue dates, 2024-06-01 and 2024-07-01,
# fall within SILO's real observed record, so this is a genuine forecast-vs-
# observed comparison, not just a syntax check).

using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_properties_model, example_soil_hydraulic_model
using Rasters, RasterDataSources, PointDataSources
using Dates, Unitful
using Plots

ENV["RASTERDATASOURCES_PATH"] = "z:/"#"c:/Spatial_Data"

include(joinpath(@__DIR__, "..", "src", "access_s2.jl"))

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:reference_temperature, :scalar),
    LayerSpec(:global_radiation, :scalar),
)

bendigo = (144.2826718, -36.7590183)
issue_date = Date(2024, 6, 1) #Date(2026, 7, 1)
member = 1
horizon_days = 214
dates = issue_date:Day(1):(issue_date + Day(horizon_days - 1))

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
        config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
    ),
    dem_source              = CRUCL2,
    weather_source          = AccessS2(issue_date, member),
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
    output_layers,
)

const SANDY_LOAM = (air_entry=1.5u"J/kg", b=3.1, Ksat=7.2e-4u"kg*s/m^3", field_capacity=0.21, wilting_point=0.10)
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
soil_profile = soil_profile_from_texture(SANDY_LOAM, depths)

println("Solving ACCESS-S2-driven microclimate at Bendigo (issue $issue_date, member $member)...")
problem = MicroVectorProblem(;
    model, points=[bendigo], dates, soil_profile,
    init = (; soil_moisture = fill(0.2, length(depths))),
)
@time access_result = solve(problem)

println("reference_temperature range: ", uconvert.(u"°C", extrema(access_result.reference_temperature[point=1])))
println("global_radiation range: ", extrema(access_result.global_radiation[point=1]))
println("soil_temperature[node=1] range: ", uconvert.(u"°C", extrema(access_result.soil_temperature[point=1, depth=1])))

# ── SILO comparison over the same real calendar dates, for sanity ──
silo_model = MicroMapModel(;
    micro_model = model.micro_model, dem_source = CRUCL2, weather_source = SILO,
    surface_albedo_source = 0.15, roughness_height_source = 0.004u"m",
    compute_terrain = false, output_layers,
)
silo_problem = MicroVectorProblem(;
    model = silo_model, points=[bendigo], dates, soil_profile,
    init = (; soil_moisture = fill(0.2, length(depths))),
)
@time silo_result = solve(silo_problem)
println("SILO reference_temperature range: ", uconvert.(u"°C", extrema(silo_result.reference_temperature[point=1])))

p = plot(uconvert.(u"°C", access_result.reference_temperature[point=1]); label="ACCESS-S2 (member $member)",
         title="Bendigo reference air temperature, issue $issue_date", xlabel="hour")
plot!(p, uconvert.(u"°C", silo_result.reference_temperature[point=1]); label="SILO (observed)")
savefig(p, joinpath(@__DIR__, "access_s2_vs_silo_bendigo.png"))
display(p)
