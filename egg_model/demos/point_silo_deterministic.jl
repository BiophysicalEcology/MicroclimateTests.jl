# Single-point deterministic egg-development run driven by real SILO microclimate
# output, following demos/silo.jl's MicroVectorProblem setup and demos/ectotherm.jl's
# underground-forcing pattern (fixed nest depth, zero radiation, near-still air).

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

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "development.jl"))
include(joinpath(@__DIR__, "..", "src", "thermal.jl"))
include(joinpath(@__DIR__, "..", "src", "hydric.jl"))
include(joinpath(@__DIR__, "..", "src", "phases.jl"))
include(joinpath(@__DIR__, "..", "src", "forcing.jl"))

# load parameters
include(joinpath(@__DIR__, "..", "params", "chortoicetes.jl"))

# all cached/serialized run output goes here, not directly in demos/ -- one
# gitignore entry (egg_model/demos/output/) instead of per-file patterns.
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

# ── microclimate: SILO point run, extended with the soil layers the egg model needs ──

site = geocode("Mildura, Australia")
site_name = split(site.display_name, ",")[1]
points = [site]

dates = Date(2010, 1, 1):Day(1):Date(2011, 12, 31)
#oviposition_dates = [Date(2010, m, 1) for m in 1:1:12]
oviposition_dates = [Date(2010, 1, 15), Date(2010, 4, 15)]

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
          20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
diapause = false

if diapause
    nest_depth = 5.0u"cm"
    cold_hour_threshold = 30u"d"
    diapause_hour_threshold = 240.0u"d"
else
    nest_depth = 10.0u"cm"
    cold_hour_threshold = 0.0u"d"
    diapause_hour_threshold = 0.0u"d"
end

# get air temperature and rainfall for plot
#pointlayers(SILO)
rainfall = getpoint(SILO, :daily_rain; lon=site.lon, lat=site.lat, date = dates, username = "m.kearney@unimelb.edu.au")

# Just what egg_nest_forcing (forcing.jl) actually reads, plus soil_moisture
# (kept for a possible future moisture threshold, cheap to retain -- see
# forcing.jl's comment for why everything else here was dropped: unused, or
# inert under this egg model's SoilTemperatureEquals/SteadyDarcyFlux config).
output_layers = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:soil_water_potential, :soil),
    LayerSpec(:soil_thermal_conductivity, :soil),
    LayerSpec(:soil_humidity, :soil),
    LayerSpec(:reference_temperature, :scalar)
)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
        # DynamicSoilMoisture lets moisture actually respond to rain/evaporation
        # -- without it (the default), moisture stays pinned at its initial
        # value all year, which drove soil_water_potential to ~0 (saturated)
        # almost everywhere and made the Darcy conductivity power-law singular.
        config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
    ),
    dem_source              = SRTM,
    weather_source          = SILO,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
    output_layers,
)

use_cache = false   # cache the SILO point query to avoid repeated downloads

# Campbell & Norman (1998) Table 9.1 texture-class hydraulic parameters
# (matches NicheMapR's CampNormTbl9_1) -- air_entry is the table's magnitude
# (J/kg); negate when building a profile, since air_entry_water_potential
# itself must be negative (see hydric.jl's sign-convention note).
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

# builds a uniform-with-depth SoilProfile from one CAMPBELL_NORMAN_TEXTURES
# entry, as an alternative to a live SLGA/SoilGrids fetch when a specific
# texture is known/assumed for the site rather than spatially derived.
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

# soil_source: :slga (live SLGA fetch + Cosby pedotransfer, spatially derived,
# matches comparisons/oznet) or one of the CAMPBELL_NORMAN_TEXTURES keys
# (e.g. :loam) for a fixed literature texture instead.
soil_source = :sandy_loam
soil_profile = if soil_source === :slga
    b = 0.05   # soil_area_buffer_deg, matching comparisons/oznet
    soil_area = Extent(X=(site.lon - b, site.lon + b), Y=(site.lat - b, site.lat + b))
    println("Fetching SLGA soil texture...")
    build_soil_profile(SLGA, soil_area; depths, pedotransfer_model=CosbyMultivariate()).soil_profile
else
    soil_profile_from_texture(CAMPBELL_NORMAN_TEXTURES[soil_source], depths)
end
# Organic litter/crust cap on the top ~2.5 cm (nodes 1:3 of depths: 0, 1.25,
# 2.5 cm), matching NicheMapR's own default `cap=1` behaviour (see
# micro_global.R) -- overrides the mineral thermal properties there with
# low-conductivity, high-heat-capacity litter values, damping the diurnal
# amplitude near the surface.
# soil_profile.mineral_conductivity[1:3] .= 0.05u"W/m/K"
# soil_profile.mineral_heat_capacity[1:3] .= 1920.0u"J/kg/K"

# direct OPeNDAP point queries instead of downloading/cropping the whole annual
# SILO grid (~400MB/variable) -- matches demos/silo.jl's documented fast path.
ENV["SILO_EMAIL"] = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")
MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()

problem = MicroVectorProblem(;
    model, points, dates, soil_profile,
    init = (; soil_moisture = fill(0.2, length(depths))),
)

using Serialization

cache_path = joinpath(output_dir, "microclimate_cache.jls")
if isfile(cache_path) && use_cache
    println("Loading cached microclimate result from $cache_path...")
    result = deserialize(cache_path)
else
    println("Solving SILO microclimate...")
    @time output = solve(problem)
    result = (;
        soil_temperature          = collect(output.soil_temperature[point=1]),
        soil_moisture             = collect(output.soil_moisture[point=1]),
        soil_water_potential      = collect(output.soil_water_potential[point=1]),
        soil_thermal_conductivity = collect(output.soil_thermal_conductivity[point=1]),
        soil_humidity             = collect(output.soil_humidity[point=1]),
    )
    serialize(cache_path, result)
end

nest_node = nearest_node(nest_depth, depths)
# result arrays are hourly (length(dates)*24 rows), not daily -- day_range must
# span the full hourly series or the forcing silently extrapolates a constant
# past the first length(dates) *hours* (~15 days), not the full year.
day_range = 1:size(result.soil_temperature, 1)

environment_pars = example_environment_pars()
forcing = egg_nest_forcing(result, day_range, nest_node, environment_pars)

soil_hydraulics = (;
    air_entry_potential   = soil_profile.hydraulics.air_entry_water_potential[nest_node],
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity[nest_node],
    campbell_b            = soil_profile.hydraulics.campbell_b_parameter[nest_node],
)

# ── egg model, real fitted parameters from seasonal_hatch_forecast2.R's egg_pars ──

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
    StagedDesiccationLimit(;
        dry_mass, early_mass_factor,
        initial_egg_mass, late_mass_factor,
        ramp_start=quiescence_windows[1][1], ramp_end=quiescence_windows[1][2],
    ),
)
egg_model = EggModel(;
    development_model=dm, arrest_model=arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), survival_model, geometry,
)
initial_state = EggState(;
    egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
    maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
)
# ── run from a given oviposition date, reusing the same solved microclimate ──
# (the expensive part -- the SILO/MicroclimateMapper solve above -- happens
# once; trying different lay dates against it is cheap and doesn't re-solve it)

max_duration = 720.0u"d"
# forcing beyond the real data holds the last hour constant (see forcing.jl's
# ExtrapolationType.Constant) rather than ending -- cap tspan at whichever
# comes first, max_duration or the real forcing running out, so a late-year
# lay date doesn't spend most of its run in that fake constant tail.
forcing_end_hr = length(day_range) * 1.0u"hr"

function run_from_oviposition(oviposition_date; save_trajectory=false)
    start_hr = oviposition_offset(oviposition_date, dates)
    tspan = (start_hr, min(start_hr + max_duration, forcing_end_hr))
    initial_state = EggState(;
        egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
        maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(arrest),
    )
    simulate_egg(egg_model, pars, initial_state, soil_hydraulics, forcing, tspan; save_trajectory)
end

println("Trying oviposition dates: ", oviposition_dates)
oviposition_results = map(oviposition_dates) do oviposition_date
    r = run_from_oviposition(oviposition_date; save_trajectory=true)
    start_hr = oviposition_offset(oviposition_date, dates)
    outcome = if r.hatched
        # r.hatch_time is absolute (hours since first(dates)), so adding it
        # straight to first(dates) gives the actual calendar hatch date.
        hatch_date = first(dates) + Day(round(Int, ustrip(u"d", r.hatch_time)))
        "hatched after $(round(typeof(1.0u"d"), r.hatch_time - start_hr, digits=1)) (on $hatch_date)"
    elseif r.died
        "died of $(r.death_cause) after $(round(typeof(1.0u"d"), r.death_time - start_hr, digits=1))"
    else
        "did not hatch in $(max_duration)"
    end
    println("  lay date $oviposition_date -> $outcome")
    r
end

# sanity check: no quiescence windows at all (arrest reduces to diapause-only)
# -- confirms an empty quiescence_windows list is handled correctly.
no_quiescence_arrest = ProportionWindowArrest(;
    cold_temperature=arrest.cold_temperature, diapause_window=arrest.diapause_window,
    quiescence_windows=(),
    cold_hour_threshold=arrest.cold_hour_threshold, diapause_hour_threshold=arrest.diapause_hour_threshold,
)
no_quiescence_model = EggModel(;
    development_model=dm, arrest_model=no_quiescence_arrest, hydric_model=SteadyDarcyFlux(),
    hydric_stage_model=stage, thermal_model=SoilTemperatureEquals(), geometry,
)
let
    start_hr = oviposition_offset(first(oviposition_dates), dates)
    tspan = (start_hr, min(start_hr + max_duration, forcing_end_hr))
    initial_state = EggState(;
        egg_mass=pars.initial_egg_mass, egg_water_potential=-709.4682u"J/kg",
        maximum_mass_achieved=pars.initial_egg_mass, arrest_state=initial_arrest_state(no_quiescence_arrest),
    )
    r = simulate_egg(no_quiescence_model, pars, initial_state, soil_hydraulics, forcing, tspan)
    outcome = if r.hatched
        hatch_date = first(dates) + Day(round(Int, ustrip(u"d", r.hatch_time)))
        "hatched after $(round(typeof(1.0u"d"), r.hatch_time - start_hr, digits=1)) (on $hatch_date)"
    else
        "did not hatch in $(max_duration)"
    end
    println("no-quiescence sanity check -> $outcome (final development_fraction=$(r.final_state.development_fraction))")
end

# trajectories for every oviposition date, overlaid on the same panels.
using Plots
p1 = plot(; ylabel="development", title="Egg development by lay date at $site_name" * ", diapause = $diapause", legend=false)#, legend=:outerright)
p2 = plot(; ylabel="egg mass", legend=false)
#p3 = plot(; ylabel="water\npotential", legend=false)
p3 = plot(; ylabel="egg\ntemperature", xlabel="date", legend=false)
p4 = plot(; ylabel="rainfall, mm/day", xlabel="date", legend=false)
ref_temp_times = collect(lookup(output.reference_temperature, Ti))
rainfall_times = DateTime.(rainfall.times) .+ Hour(12)
plot!(p3, ref_temp_times, collect(uconvert.(u"°C", collect(output.reference_temperature[point=1]))); color=:lightgrey)
plot!(p4, rainfall_times, rainfall.values, color=:black)
trajectory_colors = palette(:default, length(oviposition_results))
for (i, (oviposition_date, r)) in
    enumerate(zip(oviposition_dates, oviposition_results))
    traj = r.trajectory
    # traj.t is simulation time in hours since dates[1] -- convert to actual
    # calendar DateTimes so each lay date's trajectory sits at when it really
    # occurred through the year, rather than all overlaid from a shared zero.
    actual_times = DateTime(first(dates)) .+ Dates.Second.(round.(Int, ustrip.(u"s", traj.t)))
    label = string(oviposition_date)
    trajectory_color = trajectory_colors[i]
    plot!(p1, actual_times, traj.development_fraction, ylims = (0, 1); label, color=trajectory_color)
    plot!(p2, actual_times, collect(uconvert.(u"mg", traj.egg_mass)), ylims = (0.0u"mg", pars.initial_egg_mass * 2.5); label, color=trajectory_color)
    #plot!(p3, actual_times, collect(uconvert.(u"J/kg", traj.egg_water_potential)); label)
    plot!(p3, actual_times, collect(uconvert.(u"°C", traj.temperature)); label, color=trajectory_color)
end
combined_plot = plot(p1, p2, p3, p4; layout=(4, 1), size=(1000, 900), link=:x, xlims=(DateTime(2010, 1, 1), DateTime(2010, 12, 31, 23, 59, 59)))
savefig(combined_plot, joinpath(@__DIR__, "point_silo_deterministic.png"))
display(combined_plot)

# # development-rate response curve: constant-temperature development time
# # (1/rate) across a temperature sweep, using the same development_model.
# temps_C = -5.0:1.0:45.0
# rates = [development_rate(dm, t * u"°C") for t in temps_C]
# dev_times_days = [uconvert(u"d", 1.0 / r) for r in rates]
# response_plot = plot(
#     collect(temps_C), ustrip.(dev_times_days);
#     label=false, xlabel="constant temperature (°C)", ylabel="development time (days)",
#     title="Development response curve", ylim=(0, 100),
# )
# savefig(response_plot, joinpath(@__DIR__, "development_response_curve.png"))
# display(response_plot)
