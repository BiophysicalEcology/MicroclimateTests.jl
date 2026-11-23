# Paths and model choices; facts about the datasets and sites are in settings/

const OZFLUX_DATA_DIR = joinpath(@__DIR__, "data")
const SRTM_ELEVATION_DIR = joinpath(OZFLUX_DATA_DIR, "srtm_elevation")
const BARRA_DIR = joinpath(OZFLUX_DATA_DIR, "barra")
const MODIS_LAI_DIR = joinpath(OZFLUX_DATA_DIR, "modis_lai")
# MODIS MCD15A3H leaf area index starts in July 2002
const FIRST_MODIS_YEAR = 2003

soil_properties_model_choice = Microclimate.example_soil_properties_model()
rainfall_entry_mode_choice = RateLimitedFrontRainfall()
# Soil solver iterations
convergence_choice = FixedIterationConvergence(1)
# Soil heat: adaptive implicit Euler, steps up to an hour (fixed hourly steps fail Newton convergence on some hours)
soil_energy_model_choice = SoilHeatTransport1D(; ode_solver=ImplicitEuler(),
    ode_kwargs=(; dt=60.0u"minute", dtmax=60.0u"minute", abstol=0.05u"K", reltol=0.0))
# Canopy solve tolerances; ground temperature "soil" (integrated with the soil) or "hour_start"
canopy_convergence_model_choice = AndersonCanopyConvergence(; leaf_balance_tolerance=1.5u"W/m^2",
    air_temperature_tolerance=0.1u"K", vapor_density_tolerance=0.07u"g/m^3",
    inverse_obukhov_length_tolerance=1e-2u"m^-1", on_failure=SmallestResidual(),
    ground_temperature_model=get(ENV, "OZFLUX_GROUND_TEMPERATURE", "soil") == "hour_start" ?
        HourStartGroundTemperature() : SoilIntegration())
# OzFlux Precip is per time step, not a daily total
rainfall_schedule_choice = HourlyRainfall()
# Soil moisture simulated, so it can be compared against Sws
soil_moisture_strategy_choice = DynamicSoilMoisture()
leaf_convection_model_choice = ElaborateLeafConvection()
interception_model_choice = LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")
# Ground stability: "louis", "richardson" or "obukhov"; free-convection floor in m/s
ground_turbulence_choice = GroundTurbulence(;
    stability=Dict("richardson" => GroundRichardson, "obukhov" => CanopyObukhovLength,
        "louis" => LouisGroundRichardson)[get(ENV, "OZFLUX_GROUND_STABILITY", "louis")](),
    free_convection_velocity=parse(Float64, get(ENV, "OZFLUX_FREE_CONVECTION", "0.04")) * u"m/s")
# In-canopy air profile: "exact" or "bulk" Raupach, or "ktheory"
canopy_air_profile_model_choice = Dict(
    "exact" => () -> RaupachLTheoryAirProfile(; formulation=ExactRaupachFormulation(), ground_turbulence_model=ground_turbulence_choice),
    "bulk" => () -> RaupachLTheoryAirProfile(; formulation=BulkRaupachFormulation(), ground_turbulence_model=ground_turbulence_choice),
    "ktheory" => () -> KTheoryAirProfile(; ground_turbulence_model=ground_turbulence_choice))[get(ENV, "OZFLUX_AIR_PROFILE", "exact")]()
# Mixing at canopy top as a multiple of the surface layer's
roughness_sublayer_choice = Microclimate.CanopyTopMixingEnhancement(;
    max_enhancement=parse(Float64, get(ENV, "OZFLUX_CANOPY_TOP_ENHANCEMENT", "6.0")))
# Ground roughness height in m for every site, or settings/sites.toml's when unset
ground_roughness_override = haskey(ENV, "OZFLUX_GROUND_ROUGHNESS") ? parse(Float64, ENV["OZFLUX_GROUND_ROUGHNESS"]) * u"m" : nothing
# Soil nodes in m added to the grid, comma separated
extra_soil_depths_m = parse.(Float64, filter(!isempty, split(get(ENV, "OZFLUX_EXTRA_SOIL_DEPTHS", ""), ",")))
# Longwave: "radiosity" (LayeredRadiosityExchange), "layered" (LayeredLongwaveExchange) or "black_body" (BlackBodyLongwaveExchange, as micropoint)
longwave_model_choice = Dict("radiosity" => () -> LayeredRadiosityExchange(), "layered" => () -> LayeredLongwaveExchange(),
    "black_body" => () -> Microclimate.BlackBodyLongwaveExchange())[get(ENV, "OZFLUX_LONGWAVE", "radiosity")]()
boundary_layer_model_choice = MoninObukhov(; stability_limit=Microclimate.MomentumCorrectionLimit(),
    stability_functions=ChengBrutsaert(; momentum_dyer_constant=16.0, heat_dyer_constant=16.0))
vapour_pressure_equation_choice = VapourPressureLookup()
canopy_wind_model_choice = ExponentialCanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio(; ratio=0.5))
# Soil texture to Campbell hydraulics, for SLGA and measured horizons
pedotransfer_model_choice = CampbellParticleSize()

save_outputs = true
outputs_dir = get(ENV, "OZFLUX_OUTPUTS_DIR", joinpath(@__DIR__, "outputs", "default"))
# Litter: "site" (settings/sites.toml's) or "none"
litter_choice = get(ENV, "OZFLUX_LITTER", "site")
# Litter moisture: "dry" (DryLitter) or "wet" (WetLitter)
litter_moisture_choice = get(ENV, "OZFLUX_LITTER_MOISTURE", "dry")
figure_format = "png"
emissivity = 0.97
