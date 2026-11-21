# config.jl — paths and model algorithm choices for comparison.jl and single_site.jl. Facts
# about the datasets and sites are data, in tables/ (read by tables.jl).

const OZFLUX_DATA_DIR = joinpath(@__DIR__, "data")
const SRTM_ELEVATION_DIR = joinpath(OZFLUX_DATA_DIR, "srtm_elevation")
const BARRA_DIR = joinpath(OZFLUX_DATA_DIR, "barra")
# MODIS MCD15A3H leaf area index (sites.jl's modis_leaf_area_index) starts in July 2002, so runs
# start from FIRST_MODIS_YEAR
const MODIS_LAI_DIR = joinpath(OZFLUX_DATA_DIR, "modis_lai")
const FIRST_MODIS_YEAR = 2003

# SILO forcing choices (MicroclimateMapper) live in mapper.jl.

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MODEL ALGORITHM CHOICES — change the model setup here, in one place     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

soil_properties_model_choice = Microclimate.example_soil_properties_model()
rainfall_entry_mode_choice = RateLimitedFrontRainfall() #PoolCapacityRainfall(), ImplicitFluxRainfall(), RateLimitedFrontRainfall() example_soil_hydraulic_model CampbellSoilHydraulics
convergence_choice           = FixedIterationConvergence(1)  # soil solver
# Soil heat ODE: adaptive implicit Euler, steps up to an hour (fixed hourly steps fail Newton convergence on some hours)
soil_energy_model_choice = SoilHeatTransport1D(; ode_solver=ImplicitEuler(),
    ode_kwargs=(; dt=60.0u"minute", dtmax=60.0u"minute", abstol=0.05u"K", reltol=0.0))
# MultilayerCanopy's hourly leaf/air/vapor/stability solve, with per-kind residual tolerances
canopy_convergence_model_choice = AndersonCanopyConvergence(; leaf_balance_tolerance=1.5u"W/m^2",
    air_temperature_tolerance=0.1u"K", vapor_density_tolerance=0.07u"g/m^3",
    inverse_obukhov_length_tolerance=1e-2u"m^-1", on_failure=SmallestResidual(),
    ground_temperature_model=get(ENV, "OZFLUX_GROUND_TEMPERATURE", "soil") == "hour_start" ?
        HourStartGroundTemperature() : SoilIntegration())
# Stack API: no FullColumnCoupling / GroundTurbulence. The canopy solve converges the ground
# temperature itself (AndersonCanopyConvergence's default ground_temperature_model=SoilIntegration()),
# and ground-air exchange is RaupachLTheoryAirProfile's turbulent resistance + StantonGroundSublayer (default).
rainfall_schedule_choice     = HourlyRainfall()#DailyRainfall(; spread_hours=1)#HourlyRainfall()  # OzFlux Precip is per-timestep, not a daily total
soil_moisture_strategy_choice = DynamicSoilMoisture()  # simulate soil moisture, don't just prescribe it -- compared against Sws
leaf_convection_model_choice = ElaborateLeafConvection()  # or SimpleLeafConvection()
interception_model_choice = LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2") #NoInterception()  # or LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")
# Ground stability: "richardson", "obukhov" or "louis"; free-convection floor in m/s
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
# Ground roughness length in m for every site, or each site's own when unset
ground_roughness_override = haskey(ENV, "OZFLUX_GROUND_ROUGHNESS") ? parse(Float64, ENV["OZFLUX_GROUND_ROUGHNESS"]) * u"m" : nothing
# Soil nodes in m added to the grid, comma separated
extra_soil_depths_m = parse.(Float64, filter(!isempty, split(get(ENV, "OZFLUX_EXTRA_SOIL_DEPTHS", ""), ",")))
# "radiosity" (LayeredRadiosityExchange), "layered" (LayeredLongwaveExchange) or "black_body" (BlackBodyLongwaveExchange, as micropoint)
longwave_model_choice = Dict("radiosity" => () -> LayeredRadiosityExchange(), "layered" => () -> LayeredLongwaveExchange(),
    "black_body" => () -> Microclimate.BlackBodyLongwaveExchange())[get(ENV, "OZFLUX_LONGWAVE", "radiosity")]()
boundary_layer_model_choice = MoninObukhov(; stability_limit=Microclimate.MomentumCorrectionLimit(),
    stability_functions=ChengBrutsaert(; momentum_dyer_constant=16.0, heat_dyer_constant=16.0))
vapour_pressure_equation_choice = VapourPressureLookup()  # or GoffGratch()
# canopy_wind_model_choice = MixingLengthCanopyWindAttenuation(;
#                                         shelter_floor = 0.003,
#                                         shelter_pai_coefficient = 0.1,
#                                         mixing_length_coefficient = 2.0,
#                                         mixing_length_pai_coefficient = 0.25,
#                                         )
#canopy_wind_model_choice = ExponentialCanopyWindAttenuation(max_attenuation_coefficient = 2.879) #or MixingLengthCanopyWindAttenuation(; shelter_floor=..., ...)
canopy_wind_model_choice = ExponentialCanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio(; ratio=0.5))
# Soil texture to Campbell hydraulics, for SLGA and measured horizons (sites.toml.s soil_horizons)
pedotransfer_model_choice = CampbellParticleSize()  # or CosbyMultivariate(), CosbyUnivariate(); the only one scaling conductivity with bulk density

# Depth/height grid constants (NMR_DEP_CM, SIM_DEPTHS_M, N_CANOPY_LAYERS):
# see utils.jl.

# ── Save outputs to disk ──────────────────────────────────────────────────────
save_outputs  = true
outputs_dir   = get(ENV, "OZFLUX_OUTPUTS_DIR", joinpath(@__DIR__, "julia_outputs_stack"))
# "site" (sites.toml litter, needs a stack with LitterLayer) or "none" (NoLitter, or no litter field at all)
litter_choice = get(ENV, "OZFLUX_LITTER", isdefined(Microclimate, :LitterLayer) ? "site" : "none")
# "dry" (DryLitter) or "wet" (WetLitter, holding water and evaporating, needs a stack with WetLitter)
litter_moisture_choice = get(ENV, "OZFLUX_LITTER_MOISTURE", "dry")
figure_format = "png"
emissivity = 0.97
