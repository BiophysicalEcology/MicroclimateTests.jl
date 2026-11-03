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
    inverse_obukhov_length_tolerance=1e-2u"m^-1", on_failure=SmallestResidual())
# One canopy solve per hour at the previous hour's ground temperature
canopy_soil_coupling_choice = FullColumnCoupling(; convergence=FixedIterationConvergence(1))
ground_coupling_choice = GroundTurbulence(; stability=GroundRichardson())
rainfall_schedule_choice     = HourlyRainfall()#DailyRainfall(; spread_hours=1)#HourlyRainfall()  # OzFlux Precip is per-timestep, not a daily total
soil_moisture_strategy_choice = DynamicSoilMoisture()  # simulate soil moisture, don't just prescribe it -- compared against Sws
leaf_convection_model_choice = ElaborateLeafConvection()  # or SimpleLeafConvection()
interception_model_choice = LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2") #NoInterception()  # or LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")
canopy_air_profile_model_choice = RaupachLTheoryAirProfile(; far_field_mode=Val(:exact)) # or Val(:bulk), or KTheoryAirProfile()
longwave_model_choice = LayeredRadiosityExchange() # LayeredLongwaveExchange(), LayeredRadiosityExchange(), AllPairsLongwaveExchange()
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
outputs_dir   = joinpath(@__DIR__, "julia_outputs")
figure_format = "png"
emissivity = 0.97
