# config.jl — shared configuration for comparison.jl and debug_site.jl
#
# Site selection (site_num/sites/auto_sites) and the simulation date window
# (sim_start/sim_end/auto_date_range) stay in each entry-point script, since
# those are exactly the things you change per run/debug session. Everything
# here is the stuff that should stay identical across both scripts.

ENV["RASTERDATASOURCES_PATH"] = "Z:"

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MODEL ALGORITHM CHOICES — change the model setup here, in one place     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

infiltration_algorithm_choice  = MatricPotentialAlgorithm()   # swap to MatricFluxPotentialAlgorithm() etc. to compare
convergence_choice             = FixedSoilTemperatureIterations(3)
rainfall_schedule_choice       = DailyRainfall()
soil_moisture_strategy_choice  = DynamicSoilMoisture()
dem_source_choice              = SRTM
compute_terrain_choice         = false

# ── Weather source ────────────────────────────────────────────────────────────
# All sites here are SCAN/SNOTEL stations (continental USA), so every source
# below has coverage — this is what lets the same site set be re-run against
# each one to see how much of the Julia-vs-NicheMapR-vs-obs spread is driven
# by forcing choice rather than by the microclimate model itself.
#   GRIDMET             — daily, CONUS only (~4 km), 1979-present
#   ERA5                — hourly, global reanalysis
#   NCEP{SurfaceFlux}   — daily, global reanalysis, coarser (T62 Gaussian grid)
# CHELSA is deliberately left out for now: it's a single 1981-2010 monthly
# climatology, not a per-day/year forcing series, so it doesn't fit this
# auto-date-range multi-year pipeline the way the three above do. Worth a
# separate demo, not a variant here.
weather_source_choice = GRIDMET

# Earliest date each weather source actually has data for — used to clamp
# auto_date_range below. Only add an entry when wiring up a new source.
const WEATHER_SOURCE_START = Dict(
    GRIDMET            => Date(1979, 1, 1),
    ERA5               => Date(1940, 1, 1),
    RasterDataSources.NCEP{SurfaceFlux}  => Date(1948, 1, 1),
)

# ── NicheMapR via R ───────────────────────────────────────────────────────────
# run_nmr = true: call run_nmr.R with Julia's forcing to produce NMR outputs.
# reuse_nmr = true: skip R call if metout.csv already present.
run_nmr     = true
reuse_nmr   = true
nmr_out_dir = joinpath(@__DIR__, "nmr_outputs")

# ── Simulation cache ──────────────────────────────────────────────────────────
# Set reuse_simulation = true to load a previously saved result instead of
# re-running the model.  Results are stored as
#   sim_cache_dir/<site>_<sim_start>_<sim_end>.jls
# Delete the .jls file (or set reuse_simulation = false) to force a fresh run.
reuse_simulation = false
cache_simulation = false
sim_cache_dir    = joinpath(@__DIR__, "julia_cache")

# ── Weather forcing cache (on disk, separate from the in-memory weather_cache) ─
# The forcing fetch (DEM + years-of-daily/hourly variables, per site) is the
# slowest part of the pipeline and, unlike the simulation cache above, doesn't
# depend on any Microclimate.jl model parameter (infiltration algorithm,
# convergence, etc. are only applied later, when building the final
# MicroProblem) — only on site location/date range and the
# dem_source/weather_source/compute_terrain choices above. So it's safe to
# persist once and reuse across parameter sweeps. Results are stored as
#   weather_cache_dir/<source>_<site>_<sim_start>_<sim_end>.jls
# (the source name is in the filename because switching weather_source_choice
# must not silently reuse another source's cached forcing.) Delete the .jls
# file (or set reuse_weather = false) to force a fresh fetch — e.g. after
# changing dem_source_choice/weather_source_choice/compute_terrain_choice.
reuse_weather   = true
cache_weather   = true
weather_cache_dir = joinpath(@__DIR__, "weather_cache")

# ── Initial snow depth ────────────────────────────────────────────────────────
# init_snow_from_obs = true: read initial snow depth from the first obs row.
# init_snow_from_obs = false (default): start from 0 cm, matching NicheMapR.
init_snow_from_obs = false

# ── Save outputs to disk ──────────────────────────────────────────────────────
save_outputs  = true
outputs_dir   = joinpath(@__DIR__, "julia_outputs")
figure_format = "png"   # "png", "pdf", or "jpeg"

# GRIDMET-only forcing parameters (grid elevation diagnostic + coverage check
# below). Left as a lazy Raster load so it's a no-op unless weather_source_choice
# == GRIDMET.
const GRIDMET_ELEV_RASTER = Raster("z:/GRIDMET/metdata_elevationdata.nc"; lazy = true)
albedo     = 0.15
emissivity = 0.95

# States outside gridMET coverage and gridMET bounding box (CONUS ~4 km grid).
# Only consulted in resolve_site() when weather_source_choice == GRIDMET —
# ERA5/NCEP are global, so no equivalent check is needed for them.
const GRIDMET_EXCLUDED_STATES = Set(["AK", "HI", "PR"])
const GRIDMET_LON_MIN = -124.9;  const GRIDMET_LON_MAX = -66.8
const GRIDMET_LAT_MIN =   24.1;  const GRIDMET_LAT_MAX =  49.4

# ── Batch-prefetch spatial tiling ─────────────────────────────────────────────
# prefetch_weather_batch! groups sites sharing a date range into one
# MicroVectorProblem so the weather source/SRTM are only read once per group
# instead of once per site. But the SRTM DEM load covers the bounding box of
# *every* point in the group at once — combining sites scattered across the
# whole country (e.g. a random site_subset) blows that box up to a multi-GB
# read. Capping each tile's extent keeps every batch's bounding box close to
# what a single site's own buffer already needed, so batching only combines
# sites that are genuinely close together.
const PREFETCH_TILE_DEG          = 3.0   # max lon/lat span of sites grouped into one batch
const PREFETCH_MAX_SITES_PER_TILE = 40   # secondary cap even within one tile

# ── Soil properties (19 fine nodes) ──────────────────────────────────────────
# Values from SNOTEL 329 (Utah) — adjust per site as needed.
# Depths: [0, 1.25, 2.5, 3.75, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 40, 50, 75, 100, 150, 200] cm
depths = ([0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0] ./ 100.0) .* u"m"

bulk_density          = 1.3u"Mg/m^3"
saturation_moisture   = 0.4922u"m^3/m^3"
mineral_density       = 2.56u"Mg/m^3"
mineral_conductivity  = [0.2, 0.2, 0.2, 1.35, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5,
                          2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]u"W/m/K"
mineral_heat_capacity = [1920.0, 1920.0, 1920.0, 1395.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0,
                          870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0]u"J/kg/K"
air_entry_potential   = 1.1u"J/kg"
sat_hydraulic_cond    = 0.0037u"kg*s/m^3"
campbell_b            = 4.5
root_density          = [0.0, 0.0, 82000.0, 80000.0, 78000.0, 74000.0, 71000.0, 64000.0, 58000.0, 48000.0,
                          40000.0, 18000.0, 9000.0, 6000.0, 8000.0, 4000.0, 4000.0, 0.0, 0.0]u"m/m^3"

# Plant/root parameters for the Campbell water-uptake model (Campbell 1985 defaults).
root_resistance             = 2.5e10u"m^3/kg/s"
stomatal_closure_potential  = -1500.0u"J/kg"
leaf_resistance              = 2.0e6u"m^4/kg/s"
stomatal_stability_parameter = 10.0
root_radius                  = 0.001u"m"

# ── Julia snow model parameters ───────────────────────────────────────────────
snow_temp_threshold = 1.5u"°C"
snow_density_init   = 0.375u"g/cm^3"
snow_melt_factor    = 1.0
undercatch          = 1.0
rain_multiplier     = 1.0
rain_melt_factor    = 0.0125
density_function    = (0.5979, 0.2178, 0.001, 0.0038)
snow_conductivity   = 0.0u"W/m/K"
canopy_interception = 0.0

# ── Obs de-spiking ────────────────────────────────────────────────────────────
# Each entry: column_symbol => (thresh=…, clip_negative=…, tail_quantile=…)
# thresh: max allowed step change between any two values within ±halfwin steps.
# tail_quantile: fraction to cut from the upper tail (0 = no cut).
obs_despike_specs = [
    :SNWD_cm   => (thresh=50.0,  clip_negative=true,  tail_quantile=0.01),
    :WTEQ_cm   => (thresh=20.0,  clip_negative=true,  tail_quantile=0.01),
    :STO_5cm   => (thresh=10.0,  clip_negative=false, tail_quantile=0.0),
    :STO_10cm  => (thresh=10.0,  clip_negative=false, tail_quantile=0.0),
    :STO_20cm  => (thresh=10.0,  clip_negative=false, tail_quantile=0.0),
    :STO_50cm  => (thresh=10.0,  clip_negative=false, tail_quantile=0.0),
    :STO_100cm => (thresh=10.0,  clip_negative=false, tail_quantile=0.0),
    :SMS_5cm   => (thresh=30.0,  clip_negative=true,  tail_quantile=0.0),
    :SMS_10cm  => (thresh=30.0,  clip_negative=true,  tail_quantile=0.0),
    :SMS_20cm  => (thresh=30.0,  clip_negative=true,  tail_quantile=0.0),
    :SMS_50cm  => (thresh=30.0,  clip_negative=true,  tail_quantile=0.0),
    :SMS_100cm => (thresh=30.0,  clip_negative=true,  tail_quantile=0.0),
]

# ── Terrain / profile heights ─────────────────────────────────────────────────
roughness_height = 0.004u"m"
profile_heights  = [0.01u"m", 2.0u"m"]   # [user height, reference height]

# ── Soil thermal model constants ──────────────────────────────────────────────
# de Vries (1963) shape factor: 0.1 for mineral soil, 0.33 for organic soils.
# NicheMapR hardcodes its own shape factor internally and doesn't expose it
# via microinput, so there's no NMR-side value to match here.
de_vries_shape_factor = 0.1
recirculation_power   = 4.0
return_flow_threshold = 0.162
