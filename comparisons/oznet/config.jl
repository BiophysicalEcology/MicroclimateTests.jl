# config.jl — shared configuration for comparison.jl and debug_site.jl
#
# Site selection and the simulation date window stay in each entry-point
# script (that's what you change per run/debug session). Everything here is
# the stuff that should stay identical across both scripts.
#
# Unlike scan_snotel/, depths are NOT fixed here -- each OzNet site has its
# own SMDEP/TDEP node list (from oznetsiteinfo2.csv), so the depths array is
# built per-site in pipeline.jl's resolve_site(). Soil hydraulic properties
# are likewise not hand-specified here: they come from a live SLGA fetch +
# pedotransfer model per site (pedotransfer_model_choice below), not a fixed
# literature profile.

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/" #"Z:"

# Where the raw OzNet data lives (not part of this repo -- see README.md).
const OZNET_DATA_DIR = raw"C:\Users\mrke\Dropbox\Datasets\oznet"
const OZNET_SITEINFO = joinpath(OZNET_DATA_DIR, "oznetsiteinfo2.csv")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MODEL ALGORITHM CHOICES — change the model setup here, in one place     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

infiltration_algorithm_choice  = MatricPotentialAlgorithm()   # swap to MatricFluxPotentialAlgorithm() etc. to compare
convergence_choice             = FixedIterationConvergence(10)
rainfall_schedule_choice       = DailyRainfall()
soil_moisture_strategy_choice  = DynamicSoilMoisture()
dem_source_choice              = SRTM
compute_terrain_choice         = false

# Cosby et al. (1984) pedotransfer, same as the original NicheMapR-based
# study this pipeline reproduces (Kearney & Maino 2018) -- swap to
# CosbyUnivariate() or Campbell1985() to compare.
pedotransfer_model_choice = CosbyMultivariate()

# ── Weather source ────────────────────────────────────────────────────────────
# OzNet sites are all in the Murrumbidgee catchment (southeastern Australia),
# so every source below has coverage.
#   SILO                — daily, Australia-wide (~5 km), 1889-present; the
#                          direct analogue of the AWAP grids the original
#                          NicheMapR study used (micro_aust). Default.
#   BARRA{BARRAC2,AUST04} — hourly, ~4 km regional reanalysis, 1979-present.
#   NCEP{SurfaceFlux}   — daily, global reanalysis, coarser (T62 Gaussian
#                          grid); the direct analogue of micro_ncep.
weather_source_choice = SILO#BARRA{BARRAC2, AUST04}

const WEATHER_SOURCE_START = Dict(
    SILO                              => Date(1889, 1, 1),
    BARRA{BARRAC2, AUST04}            => Date(1979, 1, 1),
    RasterDataSources.NCEP{SurfaceFlux} => Date(1948, 1, 1),
)

# SILO's DataDrill point-query API requires an email address.
ENV["SILO_EMAIL"] = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")

# SILO/BARRA only: fetch via PointDataSources.jl OPeNDAP point queries
# (MicroclimateMapper.PointQuery, set in pipeline.jl) instead of downloading/
# cropping a whole raster per variable per site. No effect on NCEP.
use_opendap_points = true

# ── NicheMapR via R ───────────────────────────────────────────────────────────
# compare_nmr = false: skip NicheMapR entirely -- no R call, no CSV writing,
# no reading of nmr_outputs/ -- for a Julia-only run against observations.
# run_nmr = true: call run_nmr.R with Julia's forcing to produce NMR outputs.
# reuse_nmr = true: skip R call if metout.csv already present.
compare_nmr = true
run_nmr     = true
reuse_nmr   = true
nmr_out_dir = joinpath(@__DIR__, "nmr_outputs")

# ── Simulation cache ──────────────────────────────────────────────────────────
reuse_simulation = false
cache_simulation = true
sim_cache_dir    = joinpath(@__DIR__, "julia_cache")

# ── Weather forcing cache (on disk, separate from the in-memory weather_cache) ─
# Namespaced by weather source in the filename (see pipeline.jl) so switching
# weather_source_choice can't silently reuse another source's cached forcing.
reuse_weather     = true
cache_weather     = true
weather_cache_dir = joinpath(@__DIR__, "weather_cache")

# ── Save outputs to disk ──────────────────────────────────────────────────────
save_outputs  = true
outputs_dir   = joinpath(@__DIR__, "julia_outputs")
figure_format = "png"   # "png", "pdf", or "jpeg"

# Surface water pooling depth plot (Julia's surface_water vs NMR's CONDEP).
plot_surface_water = false

# REFL = 0.2 in the original NicheMapR study (Australia_Soil_Test.Rmd) --
# kept the same for comparability, rather than the 0.15 used in scan_snotel.
albedo            = 0.2
emissivity        = 0.95
roughness_height  = 0.004u"m"
profile_heights   = [0.01u"m", 2.0u"m"]   # [user height, reference height] -- matches NicheMapR's Refhyt default

# Small bounding box around each site's point for the live SLGA fetch --
# SLGA has no point-query API (see build_soil_profile's SLGA method), so a
# uniform profile is derived from this local area.
soil_area_buffer_deg = 0.05

# Cap simulation length -- matches the 4-year (2007-2010) window used in the
# original study. auto_date_range in comparison.jl/debug_site.jl anchors this
# to each site's own OzNet observation coverage instead of a fixed window.
max_sim_years = 4
