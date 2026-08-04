# config.jl — shared configuration for comparison.jl and debug_site.jl
#
# Primary forcing is each site's own tower file (30-min, aggregated to
# hourly). A SILO gridded-forcing run is also available per site-year
# (pipeline.jl's run_site_silo) both as its own comparison against the tower
# obs (gridded vs point-observed skill) and, via run_site_gapfilled, as a
# gap-filler for missing tower forcing — see README for caveats (SILO is
# daily-native and has no wind or longwave layer).

const OZFLUX_DATA_DIR = joinpath(@__DIR__, "data")

ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/"

# ── SILO weather source (secondary comparison + forcing gap-fill) ───────────
weather_source_choice = SILO
ENV["SILO_EMAIL"] = get(ENV, "SILO_EMAIL", "m.kearney@unimelb.edu.au")
use_opendap_points = true  # PointDataSources.jl OPeNDAP point query, not a whole-raster download
dem_source_choice = SRTM
compute_terrain_choice = false

reuse_weather     = true
cache_weather     = true
weather_cache_dir = joinpath(@__DIR__, "weather_cache")

# ── Soil texture (live SLGA fetch + pedotransfer, like oznet) ───────────────
# SLGA has no point-query API -- a uniform profile is derived from a small
# bounding box around each site's coordinates (same approach as oznet).
pedotransfer_model_choice = CosbyMultivariate()
soil_area_buffer_deg = 0.05

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MODEL ALGORITHM CHOICES — change the model setup here, in one place     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

soil_hydraulic_model_choice  = Microclimate.example_soil_hydraulic_model()
soil_properties_model_choice = Microclimate.example_soil_properties_model()
convergence_choice           = FixedIterationConvergence(1)  # soil solver
canopy_convergence_choice    = IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=100)  # MultilayerCanopy's Picard loop
canopy_relaxation_choice     = 0.1
rainfall_schedule_choice     = HourlyRainfall()  # OzFlux Precip is per-timestep, not a daily total
soil_moisture_strategy_choice = DynamicSoilMoisture()  # simulate soil moisture, don't just prescribe it -- compared against Sws
leaf_convection_model_choice = ElaborateLeafConvection()  # or SimpleLeafConvection()
interception_model_choice = LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2") #NoInterception()  # or LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")

# ── "legacy" canopy_mode (see pipeline.jl prepare_site's canopy_mode kwarg) ──
# NoCanopy + a PAI-derived shade fraction + a wind-speed knockdown + a large
# horizon angle (sun only reaches the ground near-overhead) -- how forest
# sites were approximated before MultilayerCanopy existed. Run alongside
# canopy_mode=:full to quantify what the full canopy model actually buys.
# Per-site, not generic -- canopy openness varies enormously across these
# sites (Calperum's sparse mallee vs Cape Tribulation's closed rainforest),
# so a single wind_multiplier/horizon_angle would misrepresent one or the
# other. Falls back to DEFAULT_LEGACY_PARAMS for any site without an entry.
const DEFAULT_LEGACY_PARAMS = (extinction_coefficient=0.5, wind_multiplier=0.5, horizon_angle=80.0u"°")
const SITE_LEGACY_PARAMS = Dict(
    "CapeTribulation" => (extinction_coefficient=0.5, wind_multiplier=0.3, horizon_angle=85.0u"°"),  # closed tropical rainforest canopy
    "Calperum"        => (extinction_coefficient=0.5, wind_multiplier=0.7, horizon_angle=15.0u"°"),  # sparse, open mallee woodland
    "Whroo"           => (extinction_coefficient=0.5, wind_multiplier=0.5, horizon_angle=65.0u"°"),  # dry sclerophyll woodland, intermediate
    "Wallaby"         => (extinction_coefficient=0.5, wind_multiplier=0.4, horizon_angle=70.0u"°"),  # regrowth ash forest, denser than woodland
    "GWW"             => (extinction_coefficient=0.5, wind_multiplier=0.6, horizon_angle=25.0u"°"),  # open eucalypt woodland
    "Longreach"       => (extinction_coefficient=0.5, wind_multiplier=0.9, horizon_angle=5.0u"°"),   # grassland, effectively no canopy shelter
    "TiTreeEast"        => (extinction_coefficient=0.5, wind_multiplier=0.8, horizon_angle=10.0u"°"),  # sparse mulga/spinifex woodland
    "AliceSpringsMulga" => (extinction_coefficient=0.5, wind_multiplier=0.8, horizon_angle=10.0u"°"),
)
legacy_params(site_name) = get(SITE_LEGACY_PARAMS, site_name, DEFAULT_LEGACY_PARAMS)

# Base simulation depth nodes (m): NicheMapR's classic 10-node scheme (cm)
# plus an arithmetic midpoint between each pair -- same 19-node default as
# oznet's NMR_DEP19_CM. pipeline.jl's _build_depths adds each site's own
# observed Ts/Sws sensor depths on top of this (see _build_heights for the
# same idea applied to heights), so simulated nodes land exactly on
# comparison depths instead of only approximately.
const NMR_DEP_CM = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]
const SIM_DEPTHS_M = let d = NMR_DEP_CM
    vcat([d[1]], reduce(vcat, [[(d[i] + d[i+1]) / 2, d[i+1]] for i in 1:(length(d) - 1)])) ./ 100.0
end

# Canopy height nodes: FIXED absolute spacing (not scaled by canopy_height,
# unlike SIM_DEPTHS_M's fraction-of-depth-range approach) -- graded like a
# soil profile up to 2m (fine right at the ground, where the wind/
# temperature log-law gradient is steepest, coarsening quickly since the
# profile flattens fast), then uniform coarse steps from 2m to canopy_height.
# Fixed/absolute rather than fraction-of-canopy-height because the real
# sensor heights being compared against (1/2/4/8/16m) are themselves fixed
# absolute values -- a fraction-based grid gives worse near-ground
# resolution at exactly the tall-canopy sites (Whroo, Wallaby) with the best
# sensor data. First node is 0.05m, not 0.01m: site.roughness_height=0.02m,
# and atmospheric_surface_profile!'s bare-ground default requires every
# evaluation height to exceed that (see _build_heights).
const CANOPY_NEAR_GROUND_M = [0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50, 2.00]
const CANOPY_COARSE_STEP_M = 1.5

# Reference/user heights for MicroModel.heights -- must bracket every
# forcing/comparison height actually used. 20m covers Cape Tribulation's
# 45m tower and Calperum's 20m tower via the model's own reference-height
# extrapolation; per-site heights arrays are built in pipeline.jl instead
# when multi-height profile comparison needs specific intermediate heights.
const DEFAULT_HEIGHTS_M = [0.01, 2.0]

# ── Reference height for forcing (Ta/RH/Ws) ──────────────────────────────────
# Microclimate.jl always treats forcing as valid AT last(heights) (see
# pipeline.jl's _build_heights/resolve_site) -- i.e. whatever height sits at
# the top of the model's heights grid. The bare "Ta"/"Ws" tower variables are
# themselves QC-merged composites across several physical instrument heights
# (see each file's own `height` attribute + `description_L3`) which often do
# NOT match the site's `tower_height` global attribute:
#   CapeTribulation, GWW: bare Ta/Ws genuinely sit at tower_height already --
#     no entry needed, falls back to tower_height/bare "Ta"/"RH"/"Ws" below.
#   Calperum, Whroo: a genuine single-sensor tower-top instrument exists
#     (Ws_SONIC_Av/Ta_SONIC_Av at Calperum, Ws_CSAT/Ta_HMP_36m at Whroo) --
#     SITE_FORCING_VARS below prefers that over the lower/blended composite.
#   Wallaby: no tower-top wind/temp sensor exists at all -- Ws_CSAT (= bare
#     Ws) genuinely sits at only 5m, well below both canopy top and
#     tower_height=12m. reference_height must be 5m here, not tower_height,
#     or the same wrong-reference-height bug applies to the *only* data
#     available (there's nothing better to switch to).
#   Longreach: bare Ta/Ws already sit at ~2m; tower is only 3m (small
#     grassland mast) -- minor mismatch, set explicitly now the mechanism
#     exists.
const SITE_REFERENCE_HEIGHT_M = Dict(
    "Calperum"          => 20.0,
    "Whroo"             => 36.0,
    "Wallaby"           => 5.0,
    "Longreach"         => 2.0,
    "TiTreeEast"        => 9.8,   # bare Ta/Ws/AH all ~9.8m; tower_height=13.7m but no sensor sits there
    "AliceSpringsMulga" => 11.6,  # bare Ta/Ws/AH all ~11.6m; tower_height=13.7m but no sensor sits there
)

# Per-site override of which tower variable actually supplies Ta/Ws/AH
# forcing, for sites where the bare merged "Ta"/"Ws" composite's real height
# (see SITE_REFERENCE_HEIGHT_M comment) doesn't match reference_height.
# `ah`: humidity always re-derived from absolute humidity (via _rh_from_ah)
# rather than trusting the bare "RH" composite's height, which is equally
# unverified for these sites. Falls back to bare "Ta"/"Ws"/"RH" (unaltered
# height-mismatch and all) for any site not listed.
const SITE_FORCING_VARS = Dict(
    "Calperum" => (ta="Ta_SONIC_Av", ws="Ws_SONIC_Av", ah="AH_IRGA_Av"),  # true 20m tower-top sensors
    "Whroo"    => (ta="Ta_HMP_36m",  ws="Ws_CSAT",      ah="Ah"),          # true 36m; bare "Ah" is already 36m
)

# ── Per-site canopy leaf area index ──────────────────────────────────────────
# Guessed, no emprirical data.
const SITE_LEAF_AREA_INDEX = Dict{String,Float64}(
    "CapeTribulation" => 7.0,
    "Calperum"        => 0.35,
    "Whroo"           => 0.7,
    "Wallaby"         => 6.0,
    "GWW"             => 1.0,   # semi-arid eucalypt woodland, canopy_height=18m per file metadata
    "Longreach"       => 1.5,   # grassland, canopy_height=0.5m per file metadata
    "TiTreeEast"        => 0.5,  # sparse mulga/spinifex woodland, canopy_height=6.5m per file metadata
    "AliceSpringsMulga" => 0.5,  # same study/template as TiTreeEast -- sparse mulga woodland, canopy_height=6.5m
)

# ── Per-site vertical PAI density shape (pipeline.jl's PAI_SHAPES/
# plant_area_index_profile) -- replaces the flat scalar split-evenly-per-
# layer default with a shaped profile matching each site's structure.
# :top_heavy -- closed rainforest crown (Cape Tribulation): dense near the
#   top, sparse near the ground.
# :bottom_heavy -- open woodland (Whroo, Wallaby): sparse emergent crown,
#   denser sub-canopy/understorey (matches the GEDI-derived profiles).
# :mid_crown -- a distinct crown layer, open above and below.
# :uniform -- flat density (same as the old scalar-PAI behaviour).
# Falls back to :uniform for any site without an entry.
const SITE_PAI_SHAPE = Dict(
    "CapeTribulation" => :top_heavy,
    "Calperum"        => :bottom_heavy,
    "Whroo"           => :bottom_heavy,
    "Wallaby"         => :bottom_heavy,
    "GWW"             => :bottom_heavy,  # open eucalypt woodland, same reasoning as Whroo/Wallaby
    "Longreach"       => :uniform,       # grassland -- no crown structure
    "TiTreeEast"        => :bottom_heavy,  # sparse mulga woodland, simple shrub/small-tree structure
    "AliceSpringsMulga" => :bottom_heavy,
)

# ── Per-site canopy leaf spectral properties ──────────────────────────────────────────
# Guessed, no emprical data.
const SITE_LEAF_REFLECTANCE = Dict{String,Float64}(
    "CapeTribulation" => 0.25,
    "Calperum"        => 0.25,
    "Whroo"           => 0.15,
    "Wallaby"         => 0.25,
    "GWW"             => 0.25,
    "Longreach"       => 0.25,
    "TiTreeEast"        => 0.25,
    "AliceSpringsMulga" => 0.20,
)
const SITE_LEAF_TRANSMITTANCE = Dict{String,Float64}(
    "CapeTribulation" => 0.25,
    "Calperum"        => 0.10,
    "Whroo"           => 0.25,
    "Wallaby"         => 0.25,
    "GWW"             => 0.10,
    "Longreach"       => 0.25,
    "TiTreeEast"        => 0.10,
    "AliceSpringsMulga" => 0.10,
)

# Per-site surface albedo -- ground cover/soil colour varies a lot across
# these sites (Calperum's pale sandy mallee soil vs Cape Tribulation's dark
# wet rainforest litter). Falls back to DEFAULT_ALBEDO for any site without
# an entry.
const DEFAULT_ALBEDO = 0.20
const SITE_ALBEDO = Dict(
    "CapeTribulation" => 0.13,  # dark, wet closed-canopy rainforest litter
    "Calperum"        => 0.15,  # pale sandy mallee soil, sparse cover
    "Whroo"           => 0.10,  # box/ironbark woodland, leaf litter + some bare soil
    "Wallaby"         => 0.12,  # dense wet sclerophyll regrowth
    "GWW"             => 0.15,  # semi-arid eucalypt woodland, pale sandy soil
    "Longreach"       => 0.20,  # dry grassland, pale soil/cured grass
    "TiTreeEast"        => 0.20,  # central Australian red sand, sparse mulga/spinifex
    "AliceSpringsMulga" => 0.20,
)

# ── Per-site multi-height forcing/validation variables ──────────────────────
# No shared naming convention across sites for height-resolved sensors (see
# utils.jl's discover_depth_series for the one case that *is* uniform:
# depth-suffixed soil variables). Only listed for sites that actually have
# multiple heights; absent from this Dict just means no profile comparison.
const SITE_HEIGHT_SERIES = Dict{String,Vector{Tuple{Float64,String}}}(
    "Calperum" => [
        (2.0,  "Ta_HMP_2m"),
        (20.0, "Ta_SONIC_Av"),
        (2.0,  "AH_HMP_2m"),
        (20.0, "AH_IRGA_Av"),
        (2.0,  "Ws_RMY2m_Av"),
        (10.0, "Ws_RMY10m_Av"),
        (20.0, "Ws_SONIC_Av"),
    ],
    # canopy_height=28m, tower=35m -- 1/2/4/8/16m are sub-canopy, 32/36m above.
    # The best available real in-canopy profile test (Calperum only has one
    # below-canopy height).
    "Whroo" => [
        (1.0,  "Ta_HMP_1m"),  (1.0,  "Ws_RMY_1m"),
        (2.0,  "Ta_HMP_2m"),  (2.0,  "Ws_RMY_2m"),
        (4.0,  "Ta_HMP_4m"),  (4.0,  "Ws_RMY_4m"),
        (8.0,  "Ta_HMP_8m"),  (8.0,  "Ws_RMY_8m"),
        (16.0, "Ta_HMP_16m"), (16.0, "Ws_RMY_16m"),
        (32.0, "Ta_HMP_32m"), (32.0, "Ws_RMY_32m"),
        (36.0, "Ta_HMP_36m"), (36.0, "Ws_CSAT"),  # no Ws_RMY_36m in the file; Ws_CSAT is the true top-of-tower sonic
    ],
    # canopy_height 8-10m (attribute is a range; _parse_length_attr takes the
    # lower bound) -- single sub-canopy Ta sensor at 5m, no Ws at height here.
    "Wallaby" => [
        (5.0, "Ta_HMP_5m"),
    ],
    # canopy_height=6.5m -- genuine sub-canopy-to-near-top mast profile.
    "AliceSpringsMulga" => [
        (2.0,  "Ta_HMP_200cm"),  (2.0,  "Ws_SENTRY_200cm_Av"),
        (4.25, "Ta_HMP_425cm"),  (4.25, "Ws_SENTRY_425cm_Av"),
        (6.57, "Ta_HMP_657cm"),  (6.62, "Ws_SENTRY_662cm_Av"),
    ],
)

# ── QC ────────────────────────────────────────────────────────────────────────
# OzFlux QCFlag convention: 0 = good data (see each file's Flag00 attribute).
const OZFLUX_GOOD_QC = 0

# Physically implausible values pass QC anyway (e.g. CapeTribulation's RH
# sits flat near 0% for weeks at a time -- a tropical rainforest, never
# actually that dry -- interspersed with short spikes). Out-of-range values
# are replaced with `missing` at read time, same as the -9999 sentinel, so
# gap-filling (not the raw sensor fault) determines what feeds the model.
const PLAUSIBLE_RANGE = Dict(
    "RH"  => (1.0, 100.0),
    "Ta"  => (-10.0, 55.0),
    "Ws"  => (0.0, 60.0),
    "Fsd" => (0.0, 1500.0),
    "Fld" => (100.0, 600.0),
)

# ── 30-min -> hourly aggregation ─────────────────────────────────────────────
# All OzFlux files here are 30-min (time_step global attribute); the model
# is hourly-only (MultilayerCanopy.CANOPY_TIMESTEP, MicroModel.hours).
const OZFLUX_NATIVE_STEP_MINUTES = 30

# ── Save outputs to disk ──────────────────────────────────────────────────────
save_outputs  = true
outputs_dir   = joinpath(@__DIR__, "julia_outputs")
figure_format = "png"
site_albedo(site_name) = get(SITE_ALBEDO, site_name, DEFAULT_ALBEDO)
emissivity = 0.97
