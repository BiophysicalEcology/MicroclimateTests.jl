# config.jl — shared configuration for comparison.jl and single_site.jl
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

soil_properties_model_choice = Microclimate.example_soil_properties_model()
convergence_choice           = FixedIterationConvergence(1)  # soil solver
# MultilayerCanopy's hourly leaf/air-temperature solve.
canopy_convergence_model_choice = PicardCanopyConvergence(;
    convergence=IterationToleranceConvergence(; tolerance=0.05u"K", max_iterations_per_day=80), relaxation=0.5) # set relaxation high for short canopies
canopy_soil_convergence_choice = IterationToleranceConvergence(; tolerance=0.05u"K", max_iterations_per_day=80) #FixedIterationConvergence(1)
canopy_soil_relaxation_choice = 1.0  # under-relaxation on canopy_soil_convergence_choice's ground_temperature update, 1.0 = none
rainfall_schedule_choice     = HourlyRainfall()  # OzFlux Precip is per-timestep, not a daily total
soil_moisture_strategy_choice = DynamicSoilMoisture()  # simulate soil moisture, don't just prescribe it -- compared against Sws
leaf_convection_model_choice = ElaborateLeafConvection()  # or SimpleLeafConvection()
interception_model_choice = LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2") #NoInterception()  # or LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")
canopy_air_profile_model_choice = RaupachLTheoryAirProfile(; far_field_mode=Val(:exact), near_field_subdivisions=20,
    relaxation=0.7, max_air_temperature_deviation=40.0u"K", aitken_omega_max=0.8, min_ground_resistance=5.0u"s/m")#RaupachLTheoryAirProfile(; far_field_mode=Val(:bulk))#RaupachLTheoryAirProfile(; far_field_mode=Val(:exact))#KTheoryAirProfile()
longwave_model_choice = LayeredRadiosityExchange() # LayeredLongwaveExchange(), LayeredRadiosityExchange(), AllPairsLongwaveExchange()
canopy_mode_choice = :full  # or :legacy
# Raised from the 0.02 m/s default -- near-calm nights were driving canopy_top_flux_boundary's ±40K clamp.
boundary_layer_model_choice = MoninObukhov(; min_friction_velocity=0.1u"m/s")
# canopy_wind_model_choice = MixingLengthCanopyWindAttenuation(;
#                                         shelter_floor = 0.003,
#                                         shelter_pai_coefficient = 0.1,
#                                         mixing_length_coefficient = 2.0,
#                                         mixing_length_pai_coefficient = 0.25,
#                                         )
canopy_wind_model_choice = ExponentialCanopyWindAttenuation(max_attenuation_coefficient = 2.879) #or MixingLengthCanopyWindAttenuation(; shelter_floor=..., ...)
#canopy_wind_model_choice = ExponentialCanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio(; ratio=0.5))
# ── "legacy" canopy_mode (see pipeline.jl prepare_site's canopy_mode kwarg) ──
# NoCanopy + a PAI-derived shade fraction + a wind-speed knockdown + a large
# horizon angle (sun only reaches the ground near-overhead) -- how forest
# sites were approximated before MultilayerCanopy existed. Run alongside
# canopy_mode=:full to quantify what the full canopy model actually buys.
# Per-site, not generic -- canopy openness varies enormously across these
# sites (Calperum's sparse mallee vs Cape Tribulation's closed rainforest),
# so a single wind_multiplier/horizon_angle would misrepresent one or the
# other. roughness_height stands in for the legacy model's Site.roughness_height
# (full mode always uses the bare-ground default of 0.004m instead -- see
# pipeline.jl -- since MultilayerCanopy's own air-profile model carries the
# canopy's aerodynamic roughness). Falls back to DEFAULT_LEGACY_PARAMS for any
# site without an entry.
const DEFAULT_LEGACY_PARAMS = (extinction_coefficient=0.5, wind_multiplier=0.5, horizon_angle=80.0u"°", roughness_height=0.004u"m")
const SITE_LEGACY_PARAMS = Dict(
    "CapeTribulation" => (extinction_coefficient=0.5, wind_multiplier=0.3, horizon_angle=85.0u"°", roughness_height=0.01u"m"),  # closed tropical rainforest canopy
    "Calperum"        => (extinction_coefficient=0.0, wind_multiplier=0.85, horizon_angle=10.0u"°", roughness_height=0.01u"m"),  # sparse, open mallee woodland
    "Whroo"           => (extinction_coefficient=1.0, wind_multiplier=0.5, horizon_angle=80.0u"°", roughness_height=0.01u"m"),  # dry sclerophyll woodland, intermediate
    "Wallaby"         => (extinction_coefficient=0.5, wind_multiplier=0.4, horizon_angle=70.0u"°", roughness_height=0.01u"m"),  # regrowth ash forest, denser than woodland
    "GWW"             => (extinction_coefficient=0.5, wind_multiplier=0.6, horizon_angle=25.0u"°", roughness_height=0.01u"m"),  # open eucalypt woodland
    "Longreach"       => (extinction_coefficient=0.0, wind_multiplier=1.0, horizon_angle=0.0u"°", roughness_height=0.01u"m"),   # grassland, effectively no canopy shelter
    "TiTreeEast"        => (extinction_coefficient=0.0, wind_multiplier=0.7, horizon_angle=10.0u"°", roughness_height=0.01u"m"),  # sparse mulga/spinifex woodland
    "AliceSpringsMulga" => (extinction_coefficient=0.0, wind_multiplier=0.8, horizon_angle=10.0u"°", roughness_height=0.01u"m"),
)
legacy_params(site_name) = get(SITE_LEGACY_PARAMS, site_name, DEFAULT_LEGACY_PARAMS)

# Depth/height grid constants (NMR_DEP_CM, SIM_DEPTHS_M, N_CANOPY_LAYERS):
# see utils.jl.

# ── Reference height for forcing (Ta/RH/Ws) ──────────────────────────────────
# Microclimate.jl always treats forcing as valid AT last(heights) (see
# pipeline.jl's _build_heights/resolve_site) -- i.e. whatever height sits at
# the top of the model's heights grid. The bare "Ta"/"Ws" tower variables are
# themselves QC-merged composites across several physical instrument heights
# (see each file's own `height` attribute + `description_L3`) which often do
# NOT match the site's `tower_height` global attribute:
#   CapeTribulation, GWW: bare Ta/Ws sit at tower_height already --
#     no entry needed, falls back to tower_height/bare "Ta"/"RH"/"Ws" below.
#   Calperum, Whroo: a single-sensor tower-top instrument exists
#     (Ws_SONIC_Av/Ta_SONIC_Av at Calperum, Ws_CSAT/Ta_HMP_36m at Whroo) --
#     SITE_FORCING_VARS below prefers that over the lower/blended composite.
#   Wallaby: no tower-top wind/temp sensor exists at all -- Ws_CSAT (= bare
#     Ws) sits at only 5m, well below both canopy top and
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
    "CapeTribulation" => 6.0,
    "Calperum"        => 0.3,
    "Whroo"           => 1.5,
    "Wallaby"         => 5.0,   # fire in 2009
    "GWW"             => 0.5,   # semi-arid eucalypt woodland, canopy_height=18m per file metadata
    "Longreach"       => 0.5,   # grassland, canopy_height=0.5m per file metadata
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
    "Calperum"        => :uniform,
    "Whroo"           => :bottom_heavy,
    "Wallaby"         => :uniform,
    "GWW"             => :uniform,  
    "Longreach"       => :uniform,       # grassland -- no crown structure
    "TiTreeEast"        => :top_heavy,  # sparse mulga woodland, simple shrub/small-tree structure
    "AliceSpringsMulga" => :top_heavy,
)

# ── Per-site canopy leaf spectral properties ──────────────────────────────────────────
# Guessed, no emprical data.
const SITE_LEAF_REFLECTANCE = Dict{String,Float64}(
    "CapeTribulation" => 0.25,
    "Calperum"        => 0.25,
    "Whroo"           => 0.25,
    "Wallaby"         => 0.15,
    "GWW"             => 0.30,
    "Longreach"       => 0.25,
    "TiTreeEast"        => 0.25,
    "AliceSpringsMulga" => 0.25,
)
const SITE_LEAF_TRANSMITTANCE = Dict{String,Float64}(
    "CapeTribulation" => 0.25,
    "Calperum"        => 0.15,
    "Whroo"           => 0.25,
    "Wallaby"         => 0.25,
    "GWW"             => 0.10,
    "Longreach"       => 0.25,
    "TiTreeEast"        => 0.15,
    "AliceSpringsMulga" => 0.15,
)

# ── Per-site leaf structural/physiological traits (Microclimate.jl's
# LeafParameters, fed into MultilayerCanopy) -- leaf_length/leaf_width feed
# HeatExchange.jl's boundary-layer convection, canopy_projection_ratio
# is Campbell's ellipsoidal `x` (1.0 spherical, 0.0 vertical, Inf horizontal).
# No site-specific literature values in hand yet, so every site starts at
# LeafParameters()'s own defaults -- free/tunable per site from here. Falls
# back to DEFAULT_LEAF_PARAMETERS for any site without an entry.
const DEFAULT_LEAF_PARAMETERS = (
    leaf_length=0.05u"m", leaf_width=0.02u"m", leaf_emissivity=0.97,
    leaf_water_potential=0.0u"J/kg", canopy_projection_ratio=1.0,
)
const SITE_LEAF_PARAMETERS = Dict(
    "CapeTribulation"   => DEFAULT_LEAF_PARAMETERS,
    "Calperum"          => DEFAULT_LEAF_PARAMETERS,
    "Whroo"             => DEFAULT_LEAF_PARAMETERS,
    "Wallaby"           => DEFAULT_LEAF_PARAMETERS,
    "GWW"               => DEFAULT_LEAF_PARAMETERS,
    "Longreach"         => DEFAULT_LEAF_PARAMETERS,
    "TiTreeEast"        => (leaf_length=0.02u"m", leaf_width=0.005u"m", leaf_emissivity=0.97,
                            leaf_water_potential=0.0u"J/kg", canopy_projection_ratio=1.0,),
    "AliceSpringsMulga" => (leaf_length=0.02u"m", leaf_width=0.005u"m", leaf_emissivity=0.97,
                            leaf_water_potential=0.0u"J/kg", canopy_projection_ratio=1.0,),
)
leaf_parameters(site_name) = LeafParameters(; get(SITE_LEAF_PARAMETERS, site_name, DEFAULT_LEAF_PARAMETERS)...)

# Per-site surface albedo -- ground cover/soil colour varies a lot across
# these sites (Calperum's pale sandy mallee soil vs Cape Tribulation's dark
# wet rainforest litter). Falls back to DEFAULT_ALBEDO for any site without
# an entry.
const DEFAULT_ALBEDO = 0.20
const SITE_ALBEDO = Dict(
    "CapeTribulation"   => 0.13,  # dark, wet closed-canopy rainforest litter
    "Calperum"          => 0.25,  # pale sandy mallee soil, sparse cover
    "Whroo"             => 0.10,  # box/ironbark woodland, leaf litter + some bare soil
    "Wallaby"           => 0.12,  # dense wet sclerophyll regrowth
    "GWW"               => 0.15,  # semi-arid eucalypt woodland, red sandy soil
    "Longreach"         => 0.25,  # dry grassland, pale soil/cured grass
    "TiTreeEast"        => 0.15,  # central Australian red sand, sparse mulga/spinifex
    "AliceSpringsMulga" => 0.15,
)

# Per-site critical soil (root) water potential at which stomata close --
# drought tolerance varies by species/vegetation type across these sites, but
# no site-specific literature values are in hand yet, so every site starts at
# the same free/tunable default as example_soil_hydraulic_model's own default
# (see Microclimate.jl/src/soil_hydraulics/campbell.jl). Falls back to
# DEFAULT_STOMATAL_CLOSURE_POTENTIAL for any site without an entry.
const DEFAULT_STOMATAL_CLOSURE_POTENTIAL = -1500.0u"J/kg"
const SITE_STOMATAL_CLOSURE_POTENTIAL = Dict{String,typeof(DEFAULT_STOMATAL_CLOSURE_POTENTIAL)}(
    "CapeTribulation"   => -1500.0u"J/kg",
    "Calperum"          => -2500.0u"J/kg",
    "Whroo"             => -1500.0u"J/kg",
    "Wallaby"           => -1500.0u"J/kg",
    "GWW"               => -1500.0u"J/kg",
    "Longreach"         => -2500.0u"J/kg",
    "TiTreeEast"        => -2500.0u"J/kg",
    "AliceSpringsMulga" => -2500.0u"J/kg",
)
stomatal_closure_potential(site_name) = get(SITE_STOMATAL_CLOSURE_POTENTIAL, site_name, DEFAULT_STOMATAL_CLOSURE_POTENTIAL)
soil_hydraulic_model(site_name) = Microclimate.example_soil_hydraulic_model(;
    stomatal_closure_potential = stomatal_closure_potential(site_name))

# Organic litter-layer soil override -- not universal (grassland/sparse
# mallee sites lack a real litter layer), off by default.
const SITE_ORGANIC_CAP = Dict(
    "CapeTribulation"   => true,
    "Calperum"          => false,
    "Whroo"             => false,
    "Wallaby"           => false,
    "GWW"               => false,
    "Longreach"         => false,
    "TiTreeEast"        => false,
    "AliceSpringsMulga" => false,
)
organic_cap(site_name) = get(SITE_ORGANIC_CAP, site_name, false)

# ── Soil source: real per-site SLGA fetch, a depth-flattened version of the
# same SLGA fetch, or a fixed literature texture class (Campbell & Norman
# 1998, Table 9.1) instead. Three choices:
#   :slga         -- real per-depth SLGA profile (default). Full resolution
#                     for Microclimate.jl; micropoint always gets this
#                     flattened anyway (single-slab solver -- see
#                     micropoint/ozflux/write_ozflux_micropoint_inputs.jl),
#                     so a :slga run is NOT like-for-like between the two
#                     models -- Microclimate.jl sees real depth variation,
#                     micropoint sees that same profile's depth-weighted mean.
#   :slga_uniform -- the same real SLGA fetch, collapsed to one depth-weighted
#                     mean slab (flatten_soil_profile, utils.jl) BEFORE either
#                     model runs -- both then see the exact same flat numbers,
#                     for a like-for-like comparison. Also avoids
#                     RunModelFull's indefinite hang on a full per-depth SLGA
#                     profile (bisected to Smax; see
#                     micropoint/ozflux/README.md) -- not that it matters
#                     here, since micropoint never saw the layered version.
#   a CAMPBELL_NORMAN_TEXTURES key -- fixed literature texture, ignores SLGA
#                     entirely. SLGA's real porosity can be low enough (some
#                     sites' Smax down around 0.15-0.3 vs the "Clay loam"
#                     placeholder's 0.46) to destabilize micropoint's solver
#                     regardless of :slga vs :slga_uniform -- this sidesteps
#                     that at the root when it's a problem.
# Per-site override; :slga is the default.
const CAMPBELL_NORMAN_TEXTURES = (
    sand             = (air_entry=0.7u"J/kg", b=1.7, Ksat=5.8e-3u"kg*s/m^3"),
    loamy_sand       = (air_entry=0.9u"J/kg", b=2.1, Ksat=1.7e-3u"kg*s/m^3"),
    sandy_loam       = (air_entry=1.5u"J/kg", b=3.1, Ksat=7.2e-4u"kg*s/m^3"),
    loam             = (air_entry=1.1u"J/kg", b=4.5, Ksat=3.7e-4u"kg*s/m^3"),
    silt_loam        = (air_entry=2.1u"J/kg", b=4.7, Ksat=1.9e-4u"kg*s/m^3"),
    sandy_clay_loam  = (air_entry=2.8u"J/kg", b=4.0, Ksat=1.2e-3u"kg*s/m^3"),
    clay_loam        = (air_entry=2.6u"J/kg", b=5.2, Ksat=6.4e-5u"kg*s/m^3"),
    silty_clay_loam  = (air_entry=3.3u"J/kg", b=6.6, Ksat=4.2e-5u"kg*s/m^3"),
    sandy_clay       = (air_entry=2.9u"J/kg", b=6.0, Ksat=3.3e-5u"kg*s/m^3"),
    silty_clay       = (air_entry=3.4u"J/kg", b=7.9, Ksat=2.5e-5u"kg*s/m^3"),
    clay             = (air_entry=3.7u"J/kg", b=7.6, Ksat=1.7e-5u"kg*s/m^3"),
)

# micropoint's own createsoilc(soiltype=...) strings, one per texture above --
# used so the micropoint export can ask createsoilc for a fully
# self-consistent table (Vq/Vm/Vo/Smin/n included, not just Smax/b/psi_e/Ksat)
# instead of mixing a placeholder texture's structure with our own overrides.
const CAMPBELL_NORMAN_MICROPOINT_NAME = Dict(
    :sand => "Sand", :loamy_sand => "Loamy sand", :sandy_loam => "Sandy loam",
    :loam => "Loam", :silt_loam => "Silt loam", :sandy_clay_loam => "Sandy clay loam",
    :clay_loam => "Clay loam", :silty_clay_loam => "Silty clay loam",
    :sandy_clay => "Sandy clay", :silty_clay => "Silty clay", :clay => "Clay",
)

const SITE_SOIL_SOURCE = Dict{String,Symbol}(
    "CapeTribulation" => :clay_loam,
    "Calperum" => :sandy_loam,
    "Wallaby" => :clay_loam,
    "Whroo" => :clay_loam,
    "GWW" => :clay_loam,
    "Longreach" => :clay,
    "TiTreeEast" => :clay_loam,
    "AliceSpringsMulga" => :clay_loam,
)  # e.g. "Whroo" => :sandy_loam
soil_source(site_name) = get(SITE_SOIL_SOURCE, site_name, :slga)

# ── UTC offset (micropoint/ozflux forcing only -- comparisons/ozflux itself
# stays in local time throughout). OzFlux records local STANDARD time
# year-round (no DST jump); minutes, since ACST sits at a half hour (mainland
# zones: AEST=+600, ACST=+570, AWST=+480). Whroo verified directly against
# the raw 30-min index around a real DST transition date; the rest come from
# each file's own time_zone attribute (Wallaby has none -- inferred AEST
# from its Victorian lat/lon). ───────────────────────────────────────────────
const SITE_UTC_OFFSET_MINUTES = Dict(
    "CapeTribulation"   => 600,  # Australia/Brisbane
    "Calperum"          => 570,  # Australia/Adelaide
    "Whroo"             => 600,  # Australia/Melbourne
    "Wallaby"           => 600,  # no time_zone attribute; Victorian site
    "GWW"               => 480,  # Australia/Perth
    "Longreach"         => 600,  # Australia/Brisbane
    "TiTreeEast"        => 570,  # Australia/Darwin
    "AliceSpringsMulga" => 570,  # Australia/Darwin
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
    # canopy_height=6.5m -- sub-canopy-to-near-top mast profile.
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
