# config.jl — ClimaLand-specific configuration.
#
# Included *after* `comparisons/scan_snotel/config.jl` (see comparison.jl), so
# the globals set there -- `depths`, `mineral_density`, `snow_temp_threshold`,
# `emissivity`, `roughness_height` -- are already in scope. K_sat/porosity are
# per-site instead (prep.site_sat_hydraulic_cond/site_saturation_moisture, set
# in run_climaland.jl); everything else here is a fixed modeling choice, not
# derived from the site.

# ── Site + date range ─────────────────────────────────────────────────────────
site_num  = 329
sim_start = Date(2015, 1, 1)
sim_end   = Date(2017, 12, 31)

# Both models actually run from run_start (sim_start - spinup_years) through
# sim_end, so initial conditions have a year to settle -- but every stat/plot
# below only covers sim_start..sim_end, the requested year.
spinup_years = 1
run_start    = sim_start - Year(spinup_years)

plot_start = nothing   # nothing = full sim_start..sim_end
plot_end   = nothing

# ── ClimaLand floating point precision ────────────────────────────────────────
const FT_CL = Float64   # matches Microclimate.jl's Float64 outputs.

# ── ClimaLand soil domain ─────────────────────────────────────────────────────
# 40 elements over 2 m -> 5 cm per cell, matching the deepest Microclimate.jl
# comparison depth (200 cm) used in comparisons/scan_snotel.
const ZMIN_CL  = FT_CL(-2.0)
const ZMAX_CL  = FT_CL(0.0)
const NELEM_CL = 40

# ── van Genuchten retention curve shape (ClimaLand) ────────────────────────────
# alpha/n fixed at site 329's hand-derived values -- NOT a general Campbell(b,
# psi_e) -> van Genuchten conversion, and not currently re-derived per site.
# A numerical curve match (matching Campbell's own retention curve) does not
# reproduce these values, so the original derivation likely used a different
# method (e.g. a texture-based Rawls & Brakensiek regression) -- unconfirmed.
# K_sat/porosity (ν) ARE per-site now -- see run_climaland.jl, which reads
# prep.site_sat_hydraulic_cond/site_saturation_moisture directly.
const VG_ALPHA = FT_CL(1.22)    # m^-1
const VG_N     = FT_CL(1.56)
const THETA_R  = FT_CL(0.04)
const S_S_CL   = FT_CL(1e-3)    # m^-1

# Soil thermal (Johansen, ClimaLand) approximated to match the layered
# mineral_conductivity profile from scan_snotel/config.jl: low k (~0.2 W/m/K,
# organic) in the top 5 cm, high k (~2.5 W/m/K, mineral) below.
const NU_SS_OM_TOP     = FT_CL(0.30)
const NU_SS_QUARTZ_TOP = FT_CL(0.02)
const NU_SS_OM_DEEP    = FT_CL(0.01)
const NU_SS_QUARTZ_DEEP = FT_CL(0.45)
const NU_SS_GRAVEL     = FT_CL(0.0)
const THERMAL_TRANSITION_DEPTH_M = 0.05   # organic/mineral transition, matches scan_snotel's mineral_conductivity break

println("ClimaLand config: site $site_num, $run_start to $sim_end run " *
        "($spinup_years yr spin-up), reporting $sim_start to $sim_end, " *
        "vG(alpha=$VG_ALPHA, n=$VG_N)")
