# config.jl — site/soil config for the micropoint point comparison
# (comparison.jl) and the microclimf-based grid/topo comparisons
# (grid_scaling.jl, topo_scaling.jl — micropoint is point-only, no grid
# mode). `MCF_*` globals feed the microclimf scripts; `MP_*` globals feed
# micropoint.
#
# Included after comparisons/scan_snotel/{utils,config,pipeline}.jl, so
# those globals — depths, saturation_moisture, campbell_b,
# air_entry_potential, sat_hydraulic_cond, albedo, emissivity,
# roughness_height, profile_heights — are already in scope.
#
# Not `const`: plain globals so re-`include` during REPL iteration works.

# ── Site + date range ─────────────────────────────────────────────────────────
# Any comparisons/scan_snotel site ID works (lat/lon/elevation/obs flow
# through prepare_site(site_num, ...)); soil calibration below stays fixed to
# site 329's values regardless of site chosen.
site_num  = 2197
sim_start = Date(2015, 1, 1)
sim_end   = Date(2015, 12, 31)

# Both models run from run_start (sim_start - spinup_years) through sim_end
# so initial conditions settle; stats/plots below cover sim_start..sim_end.
spinup_years = 1
run_start    = sim_start - Year(spinup_years)

plot_start = Date(2015, 1, 1)
plot_end   = Date(2015, 12, 31)

# ── Comparison depths ─────────────────────────────────────────────────────────
# Both R models re-run their point model per depth. Matches ClimaLand
# comparison's depths for consistency.
MCF_DEPTHS_CM = [5.0, 20.0, 50.0]
# Above-ground heights reuse scan_snotel/config.jl's `profile_heights`.

# ── Soil parameter translation (Campbell-family, microclimf grid/topo only) ───
# Fixed values (site 329) -- grid_scaling.jl/topo_scaling.jl still use these.
# The point comparison (write_micropoint_inputs.jl) uses prep's own per-site
# SoilGrids values instead (see scan_snotel/pipeline.jl's prepare_site).
#   Smax = saturation_moisture (m^3/m^3)
#   b    = campbell_b (dimensionless)
#   psi_e [J/m^3] = air_entry_potential [J/kg] * rho_water [1000 kg/m^3]
#   gref/groundr = albedo (shortwave reflectance)
#   groundem = emissivity (longwave)
# Ksat is a best-effort conversion, not verified: R models want
# "kg/m^3/day" (Campbell 1974), Microclimate.jl's sat_hydraulic_cond is
# "kg*s/m^3" — not obviously the same dimensional family. Uses g/rho scaling
# (kg*s/m^3 -> m/s via *9.81/1000) then seconds->days. Check here first if
# soil moisture drainage looks wrong.
MCF_SMAX    = ustrip(saturation_moisture)
MCF_B       = campbell_b
MCF_PSI_E   = ustrip(u"J/kg", air_entry_potential) * 1000.0          # J/m^3
MCF_KSAT    = ustrip(sat_hydraulic_cond) * 9.81 / 1000.0 * 86400.0   # kg/m^3/day, best-effort
MCF_GROUNDR = albedo

# ── micropoint soil column structure + bare-ground surface ────────────────────
# Soil column depth/layers aren't in Microclimate.jl's own parameters — set
# to the micropoint vignette's own example values (2 m comfortably covers
# our deepest comparison depth, 50 cm; createsoilc internally multiplies by
# 1.5 for a boundary layer, so the real column goes to 3 m).
MP_NLAYERS    = 15
MP_TOTALDEPTH = 2.0   # m
MP_GROUNDEM   = emissivity
MP_RHO_KGM3   = ustrip(u"kg/m^3", uconvert(u"kg/m^3", bulk_density))

# micropoint's createsoilc uses different units than microclimf's soilc for
# two of the shared parameters (confirmed via `?createsoilc`, not guessed):
#   psi_e: air-entry water potential in METRES (head), not J/m^3 -- convert
#     J/kg -> m via /g (9.81 m/s^2), since J/kg = m^2/s^2.
#   Ksat: "kg s / m^3" -- an exact unit match to Microclimate.jl's own
#     sat_hydraulic_cond, so no conversion at all (unlike microclimf's
#     MCF_KSAT, which needs the kg/m^3/day best-effort conversion above).
MP_PSI_E = ustrip(u"J/kg", air_entry_potential) / 9.81   # m
MP_KSAT  = ustrip(sat_hydraulic_cond)                     # kg*s/m^3, exact

# Bare ground (vegp = NA): micropoint's vegetated mode needs a full
# stomatal/photosynthesis parameter set (Vcmx25, Tlow, Tup, alpha, Dcrit,
# Kxmx, psi50, ...) with no Microclimate.jl equivalent and no calibration
# data for this site — fabricating it would be pure invention. zm (bare
# surface roughness length) maps onto Microclimate.jl's roughness_height.
MP_ZM = ustrip(u"m", roughness_height)

# ── Output directory (Julia writes CSVs, R reads/writes back) ─────────────────
# Per-site function, not a fixed global -- comparison.jl runs many sites, each
# needs its own dir (also required for concurrent Rscript runs across sites,
# which a shared mutable global would race on). grid_scaling.jl/topo_scaling.jl
# still use the single default site_num.
mcf_outdir(site_num) = joinpath(@__DIR__, "outputs", "run_$(site_num)_$(sim_start)_$(sim_end)")
MCF_OUTDIR = mcf_outdir(site_num)   # single-site default, used by grid/topo scripts

println("micropoint config: default site $site_num, $run_start to $sim_end run " *
        "($spinup_years yr spin-up), reporting $sim_start to $sim_end, " *
        "depths=$MCF_DEPTHS_CM cm, bare ground")
