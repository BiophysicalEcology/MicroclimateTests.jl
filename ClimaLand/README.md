# ClimaLand

Validates CliMA's [ClimaLand.jl](https://github.com/CliMA/ClimaLand.jl)
(`Soil.EnergyHydrology`, standalone soil mode) against Microclimate.jl —
driven via MicroclimateMapper.jl/GRIDMET, the same way
[`comparisons/scan_snotel`](../comparisons/scan_snotel) does — at SNOTEL
site 329 (Beaver Dams, Utah), against real SCAN/SNOTEL observations.

## Scope: soil only, not canopy

Two prior, unfinished efforts live in
`C:\Users\mrke\Dropbox\Current Research Projects\julia_projects\`:

- `NMR_ClimaLSM_comparison/` — NicheMapR vs (presumably) ClimaLSM at a
  FLUXNET tower site (Ozark, half-hourly real forcing + LAI), comparing
  energy fluxes (LE, H, G, Rn) and soil T/moisture. That needs a canopy
  model with leaf gas exchange on both sides — Microclimate.jl doesn't do
  canopy energy balance or photosynthesis, so it isn't a fair
  "Microclimate.jl vs ClimaLand" test. Kept as background only, not ported.
- `ClimaLand/scan329_comparison.jl` — the real precedent for this folder:
  `Soil.EnergyHydrology` vs Microclimate.jl at this same site, soil-only.
  Its Campbell↔van Genuchten and de Vries↔Johansen soil parameter
  translations are reused here (see `config.jl`) — they're already
  calibrated to this exact site's soil properties, which match
  `comparisons/scan_snotel/config.jl`'s values.

## What's different from the old script

1. **Forcing comes from MicroclimateMapper.jl**, not a hand-built one-off CSV.
   `comparison.jl` calls `comparisons/scan_snotel/pipeline.jl`'s
   `prepare_site()` directly (via `include`, not duplicated) to fetch GRIDMET
   forcing and assemble Microclimate.jl's `MicroProblem`.
2. **ClimaLand is driven off Microclimate.jl's own solved hourly output**,
   not an independent reimplementation of its internal diurnal-interpolation
   scheme. `Microclimate.solve()`'s result (`MicroResult`) already contains
   the exact hourly `reference_temperature`/`reference_humidity`/
   `reference_wind_speed`/`pressure`/`global_radiation`/`sky_temperature`
   series Microclimate.jl used internally — `forcing_climaland.jl` reads
   those straight out, so both models see literally the same forcing except
   for the rain/snow phase split (which Microclimate.jl doesn't expose as a
   realized hourly water flux, so that part is still derived by hand from
   the daily rainfall total + an hourly-temperature threshold split).
3. **Fixed a version-drift bug.** The old script's last logged run
   (`full_output.txt` in the Dropbox folder) crashed at `sol_cl.u` —
   `solve!(simulation)` mutates in place and returns `nothing` in the
   installed ClimaLand version (v0.16.1), so the old `sol_cl = solve!(...)`
   pattern doesn't work anymore. `run_climaland.jl` uses `ClimaDiagnostics`
   instead (`DictWriter` + `default_diagnostics(..., output_vars=["tsoil",
   "swc", "si"])`) to capture the hourly time series — short names confirmed
   directly from ClimaLand's source
   (`src/diagnostics/land_compute_methods.jl`): `"tsoil"` → `p.soil.T` (K),
   `"swc"` → `Y.soil.ϑ_l` (liquid water content, m³/m³), `"si"` →
   `Y.soil.θ_i` (ice content, m³/m³).

## Soil moisture gets equal billing with temperature

Every stats row and every plot panel comes in both a temperature and a
moisture version — this isn't a temperature comparison with moisture bolted
on. See `comparison.jl`'s stats table (`D*cm` / `WC*cm` pairs) and the 3×3
plot grid (row 1 = temperature, row 2 = moisture, at 5/20/50cm).

## Changing site

`config.jl`'s `site_num = 329` is a plain top-level variable — swap it for
any other `comparisons/scan_snotel` site (any ID with a row in `Map metadata
export.csv` and a file in `observations/`) and everything else (lat/lon/
elevation, obs) flows through automatically via `prepare_site(site_num,
...)`. One caveat, inherited from `scan_snotel/config.jl` itself: the
Campbell soil parameters (and this folder's derived van Genuchten/Johansen
translation in `config.jl`) are one fixed profile calibrated to site 329,
applied globally regardless of which site you pick — `scan_snotel/config.jl`
already documents this as "adjust per site as needed."

## Spin-up

Both models actually run for `spinup_years` (default 1, in `config.jl`)
before `sim_start`, so initial conditions have time to settle before the
reported year starts — `run_start = sim_start - Year(spinup_years)` is what
gets passed to `prepare_site`, not `sim_start` itself. Every stat, plot, and
saved CSV still only covers `sim_start`..`sim_end`, the requested year; the
spin-up period is solved but never reported. `plot_start`/`plot_end` (also
in `config.jl`) narrow the plotted range further, within that reported year.

## Setup

ClimaLand's dependency stack (ClimaCore, ClimaTimeSteppers, ClimaParams,
ClimaUtilities, ClimaDiagnostics) isn't in the shared
`c:/git/BiophysicalEcologyEnv` used by every other script in this repo. It's
installed into your **global** Julia environment instead — Julia stacks the
global environment underneath whatever `--project=` is active, so this works
without touching BiophysicalEcologyEnv's own `Project.toml`/`Manifest.toml`:

```
julia                      # no --project flag
julia> include("setup.jl")
```

Then run the comparison the normal way for this repo:

```
cd c:/git/BiophysicalEcologyEnv
julia --project=. c:/git/MicroclimateTests.jl/ClimaLand/comparison.jl
```

If `using ClimaLand` errors under `--project=BiophysicalEcologyEnv` (a
version conflict between ClimaLand's dependency tree and something already
pinned there), fall back to an isolated environment scoped to this folder
instead — see the comment at the top of `setup.jl`.

## Files

- `setup.jl` — one-time global ClimaLand install.
- `config.jl` — site/date range, ClimaLand-specific parameter translation
  (van Genuchten, Johansen). Soil/snow calibration itself lives in
  `comparisons/scan_snotel/config.jl` (included by `comparison.jl`) — this
  file only adds ClimaLand-side settings on top.
- `forcing_climaland.jl` — builds ClimaLand's `TimeVaryingInput`s from
  Microclimate.jl's solved hourly output.
- `run_climaland.jl` — builds and solves `Soil.EnergyHydrology`, extracts
  results via `ClimaDiagnostics`.
- `comparison.jl` — entry point: runs both models, prints/saves stats,
  produces the comparison plot.
- `outputs/` — stats CSV + comparison PNG (gitignored, like every other
  comparison's outputs in this repo).

## Known risk areas (unverified — I couldn't execute Julia to confirm)

- `Soil.EnergyHydrologyParameters`/`AtmosDrivenFluxBC`/`WaterHeatBC` keyword
  names in `run_climaland.jl` are lifted from the old (possibly
  older-API-targeting) `scan329_comparison.jl` — if the constructor rejects
  a keyword, check `methods(Soil.EnergyHydrologyParameters)` /
  `?Soil.EnergyHydrology` interactively against the installed v0.16.1.
- `run_climaland.jl`'s `unpack()` assumes `DictWriter`'s dictionary keys
  convert to `Float64` seconds via `Float64(t)`. If ClimaLand's diagnostics
  use `ClimaUtilities.ITime` for `t` instead of plain numbers, this may need
  adjusting — see the comment right above `unpack()`.
