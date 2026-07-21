# micropoint

Validates Microclimate.jl (driven via MicroclimateMapper.jl/GRIDMET, same as
[`comparisons/scan_snotel`](../comparisons/scan_snotel)) against
[micropoint](https://github.com/ilyamaclean/micropoint) — Ilya Maclean's R
package for point-based mechanistic microclimate modelling — at SNOTEL site
1081 by default, against real SCAN/SNOTEL observations. Covers soil
temperature, soil moisture, downward shortwave/longwave radiation, and air
temperature at two heights.

This folder also keeps a separate, secondary set of grid-mode comparisons
(`grid_scaling.jl`, `topo_scaling.jl`) against
[microclimf](https://github.com/ilyamaclean/microclimf) — a related but
distinct package by the same author — since micropoint is point-only and has
no grid/DTM mode. See "Grid/topo comparisons" below.

## Changing site

`config.jl`'s `site_num` is a plain top-level variable — swap it for any
other `comparisons/scan_snotel` site and lat/lon/elevation/obs all flow
through automatically via `prepare_site(site_num, ...)`. The Campbell soil
parameters (and this folder's derived micropoint/microclimf translations)
are one fixed profile calibrated to site 329, applied globally regardless of
which site you pick — see `scan_snotel/config.jl`'s own comment.

## Spin-up

Both models run for `spinup_years` (default 1, in `config.jl`) before
`sim_start`, so initial conditions settle before the reported year starts —
`run_start = sim_start - Year(spinup_years)` is what gets passed to
`prepare_site`, not `sim_start` itself. Every stat, plot, and saved CSV
still only covers `sim_start`..`sim_end`; the spin-up period is solved but
never reported. `plot_start`/`plot_end` narrow the plotted range further,
within that reported year.

## How this works

Cross-language pattern matches this repo's existing NicheMapR integration
(`scan_snotel/run_nmr.R` + `pipeline.jl`'s `write_nmr_inputs`/
`run_nmr_batch!`): Julia writes `climdata.csv`/`params.csv`
(`write_micropoint_inputs.jl`), `Rscript run_micropoint.R` runs the model,
Julia reads `micropoint_out.csv` back (`comparison.jl`).

`climdata` (micropoint's hourly forcing input) is built from Microclimate.jl's
own solved hourly output (`MicroResult`), the same forcing-sharing approach
`ClimaLand/forcing_climaland.jl` uses — both models see the same air
temperature/humidity/wind/pressure/shortwave/sky-temperature series. One
consequence: the shortwave/longwave comparison in this folder is really
"does micropoint's (bare-ground) radiative transfer reproduce the input it's
given", not an independent test of atmospheric radiation physics — both
models' *downward* radiation at the site is ultimately driven from the same
source.

micropoint runs in bare-ground mode (`vegp = NA`) rather than with an
invented canopy — see `config.jl`'s comment on why (its vegetated mode needs
a full stomatal/photosynthesis parameter set with no Microclimate.jl
equivalent and no calibration data for this site).

`run_micropoint.R` also gives micropoint's soil column the same organic
top layer Microclimate.jl has (`scan_snotel/config.jl`'s depth-varying
`mineral_conductivity`/`mineral_heat_capacity`, low-conductivity/
high-heat-capacity above ~2.5cm): the top nodes' `Vq`/`Vm`/`Vo` composition
is set to (near-)pure organic, tapering to texture-default mineral by ~5cm,
calibrated against micropoint's own `thermalConductivityCpp` formula (read
from source, not guessed) to land close to Microclimate.jl's target
conductivities — see that script's comment for the exact reasoning and its
limits (deeper "mineral" nodes aren't forced to numerically match
Microclimate.jl's 2.5 W/m/K, since Clay loam's own composition doesn't
reach that value under micropoint's formula without an implausible organic
fraction — it's a structural match, not an exact one).

## Known model-capability mismatches (not bugs)

- **No snow submodel.** micropoint has no equivalent of microclimf's
  `runsnowmodel` — snow depth/SWE aren't compared here at all.
- **Soil moisture**: micropoint's soil moisture submodel returns a single
  bulk volumetric water content at the requested depth, not a depth-resolved
  profile like Microclimate.jl's per-node moisture. The comparison uses that
  value against Microclimate.jl's 5cm node.
- **Below-ground re-runs**: micropoint's point model is re-run from scratch
  per height/depth (`run_micropoint.R` does this for 2 heights + 5/20/50cm).
- **psi_e/Ksat/rho units differ from microclimf's soilc** — verified against
  `?createsoilc` directly (not guessed): psi_e is metres of head (not
  J/m^3), Ksat is `kg s / m^3` (an exact match to Microclimate.jl's
  `sat_hydraulic_cond`, so no conversion needed — contrast microclimf's
  best-effort Ksat conversion), rho is `kg/m^3` (not Mg/m^3). See
  `config.jl`'s comment.
- **Soil composition fractions (Vq/Vm/Vo/Mc)** aren't translated from
  Microclimate.jl at all — it has no equivalent parameters. `run_micropoint.R`
  uses `createsoilc`'s texture-based defaults for these.
- **`RunMicro`'s per-layer soil vectors** (Smax/b/psi_e/Ksat/rho) are
  replicated across all `nlayers + 1` nodes in `run_micropoint.R` — verified
  via `str(createsoilc(...))` that these are vectors, not scalars.

## Setup

micropoint isn't on CRAN (`remotes::install_github` only):

```
Rscript install_micropoint.R
```

Then run the comparison:

```
cd c:/git/BiophysicalEcologyEnv
julia --project=. c:/git/MicroclimateTests.jl/micropoint/comparison.jl
```

## Files

- `install_micropoint.R` — one-time R package install.
- `config.jl` — site/date range, comparison depths/heights, the
  Campbell-to-Campbell soil parameter translation shared by micropoint and
  microclimf, micropoint's bare-ground/soil-column settings.
- `write_micropoint_inputs.jl` — builds `climdata.csv`/`params.csv` from
  Microclimate.jl's solved output.
- `run_micropoint.R` — builds `soilc` (bare ground, no `vegp`), runs
  `RunMicro` per height/depth, writes `micropoint_out.csv` +
  `micropoint_timing.csv`.
- `comparison.jl` — entry point: runs both models, prints/saves stats
  (temperature, moisture, radiation, air temp) and a timing comparison,
  produces the comparison plot.
- `outputs/` — stats/timing CSVs + comparison PNG (gitignored).

## Grid/topo comparisons (microclimf, not micropoint)

`grid_scaling.jl` / `run_microclimf_grid.R` and `topo_scaling.jl` /
`run_microclimf_topo.R` compare MicroclimateMapper.jl's raster mode against
microclimf's native grid mode (`runmicro`) — solve-time scaling and, for
`topo_scaling.jl`, a same-forcing topographic-correction comparison on a
real DEM. These use microclimf, not micropoint, since micropoint is
point-only. Install separately:
`Rscript -e 'remotes::install_github("ilyamaclean/microclimf")'`.

`topo_scaling.jl`'s R stage (`run_microclimf_topo.R`) currently crashes with
a native access violation inside `runpointmodela`, likely from
`resampleclimdata`'s internal aggregation producing a degenerate small grid
at extreme coarse-to-fine resolution ratios — unresolved.
