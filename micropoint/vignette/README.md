# micropoint vignette comparison

Reproduces [micropoint](https://github.com/ilyamaclean/micropoint)'s own published
vignette forest example (`createvegp("BET.Te")`, `PAIgeometry`, `createsoilc("Clay
loam")`, `RunModelFull`, driven by the package's own built-in `climdata`) against an
equivalent Microclimate.jl `MultilayerCanopy` run — using micropoint's own canonical
example parameters, not a fabricated real-site calibration.

This is separate from the sibling [`micropoint/`](../README.md) comparison, which
validates **bare-ground** mode against real SNOTEL sites. That comparison explicitly
skips vegetated mode because micropoint's vegetated mode needs a full Eller-et-al.
photosynthesis/stomatal parameter set with no Microclimate.jl equivalent and no
calibration data for a real site. This comparison sidesteps that: it uses micropoint's
own built-in demo parameters (not a real site), so there's nothing to fabricate — and it
lets the two known-mismatched physics packages run side by side anyway, documented below,
since a full physics match was never the point.

Per-user scope: "a reasonable number of comparisons ... to check that similar results
are coming from them, and to benchmark as well" — not an exhaustive audit.

## How this works

Opposite data-flow direction from the sibling `micropoint/` comparison: there, Julia
writes inputs and R reads them. Here, R's own canned `climdata` is the source of truth,
so `run_micropoint_vignette.R` writes `climdata.csv` (plus every mappable parameter and
`RunModelFull`'s full output matrices), and `vignette_comparison.jl` reads it back to
drive an identical-forcing `MultilayerCanopy` run.

## Setup

```
Rscript install_micropoint.R   # (in ../micropoint, one-time, if not already installed)
```

Then:

```
cd c:/git/MicroclimateTests.jl/micropoint/vignette
Rscript run_micropoint_vignette.R outputs

cd c:/git/BiophysicalEcologyEnv
julia --project=. c:/git/MicroclimateTests.jl/micropoint/vignette/vignette_comparison.jl
```

## Parameter translation

Verified directly (installed micropoint 0.1.0, live runs) — not guessed:

| micropoint (`vegp`/`soilc`) | Microclimate.jl | Notes |
|---|---|---|
| `h` (25 m) | `MultilayerCanopy.canopy_height` | |
| `PAIgeometry(pai,7,70,n=20)` | `plant_area_index` (vector) | **bottom-to-top in R, top-to-bottom in Microclimate.jl** — reversed once at translation time |
| `x` (1.6) | `LeafParameters.leaf_angle_distribution_parameter` | |
| `lref`/`ltra` (0.4/0.2) | `TwoStreamRadiation.leaf_reflectance`/`leaf_transmittance` | |
| `len`/`wid` (0.12/0.05) | `LeafParameters.leaf_length`/`leaf_width` | assumed metres (cm would be an implausibly tiny leaf) |
| `vegem` (0.97) | `LeafParameters.leaf_emissivity` | |
| `mwft` (0.2) | `LayeredRainInterception.leaf_water_storage_capacity` | mm ≈ kg/m² (1:1), matching this model's own docstring convention |
| `soilc$psi_e` (metres head) | `air_entry_water_potential` | `J/kg = m × 9.81` — same conversion already verified in `../micropoint/config.jl` |
| `soilc$Ksat` (kg·s/m³) | `saturated_hydraulic_conductivity` | exact unit match, no conversion (same as `../micropoint/config.jl`) |
| `soilc$b` | `campbell_b_parameter` | |
| `soilc$gref`/`groundem` | `Site.albedo` / `environment_daily.surface_emissivity` | |

## Gotcha: `HourlyTimeseries.longwave_radiation` is not actually used by the canopy

`canopy_longwave!`'s sky term (`precompute_longwave_sky`) never reads
`environment_instant.longwave_radiation` — it always *derives* its own clear-sky
estimate (`CampbellNormanAtmosphericRadiation`) blended with a synthetic cloud term via
`cloud_cover`, regardless of what's supplied as measured downward longwave. A flat/wrong
`cloud_cover` therefore silently overrides the real sky condition every hour — most
consequentially on clear nights, where an overstated cloud fraction suppresses radiative
cooling that should otherwise happen (this produced a real, confusing bug during
development: leaf/air temperature never dropped fully to ambient at night). Fixed here
by inverting that same clear/cloud blend against climdata's real `lwdown` to back out an
equivalent per-hour `cloud_cover` (see `vignette_comparison.jl`) — not a perfect
substitute for a direct longwave input, but far better than a flat placeholder.

## Known capability mismatches (not bugs)

- **Stomatal conductance is the single biggest expected source of divergence.**
  Microclimate.jl uses `PrescribedStomatalConductance` (day/night gating only — the
  default, kept as-is). micropoint uses Eller et al.'s hydraulic-optimization
  photosynthesis model (`Vcmx25, Tup, Tlow, Dcrit, alpha, Kxmx, hv, f0, fd, psi50, apsi`)
  plus per-layer `Lfrac` (living-leaf fraction, feeds Eller's respiration term). None of
  this has a Microclimate.jl equivalent — expect real, structural differences in
  transpiration/latent heat dynamics, not just numerical noise.
- **No PAR-band radiation.** micropoint's `lrefp`/`ltrap`/`grefPAR` (photosynthetically-
  active-radiation-specific optical properties) have no Microclimate.jl equivalent —
  Microclimate.jl's two-stream model only handles broadband shortwave.
- **Soil composition unmapped.** `Smax`/`Smin`/`n`/`Vq`/`Vm`/`Vo`/`Mc`/`rho` have no
  Microclimate.jl equivalent — its Campbell hydraulics use a different retention-curve
  formulation (`psi_e`/`b`/`Ksat` only, no van-Genuchten-style `Smax`/`Smin`/`n`). Texture
  match is approximate.
- **Canopy layer heights are assumed, not verified.** `RunModelFull` doesn't return
  per-layer heights — both models here assume even spacing (`canopy_height × i/20`) from
  ground to canopy top. If micropoint's own internal spacing differs, the height axis on
  the profile plot will be systematically off (values would still pair up correctly
  layer-by-layer for the top/bottom time-series comparisons, which don't depend on the
  height assumption).
- **Diffuse-fraction handling differs.** micropoint's `climdata` supplies a direct
  measured `difrad`/`swdown` split. Microclimate.jl's `HourlyTimeseries` only takes total
  `global_radiation` and derives its own diffuse fraction internally
  (`AbstractDiffuseFractionModel`) — so even with identical total irradiance, the two
  models' direct/diffuse partition (and therefore canopy radiation penetration) can differ
  slightly.
- **Soil moisture not compared.** Kept out of scope here (surface/deepest *temperature*
  only) — the sibling `micropoint/README.md` already documents the same bulk-vs-profile
  soil moisture representation mismatch, which would apply here too.

## Files

- `run_micropoint_vignette.R` — builds micropoint's own vignette forest example
  (`createvegp("BET.Te")`, `PAIgeometry`, `createsoilc("Clay loam")`), runs
  `RunModelFull` on the package's built-in `climdata`, times it, writes every
  input/output as CSV.
- `vignette_comparison.jl` — reads those CSVs, builds an equivalent `MultilayerCanopy`
  run from the identical hourly forcing, times it, and compares: canopy-top/bottom
  air/leaf temperature, relative humidity, and wind speed (stats, direct hourly, no
  observations involved); surface/deepest soil temperature; a full height-profile plot
  (temperature, humidity, downward/upward shortwave+longwave, wind) at the year's
  peak-shortwave hour; a canopy-top air/leaf temperature time series over the
  surrounding week. Prints a timing comparison (R `RunModelFull` vs. Julia `solve`).
- `outputs/` — CSVs + comparison PNGs (gitignored).
