# ozflux

Validates the canopy model (`MultilayerCanopy`) against real TERN/OzFlux
eddy-covariance flux-tower observations ([ozflux.org.au](https://ozflux.org.au/)):
net/upward radiation, sensible/latent/ground heat flux, soil
temperature/moisture at depth, and — where a site has multi-height sensors —
the shape of the wind/temperature profile through and above the canopy.

Unlike `oznet`/`scan_snotel`, forcing comes straight from each site's own
tower file (30-min, aggregated to hourly), not a fetched grid — this is a
canopy-model-focused pipeline, not a soil-moisture-grid one. A SILO
gridded-forcing run is also available per site-year (`run_site_silo`), both
as its own comparison against the tower obs (gridded vs point-observed skill)
and, via `run_site_gapfilled`, as a real-data gap-filler for missing tower
forcing (see "SILO forcing" below).

## Files

- `config.jl` — model/SILO/legacy-mode choices, per-site LAI and multi-height
  sensor lists, shared by both entry points.
- `utils.jl` — generic stats helpers (shared with `oznet`) plus the OzFlux
  NetCDF reader/aggregator (`read_ozflux_nc`, `discover_depth_series`,
  `discover_site_years`).
- `pipeline.jl` — the mechanics: `resolve_site` (metadata) ->
  `prepare_site`/`prepare_site_silo` (read/aggregate/build the
  `MicroProblem`) -> `run_site`/`run_site_silo`/`run_site_gapfilled` (solve).
- `report.jl` — stats + plots for every comparison target
  (`report_site_results`).
- `comparison.jl` — entry point: batch over a site x year(s) x forcing-mode x
  canopy-mode grid, writes stats/plots, prints a cross-site summary.
- `single_site.jl` — entry point: one site-year, no try/catch, for stepping
  through in the debugger.

## Data

`<Site>_<Year>_L3.nc` files in `OZFLUX_DATA_DIR` (`data/` by default) — one
OzFlux L3 NetCDF per site-year, downloaded from ozflux.org.au. Not committed
(gitignored). `discover_site_years` scans this directory, so `comparison.jl`
picks up whatever's actually on disk rather than a hand-maintained list.

Current sites: Cape Tribulation (tropical rainforest, canopy ~25 m, tower
45 m), Calperum (mallee, canopy ~3 m, tower 20 m, multi-height Ta/AH/Ws),
Whroo (dry sclerophyll woodland, canopy 28 m, tower 35 m, Ta/Ws at 1/2/4/8/16 m
sub-canopy and 32/36 m above — the richest in-canopy profile available), and
Wallaby (regrowth ash forest, canopy 8-10 m, tower 12 m, one sub-canopy Ta
sensor at 5 m).

Variable-naming conventions are **not** consistent across sites (e.g.
Calperum's `Ws_RMY2m_Av` vs Whroo's `Ws_RMY_2m`), which is why
`SITE_HEIGHT_SERIES` in `config.jl` is an explicit per-site list rather than
an auto-detected one.

### LAI placeholder

`SITE_LEAF_AREA_INDEX` in `config.jl` is a per-site, hand-supplied constant —
no OzFlux L3 file carries LAI in its metadata, and no phenology/dynamic LAI
is modelled. `resolve_site` errors for any site without an entry rather than
falling back to a generic literature guess, which would misrepresent the
comparison. Less of a concern at evergreen sites (Cape Tribulation), more so
anywhere with strong seasonal LAI variation.

## SILO forcing (secondary comparison + gap-fill)

`weather_source_choice = SILO` (`config.jl`) drives a second, independent run
per site-year via `MicroclimateMapper`'s usual point-fetch machinery (same
pattern as `oznet`). Two uses:

1. **`run_site_silo`** — the canopy model driven entirely by SILO's daily
   min/max grids (+ CRUCL2 monthly wind climatology, SILO's own wind
   fallback) instead of tower obs, evaluated against the same tower targets.
   An interesting comparison in its own right: how much skill is lost going
   from point-observed hourly forcing to gridded daily forcing?
2. **`run_site_gapfilled`** — fills missing tower `Ta`/`RH`/`Ws`/`Fsd` hours
   from a solved SILO run's own recorded per-hour forcing (`MicroResult`'s
   `reference_*`/`global_radiation` columns — the diel-synthesized values the
   solver actually used, not simulated fluxes), instead of `prepare_site`'s
   default plain linear interpolation across gaps.

**Caveats**: SILO has no wind-speed layer at all (CRUCL2 climatology
substitutes) and no downward-longwave layer either -- `Fld` gaps are instead
backed out of the SILO run's own solved sky_temperature (Stefan-Boltzmann:
incoming longwave = σ·T_sky⁴), a physically-grounded estimate rather than a
straight-line interpolation, but still model-derived, not observed. SILO is
daily-native — `environment_hourly` for a SILO fetch is a stub; the "hourly"
values used for gap-filling are diel-synthesized from the daily min/max, not
truly observed sub-daily dynamics.

## Legacy canopy_mode

`canopy_mode = :legacy` (`prepare_site`/`run_site`/`run_site_gapfilled`) runs
`NoCanopy()` with a PAI-derived shade fraction (Beer-Lambert,
`legacy_extinction_coefficient`), a wind-speed knockdown
(`legacy_wind_multiplier`), a large horizon angle
(`legacy_horizon_angle`, so direct beam only reaches the ground near
solar noon), and a ground-surface roughness height (`legacy_roughness_height`,
replacing the 0.004m bare-ground default `:full` uses) — approximating a
forest site the way it would have been done before `MultilayerCanopy` existed.
Run alongside `canopy_mode = :full` to see what the full canopy model
actually buys over that approximation.
`:legacy` has no below-canopy resolution and no canopy-summed Fh/Fe/Fsu/Flu/Fn
(`report.jl` skips those comparisons and falls back to the ground-level MOST
estimate for Fh only); it isn't wired up for the SILO forcing path.

## Comparison targets

Per site-year: forcing sanity plots (Ta/RH/Ws/Fsd/Fld/Precip, obs only);
`Fsu`/`Flu` (canopy-top boundary shortwave/longwave) and `Fn`; `Fh`/`Fe`
(canopy-summed sensible/latent heat) and `Fg` (Fourier's-law ground heat
flux); `Ts`/`Sws` at every observed depth; and, where `SITE_HEIGHT_SERIES`
has entries, Ta/Ws at each sensor height, routed to the canopy-resolved
output below `canopy_height` or the free-atmosphere MOST profile above it —
the sharpest available test of whether the vertical profile shape is right,
not just the aggregate flux.

Every comparison: `ModelStats` (r/RMSE/bias/n) via `compute_stats`, an
overlaid obs-vs-model timeseries plot, and a scatter/1:1 plot, saved under
`outputs_dir` (`julia_outputs/` by default, gitignored) per site-year and
forcing/canopy mode.

### Vertical profile snapshots

`report.jl`'s `plot_canopy_profiles_all(result, result.output, profile_times)`
(called from `single_site.jl`) plots height-vs-value profiles (Ta and Ws, one
panel per hour in `profile_times`) instead of a fixed-height time series: the
model's full profile as a line (canopy-resolved below `canopy_height`,
free-atmosphere MOST above), tower obs at each `SITE_HEIGHT_SERIES` sensor
height as points. Most useful at Whroo (5 real sub-canopy heights) and
Wallaby (1); saved under `outputs_dir/profiles_vertical/`.

## Known caveats

- **No in-canopy profile obs** at most sites (Whroo excepted) — validates
  aggregate flux partitioning only, not internal transport-scheme
  correctness; several different internal parameterizations could plausibly
  match the same aggregate flux (equifinality).
- **EC-tower energy-balance-closure gap** — measured `Fn - Fg` typically
  exceeds measured `Fh + Fe` by 10-30%. Microclimate's own balance closes by
  construction, so a residual mismatch against the tower's *unclosed*
  observations isn't necessarily model error.
- **No CO2/photosynthesis module** — `Fc` isn't usable.
- **Soil profile is a literature default** (`Microclimate.example_soil_profile`),
  not a per-site texture fetch — unlike `oznet`'s live SLGA pull, since this
  pipeline doesn't otherwise need `MicroclimateMapper`'s gridded machinery for
  the primary (tower-forced) run.
- **30-min OzFlux data is paired-averaged to hourly**; sub-hourly dynamics
  are lost. `MultilayerCanopy`/`MicroModel` are hourly-only in any case.

## Usage

```
julia comparison.jl     # batch run, see SITE_YEARS/FORCING_MODES/CANOPY_MODES in comparison.jl
julia single_site.jl     # single site-year, full stack trace on failure
```
