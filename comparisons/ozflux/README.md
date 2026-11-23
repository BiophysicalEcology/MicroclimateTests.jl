# ozflux

Validates the canopy model (`MultilayerCanopy`) against real TERN/OzFlux
eddy-covariance flux-tower observations ([ozflux.org.au](https://ozflux.org.au/)):
net/upward radiation, sensible/latent/ground heat flux, soil
temperature/moisture at depth, and — where a site has multi-height sensors —
the shape of the wind/temperature profile through and above the canopy.

Forcing comes from each site's highest above-canopy tower sensors (30-min,
aggregated to hourly), with BARRA-R2 reanalysis for whatever the tower can't
supply (`forcing_selection.jl`). A SILO gridded-forcing run is also available
per site-year (`run_site_silo`) as its own comparison against the tower obs.

## Files

- `config.jl` — paths and model algorithm choices, shared by both entry points.
- `settings/` — every fact about the datasets and sites (see "Settings").
- `utils.jl` — generic stats helpers (shared with `oznet`) plus the OzFlux
  NetCDF reader/aggregator (`read_ozflux_nc`, `discover_depth_series`,
  `discover_site_years`).
- `pipeline.jl` — the mechanics: `resolve_site` (metadata) ->
  `prepare_site` (read/aggregate/build the `MicroProblem`) -> `run_site`
  (solve).
- `mapper.jl` — everything that needs MicroclimateMapper: SILO forcing
  (`run_site_silo`) and SLGA soil (`soil_source` `:slga`/`:slga_uniform`).
- `report.jl` — stats + plots for every comparison target
  (`report_site_results`).
- `comparison.jl` — entry point: every site's best year (or whole record
  where shorter), writes stats/plots, prints a cross-site summary.
- `single_site.jl` — entry point: one site-year, no try/catch, for stepping
  through in the debugger.
- `site_plots.jl` — entry point: Makie plots of every site's stats from
  `all_stats.csv`, into `julia_outputs/site_plots`. Needs no solve.

## Data

`<Site>[<n>]_L3.nc` files in `OZFLUX_DATA_DIR` (`data/` by default) — one
OzFlux L3 NetCDF per site holding its whole record (TERN OzFlux releases,
e.g. `Whroo1_L3.nc`). `read_ozflux_nc` reads only the requested years. Not
committed (gitignored). `discover_site_years` scans this directory and each
file's time axis, so `comparison.jl` picks up whatever's actually on disk
rather than a hand-maintained list.

Current sites: Cape Tribulation (tropical rainforest, canopy ~25 m, tower
45 m), Calperum (mallee, canopy ~3 m, tower 20 m, multi-height Ta/AH/Ws),
Whroo (dry sclerophyll woodland, canopy 28 m, tower 35 m, Ta/Ws at 1/2/4/8/16 m
sub-canopy and 32/36 m above — the richest in-canopy profile available), and
Wallaby (regrowth ash forest, canopy 8-10 m, tower 12 m, one sub-canopy Ta
sensor at 5 m).

Variable-naming conventions are **not** consistent across sites (e.g.
Calperum's `Ws_RMY2m_Av` vs Whroo's `Ws_RMY_2m`); the naming rules the code
relies on are data, in `settings/ozflux_conventions.toml`.

### Settings

Every fact about the datasets and sites is data in `settings/`, never code:
OzFlux naming conventions, data corrections (mistakes only, each with its
evidence), per-site settings with their defaults, soil texture classes and peat (`peat.toml`).

### Leaf area index

MODIS MCD15A3H's median over the run year, unless `settings/sites.toml` gives a
site a measured value. No phenology/dynamic LAI is modelled.

## SILO forcing (secondary comparison)

**`run_site_silo`** runs the canopy model driven entirely by SILO's daily
min/max grids (+ CRUCL2 monthly wind climatology) via `MicroclimateMapper`,
evaluated against the same tower targets: how much skill is lost going from
point-observed hourly forcing to gridded daily forcing?

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
- **Soil is SLGA at the site's 90 m cell** (texture converted to USDA classes),
  except at sites whose file names a peat or Organosol, which SLGA misses (cited
  values in `settings/peat.toml`), and Whroo, with measured horizons. Downloads need
  a TERN API key in `ENV["RASTERDATASOURCES_SLGA_TERN_API_KEY"]`.
- **30-min OzFlux data is paired-averaged to hourly**; sub-hourly dynamics
  are lost. `MultilayerCanopy`/`MicroModel` are hourly-only in any case.

## Usage

```
julia comparison.jl     # batch run, see SITE_YEARS/FORCING_MODES/CANOPY_MODES in comparison.jl
julia single_site.jl     # single site-year, full stack trace on failure
julia site_plots.jl     # plots of every site's stats, after comparison.jl
```
