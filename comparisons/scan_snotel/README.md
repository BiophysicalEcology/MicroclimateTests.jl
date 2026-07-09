# scan_snotel

Validates MicroclimateMapper.jl's soil temperature/moisture and snow output
against real USDA NRCS SCAN/SNOTEL station observations, with NicheMapR run
on identical forcing as a secondary comparison. `config.jl`'s
`weather_source_choice` controls which forcing dataset drives both models
(`GRIDMET`, `ERA5`, or `NCEP{SurfaceFlux}` — see the comment there for why
CHELSA isn't included).

## Files

- `config.jl` — model/soil/snow parameters and weather-source choice, shared by both entry points.
- `utils.jl` — generic stats/interpolation/de-spiking helpers.
- `pipeline.jl` — the actual mechanics (fetch forcing, prepare a site, run Julia, run NicheMapR, report results). No site list or date range here — that's in the entry points.
- `comparison.jl` — entry point: runs a batch of sites, writes stats/plots, prints a cross-site summary.
- `debug_site.jl` — entry point: runs one site with no try/catch, for stepping through in the debugger.
- `run_nmr.R` — reads the CSVs `pipeline.jl` writes and calls `NicheMapR::microclimate()` directly.

## Data (not committed)

Two inputs are required and are gitignored — nothing here should end up
tracked in git:

- `Map metadata export.csv` — site ID/lat/lon/elevation/state/network. Export
  from the [NRCS National Water and Climate Center site list](https://wcc.sc.egov.usda.gov/nwcc/).
- `observations/<site_id>.csv` — one CSV per site with columns `DateTime,
  SNWD.I, WTEQ.I, STO.I_2, STO.I_4, STO.I_8, STO.I_20, STO.I_40, SMS.I_2,
  SMS.I_4, SMS.I_8, SMS.I_20, SMS.I_40` (NRCS SCAN/SNOTEL report codes for
  snow depth, SWE, and soil temperature/moisture at 2/4/8/20/40 inch depths).
  Pulled from the same NWCC reporting API — no fetch script exists yet; if
  you write one, it belongs in this folder (e.g. `fetch_observations.jl`) so
  it stays next to the pipeline that consumes its output.

Everything else this pipeline produces — `weather_cache/`, `julia_cache/`,
`julia_outputs/`, `nmr_outputs/` — is also gitignored; delete any of them
freely to force a re-fetch/re-run (see the `reuse_*`/`cache_*` toggles in
`config.jl`).

## Usage

```
julia comparison.jl     # batch run, see site_subset/auto_date_range in comparison.jl
julia debug_site.jl     # single site, full stack trace on failure
```
