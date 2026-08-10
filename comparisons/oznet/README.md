# oznet

Validates MicroclimateMapper.jl's soil temperature/moisture output against
real OzNet Murrumbidgee catchment soil-monitoring observations (Smith et al.
2012, [oznet.org.au](http://www.oznet.org.au/)), with NicheMapR run on
identical forcing as a secondary comparison. This is the same dataset used
in the original NicheMapR-based study (Kearney & Maino 2018, *J. Hydrology*)
— the R prototype from that study lives alongside this folder
(`Australia_Soil_Test.Rmd`, `compare.R`, `plot_moist.R`, `plot_temp.R`,
`soilprops.txt`) for reference/provenance; this Julia pipeline doesn't read
any of them (see "Why not soilprops.txt" below).

`config.jl`'s `weather_source_choice` controls which forcing dataset drives
both models (`SILO`, analogous to the original study's AWAP grids —
`BARRA{BARRAC2,AUST04}` — the default, hourly ~4 km reanalysis — or
`NCEP{SurfaceFlux}`, analogous to the original study's `micro_ncep`). SILO and
BARRA are fetched via true OPeNDAP point queries
(`MicroclimateMapper.PointQuery`, toggled by `use_opendap_points` in
`config.jl`) rather than downloading/cropping a whole raster per site.

## Files

- `config.jl` — weather-source/pedotransfer-model choices, shared by both entry points.
- `utils.jl` — generic stats/interpolation helpers, plus `read_oznet_obs` (the OzNet file-format parser).
- `pipeline.jl` — the actual mechanics (fetch forcing, fetch SLGA soil texture, prepare a site, run Julia, run NicheMapR, report results).
- `comparison.jl` — entry point: runs a batch of sites, writes stats/plots, prints a cross-site summary.
- `single_site.jl` — entry point: runs one site with no try/catch, for stepping through in the debugger.
- `single_site.R` — same idea on the R side: sets `outdir` and sources `run_nmr.R` directly in an interactive R session, bypassing the output-swallowing you get from `run_nmr_batch!`'s concurrent `Rscript` subprocesses.
- `run_nmr.R` — reads the CSVs `pipeline.jl` writes and calls `NicheMapR::microclimate()` directly.

## How this differs from scan_snotel/

- **Depths are fixed, and shared with NicheMapR.** Both models run on the
  same 19-node scheme (`NMR_DEP19_CM` in `pipeline.jl`): NicheMapR's classic
  10-node `[0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]` cm grid plus an
  arithmetic midpoint between each pair, matching NicheMapR's own internal
  19-node hydraulic scheme (`pedotransfer()`/`BD[seq(1,19,2)]`). Julia's
  `MicroModel` runs on all 19 nodes; NicheMapR's Fortran solver still only
  takes the 10 real (odd-indexed) nodes. Comparison stats use the *nearest*
  node to each observation depth on both sides.
- **Soil profile is a live fetch, not a literature constant.** Each site's
  hydraulic profile comes from `build_soil_profile(SLGA, area; depths,
  pedotransfer_model)` (Cosby et al. 1984, matching the original study) at a
  small buffer around the site's coordinates, not a hand-specified constant
  profile.
- **Organic litter/crust cap on the top ~2.5 cm.** `prepare_site` overrides
  the SLGA-derived mineral thermal properties on nodes 1:3 (0, 1.25, 2.5 cm)
  with low-conductivity/high-heat-capacity litter values, matching
  NicheMapR's own default `cap=1` behaviour (`micro_global.R`) — without it,
  both models show unrealistically large diurnal temperature swings at 4 and
  15 cm. This is set once, in Julia; `write_nmr_inputs` writes Julia's
  already-capped profile straight into `nmr_soil.csv`, so NicheMapR picks it
  up automatically rather than needing its own separate override.
- **No snow model** (`NoSnow()`) — the Murrumbidgee catchment doesn't see
  persistent snow, so `run_nmr.R` runs with `snowmodel = 0`, which uses a
  simpler `soilinit` layout than scan_snotel's snow-enabled one (see the
  comment in `run_nmr.R`).

### Why not soilprops.txt

The original study's `soilprops.txt` is a static SLGA snapshot from
2014-08-01, pre-extracted for only 35 of the 38 OzNet sites in
`oznetsiteinfo2.csv`. Since this repo already has a live SLGA fetch path
(`build_soil_profile`, used elsewhere in `datasources/soiltexture.jl` and
`demos/barra.jl`), using it here covers all 38 sites uniformly with current
data rather than reproducing the original study's exact historical numbers.

## Data (not committed)

The OzNet dataset lives at `OZNET_DATA_DIR` in `config.jl`
(`C:\Users\mrke\Dropbox\Datasets\oznet`), not in this repo:

- `oznetsiteinfo2.csv` — per-site metadata: name, lon/lat, data filename,
  sampling interval, observation depths, and the 10-node SMDEP/TDEP
  simulation depth lists.
- `<site>_<freq>min_sm.txt` — one file per site (2 header lines + whitespace-
  delimited data; see `read_oznet_obs` in `utils.jl` for the exact format,
  including the one known site-level quirk, y3, where the metadata lists one
  fewer temperature depth than the file has columns for).

Everything this pipeline produces — `weather_cache/`, `julia_cache/`,
`julia_outputs/`, `nmr_outputs/` — is gitignored; delete any of them freely
to force a re-fetch/re-run (see the `reuse_*`/`cache_*` toggles in
`config.jl`).

## Usage

```
julia comparison.jl     # batch run, see site_subset/auto_date_range in comparison.jl
julia single_site.jl     # single site, full stack trace on failure
```

## References

Kearney, M. R., & Maino, J. L. (2018). Can next-generation soil data
products improve soil moisture modelling at the continental scale? An
assessment using a new microclimate package for the R programming
environment. *Journal of Hydrology*, 561, 662–673.

Smith, A. B., Walker, J. P., Western, A. W., Young, R. I., Ellett, K. M.,
Pipunic, R. C., … Richter, H. (2012). The Murrumbidgee soil moisture
monitoring network data set. *Water Resources Research*, 48(7).
