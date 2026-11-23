# ozflux

Runs Microclimate.jl at the TERN OzFlux flux towers ([ozflux.org.au](https://ozflux.org.au/)) and scores it against their observations.

## Setup

- `./download_ozflux.sh` downloads every site's OzFlux L3 file into `data/`.
- BARRA-R2, MODIS leaf area and SRTM elevation download on first use, into `data/`.
- SLGA soil needs a TERN API key in `RASTERDATASOURCES_SLGA_TERN_API_KEY`; SILO forcing (`mapper.jl`) needs `SILO_EMAIL`.

## Run

```
julia comparison.jl    # every site, stats and plots into the outputs folder
julia single_site.jl   # one site, for debugging
julia site_plots.jl    # cross-site plots from comparison.jl's all_stats.csv
```

## Folders

- `settings/` — site settings, data corrections, OzFlux naming conventions, land-cover classes, soil textures
- `data/` — downloaded inputs (not tracked)
- `outputs/<run>/` — results (not tracked)

## Environment variables

| Variable | Default | Choice |
|---|---|---|
| `OZFLUX_OUTPUTS_DIR` | `outputs/default` | output folder |
| `OZFLUX_AIR_PROFILE` | `exact` | in-canopy air profile: `exact` or `bulk` Raupach, `ktheory` |
| `OZFLUX_CANOPY_TOP_ENHANCEMENT` | `6.0` | canopy-top mixing enhancement |
| `OZFLUX_GROUND_STABILITY` | `louis` | ground stability: `louis`, `richardson`, `obukhov` |
| `OZFLUX_FREE_CONVECTION` | `0.04` | free-convection velocity floor, m/s |
| `OZFLUX_GROUND_TEMPERATURE` | `soil` | ground temperature in the canopy solve: `soil` or `hour_start` |
| `OZFLUX_GROUND_ROUGHNESS` | from `settings/sites.toml` | ground roughness height for every site, m |
| `OZFLUX_EXTRA_SOIL_DEPTHS` | none | extra soil nodes, m, comma separated |
| `OZFLUX_LONGWAVE` | `radiosity` | longwave: `radiosity`, `layered`, `black_body` |
| `OZFLUX_LITTER` | `site` | litter: `site` or `none` |
| `OZFLUX_LITTER_MOISTURE` | `dry` | litter moisture: `dry` or `wet` |
| `MICROCLIMATETESTS_RESCORE` | `false` | rescore saved results without solving |
