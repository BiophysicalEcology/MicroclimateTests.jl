# ozflux_below_canopy — micropoint (R) vs Microclimate.jl vs real OzFlux obs

Three-way comparison at Whroo: drives micropoint's vegetated
`RunModelFull` with the same real, gap-filled OzFlux tower forcing
`comparisons/ozflux/pipeline.jl`'s `run_site_gapfilled` already uses for
Microclimate.jl's `MultilayerCanopy`, then compares all three (tower/
sub-canopy sensor obs, Microclimate.jl, micropoint) for canopy air/wind/leaf
temperature, soil temperature, and soil moisture, plus a solve-time
comparison.

Julia writes input CSVs, `Rscript <script>.R <outdir>` runs as a
subprocess, Julia reads the output CSVs back.

## Files

- `write_ozflux_micropoint_inputs.jl` — writes `climdata.csv`, `params.csv`,
  `canopy_params.csv`, `canopy_layers.csv` from a solved `comparisons/ozflux`
  result.
- `run_micropoint_ozflux_vegetated.R` — reads those four CSVs, runs
  `RunModelFull` in chunks (see below), writes `veg_*.csv` + timing +
  `soil_depths.csv`.
- `ozflux_below_canopy_comparison.jl` — entry point: `run_site_gapfilled` ->
  write inputs -> `Rscript` -> read `veg_*.csv` back -> 3-way stats/plots +
  timing comparison. Run with `julia ozflux_below_canopy_comparison.jl`.
- `outputs/` — gitignored.

## Known limitation: partial-year coverage (read before trusting the stats)

`RunModelFull` can hit a permanent, non-recovering NaN cascade under real
(as opposed to smooth/synthetic) forcing — once triggered, every hour after
it in that call is NaN, with no recovery, confirmed deterministic (repeat
runs on identical input give identical results). Two mitigations are
already applied:

1. **Wind speed uses Microclimate.jl's own canopy-top wind**
   (`output.canopy.wind_speed[:,1]`), not the raw tower `Ws`. A/B tested
   directly (same rows, same R session): raw tower wind's real gust-to-gust
   variability cascades 453/744 hours in a one-month test; the
   canopy-attenuated series (same underlying meteorology, passed through
   the canopy's own momentum sink) drops that to 1/744. Ta/RH stay raw —
   tested and confirmed they don't need this (raw vs. canopy-top differ by
   <0.05°C / <0.3 percentage points at Whroo's PAI=0.7) — and `zref` stays
   at the true tower reference height (36m for Whroo), not canopy_height:
   switching Ta/RH to canopy-top too *and* setting `zref=canopy_height` for
   full height-consistency was tried and made the cascade dramatically
   *worse* (canopy_height as `zref` is a boundary case `RunModelFull`'s own
   docs treat differently from genuinely-above-canopy forcing) — reverted.
   The height mismatch this leaves for wind specifically is a smaller,
   accepted imprecision against a much larger, confirmed stability gain.
2. **The R run is chunked** into ~monthly blocks (fixed row count, *not*
   calendar months parsed from `obs_time` — `obs_time` is UTC-shifted from
   `climdata`'s own local-time-ordered rows, so a UTC-calendar-month split
   creates a tiny leading fragment whose own instability poisons the real
   first month's initial soil state). A cascade is contained to the block
   it starts in instead of losing the rest of the year. Soil state (weeks
   to months of memory) carries forward between blocks via
   `SoilTempIni`/`ThetaIni`, but only from a block that was itself mostly
   clean (>=90% non-NaN) — a heavily-cascaded block's tail state isn't a
   trustworthy initial condition and gets reset instead of propagating the
   corruption indefinitely (confirmed empirically: a clean-in-isolation
   month came out 100% NaN when initialised from the previous month's own
   99.9%-NaN ending state).

Even with both mitigations, roughly half the year still cascades somewhere
(varies by block). **The obs-vs-micropoint stats/plots only cover the hours
that came through valid** — `n` in the printed stats table shows exactly how
much of the year that is per variable (compare against the obs-vs-Julia
`n`, which covers the full year modulo real data gaps). This is a genuine,
only-partially-resolved fragility in micropoint's own solver under real
tower-data noise, not a Microclimate.jl issue — the obs-vs-Julia leg is
unaffected and matches `comparisons/ozflux`'s own independently-run stats
for the same site-year exactly (checked directly, not assumed).

## Other scope notes

- `RunModelFull` returns no sensible/latent heat flux — Fh/Fe comparison
  against micropoint is out of scope (Fsu/Flu/Fn/Ta/Ws/leaf temp/soil are
  in scope).
- micropoint's layered output tops out at `canopy_height` — Whroo's 32m/36m
  sensors (above canopy) stay a 2-way comparison (obs vs Julia), already
  covered by `comparisons/ozflux/report.jl`.
- `canopy_layers.csv` uses micropoint's own required evenly-spaced
  bottom-to-top grid (confirmed via `?RunModelFull`/`?PAIgeometry`), shaped
  with the same `PAI_SHAPES`/`SITE_PAI_SHAPE` density function
  `comparisons/ozflux/pipeline.jl` uses on its own uneven grid, rescaled to
  the site's total LAI — both models start from the same PAI *shape*, not
  just the same total.
- Soil hydraulic params (`psi_e`/`Ksat`/`rho`/`smax`/`b`) come from one
  representative near-surface (~10cm) node of `comparisons/ozflux/config.jl`'s
  `soil_source(site_name)` profile (`:slga` by default, or a
  `CAMPBELL_NORMAN_TEXTURES` key for a fixed literature texture class),
  applied flat across all of micropoint's layers, using the same unit
  conversions documented in `micropoint/config.jl`/`micropoint/README.md`.
  **Not a full per-depth profile** — that was tried (each param interpolated
  onto micropoint's own node depths) and made `RunModelFull` hang
  indefinitely (not a NaN cascade — a genuine non-returning call, confirmed
  via isolated Rscript reproduction and bisected to `Smax` specifically: a
  flat, single, self-consistent `Smax` well below micropoint's own default
  soil-type table (SLGA's real bulk density at Whroo gives ~0.15-0.28 vs.
  micropoint's own default of 0.46) reproduces the hang on its own, depth
  variation or not). `soil_source` set to a texture class instead of `:slga`
  sidesteps this at the root, since e.g. `:sandy_loam`'s default bulk/mineral
  density gives Smax≈0.49, in the same range as micropoint's own defaults.
  `createsoilc`'s texture string is left at `"Clay loam"` (a placeholder) —
  the explicit `smax`/`b`/`psi_e`/`ksat` overrides dominate its effect anyway.
- Organic litter-layer soil override is a per-site toggle
  (`SITE_ORGANIC_CAP`/`organic_cap` in `comparisons/ozflux/config.jl`),
  applied to both the Julia solve and (via `params.csv`'s `organic_cap`
  column) the R side.
- `obs_time`'s UTC offset comes from `SITE_UTC_OFFSET_MINUTES`
  (`comparisons/ozflux/config.jl`), one fixed per-site value (OzFlux records
  in local standard time year-round, no DST jump). Whroo's is the only one
  directly verified against the raw 30-min index around a real DST
  transition date; the rest come from each file's `time_zone` attribute. Not
  the longitude-based solar-time correction
  `micropoint/write_micropoint_inputs.jl` uses — that's a different
  (Microclimate.jl-internal solar-time) problem that doesn't apply to real
  recorded clock time.

## Running other sites/years

Set `site_name`/`years` at the top of `ozflux_below_canopy_comparison.jl` —
any `comparisons/ozflux/data` site works (all have a
`SITE_UTC_OFFSET_MINUTES` entry). Whroo/2015 was the first tried and has the
best real sub-canopy sensor coverage of any site here; other site-years may
cascade more (see above) or need `profile_times` moved to a cleanly-solved
window.
