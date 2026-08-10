# MicroclimateTests.jl

Demo, benchmark, and validation scripts for
[MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) —
not a registered package, just something to clone and run. Uses the shared
`BiophysicalEcologyEnv` Julia environment (see `.vscode/settings.json`)
rather than its own `Project.toml`.

## Layout

- `demos/` — get-started runnable examples of MicroclimateMapper.jl itself (point and raster mode).
- `datasources/` — scripts exercising the underlying data-source packages (RasterDataSources.jl, PointDataSources.jl) that MicroclimateMapper.jl fetches forcing from: BARRA, SILO, CRU CL, soil texture, etc.
- `benchmarks/` — like-for-like timing/memory comparisons between MicroclimateMapper.jl and other microclimate models, from reading in data through to a finished simulation.
- `comparisons/` — validation against real station/field observations, with other models (NicheMapR, microclimf, CliMAland.jl) run on identical inputs as a secondary comparison. Each subfolder is a self-contained study (its own `config.jl`/`pipeline.jl`/`comparison.jl`/`single_site.jl`) — see `comparisons/scan_snotel/README.md` for the first one.

Nothing fetched or generated (weather caches, model outputs, station
observations) is committed — see `.gitignore`. Where a comparison needs real
data, its own README says where to get it.
