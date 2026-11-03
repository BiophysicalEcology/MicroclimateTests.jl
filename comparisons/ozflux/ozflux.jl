# ozflux.jl — the OzFlux comparison pipeline, in dependency order. Entry points include this,
# with MicroclimateMapper, Rasters, ArchGDAL, RasterDataSources, PointDataSources, GeoInterface
# and Serialization loaded for mapper.jl.
include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "tables.jl"))
include(joinpath(@__DIR__, "data_corrections.jl"))
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "sites.jl"))
include(joinpath(@__DIR__, "barra.jl"))
include(joinpath(@__DIR__, "forcing_selection.jl"))
include(joinpath(@__DIR__, "pipeline.jl"))
include(joinpath(@__DIR__, "report.jl"))
include(joinpath(@__DIR__, "mapper.jl"))
