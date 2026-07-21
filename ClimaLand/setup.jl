# setup.jl — one-time step: install ClimaLand into your *global* Julia
# environment (not a project-specific one).
#
# Julia stacks the global environment (`~/.julia/environments/v1.12`)
# underneath whatever `--project=` is active, so ClimaLand becomes `using`-able
# from scripts run with `--project=c:/git/BiophysicalEcologyEnv` (the shared
# environment every other script in this repo already uses) without touching
# that environment's own Project.toml/Manifest.toml at all.
#
# Run this once from a plain `julia` REPL (no --project flag, so the default
# global environment is active):
#
#   julia
#   julia> include("setup.jl")
#
# then run the comparison the normal way for this repo:
#
#   cd c:/git/BiophysicalEcologyEnv && julia --project=. c:/git/MicroclimateTests.jl/ClimaLand/comparison.jl
#
# If `using ClimaLand` fails when comparison.jl runs with BiophysicalEcologyEnv
# active (version conflict between ClimaLand's dependency tree and something
# already pinned there — e.g. StaticArrays, NCDatasets, DataInterpolations),
# fall back to an isolated environment instead: `Pkg.activate(@__DIR__)` in
# this folder, `Pkg.add("ClimaLand")` there, and run comparison.jl with
# `--project=c:/git/MicroclimateTests.jl/ClimaLand` instead.

using Pkg

if !occursin("environments", Base.active_project())
    @warn "Active project is $(Base.active_project()) — this doesn't look like " *
          "the default global environment. Run plain `julia` (no --project flag) " *
          "before including this script, or ClimaLand will end up installed into " *
          "whatever project *is* active instead of globally."
end

Pkg.add("ClimaLand")

println("\nClimaLand installed globally. Run the comparison with:")
println("  cd c:/git/BiophysicalEcologyEnv && julia --project=. c:/git/MicroclimateTests.jl/ClimaLand/comparison.jl")
