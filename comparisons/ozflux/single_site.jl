# Entry point for debugging one site: no try/catch, so errors show their full stack trace

using Microclimate, MicroclimateMapper, Unitful, FluidProperties, SolarRadiation, OrdinaryDiffEqSDIRK
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources
using Serialization, Measures

include("ozflux.jl")

# Any site in data/
site_name = "LongreachMitchellGrassRangeland"
years = [2020]
# Days to run, nothing for all
max_days = nothing

result = run_site(site_name, years; max_days)

plot_start = DateTime(years[1], 1, 1, 0)
plot_end = DateTime(years[1], 12, 30, 23)
plot_forcing = false
# nothing plots every variable, or one key such as "Fn", "Ts" or "Ta_HMP_2m"
plot_variable = nothing

stats = report_site_results(result, result.output; plot_start, plot_end, display_plots=true,
    plot_forcing, plot_variable)

# Vertical canopy profiles, obs against model
profile_times = DateTime(years[1], 11, 1, 1):Hour(1):DateTime(years[1], 11, 1, 16)
plot_canopy_profiles_all(result, result.output, profile_times)
