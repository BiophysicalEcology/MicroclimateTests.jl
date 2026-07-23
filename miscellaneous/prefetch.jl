include("c:/git/MicroclimateTests.jl/miscellaneous/prefetch_silo.jl")
ENV["RASTERDATASOURCES_PATH"] = "z:/"

# default layers (all 6 points_australia.jl needs): max_temp, min_temp,
# daily_rain, rh_tmax, rh_tmin, radiation
prefetch_and_rechunk_silo(2010:2011)

# or just a subset of layers:
prefetch_and_rechunk_silo(2025:2026; layers = (:daily_rain,))