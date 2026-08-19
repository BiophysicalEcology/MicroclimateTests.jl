# compare_kimba.jl
#
# Dew/frost model validation, phase 1: feed NicheMapR's own hourly forcing
# (TAREF/RH/VREF/SOLR/ZEN/TSKYC from micro_global(loc=c(136.4,-33.1), cap=0),
# see generate_reference.R) directly into Microclimate.jl via HourlyTimeseries
# -- bypassing both models' own diel synthesis and any climate-normals-source
# difference (CRUCL2 vs whatever micro_global's built-in database uses)
# entirely. test/micro_testrun_monthly.jl already shows the engine reproduces
# NicheMapR closely given identical MIN/MAX-derived forcing; this isolates
# whether it still holds given the raw hourly values NicheMapR itself used,
# and whether GarrattSegalCondensation reproduces DEW/FROST from osub.f.
#
# If this matches well, the CRUCL2-driven monthly demo's sparser dew than
# micro_global's own Kimba run is a forcing/climate-normals difference, not a
# physics or dew-model bug.

using Microclimate
using FluidProperties
using Unitful
using CSV, DataFrames
using Plots
using Statistics

datadir = joinpath(@__DIR__, "data")

metout = DataFrame(CSV.File(joinpath(datadir, "metout.csv")))
soil_nmr = DataFrame(CSV.File(joinpath(datadir, "soil.csv")))
soilprops = DataFrame(CSV.File(joinpath(datadir, "soilprops.csv")))
site_params = DataFrame(CSV.File(joinpath(datadir, "site_params.csv")))[1, :]
dep = DataFrame(CSV.File(joinpath(datadir, "DEP.csv")))

days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
n_hours = length(days) * 24
@assert nrow(metout) == n_hours

# micro_global's own 10-node output grid, not Julia's finer DEFAULT_DEPTHS --
# so soil_temperature's columns line up 1:1 with soil_nmr's for a full
# depth-resolved comparison, not just the surface node.
depths = dep.DEP .* 1.0u"cm"
n_depths = length(depths)
heights = [site_params.Usrhyt, site_params.Refhyt]u"m"

site = Site(;
    latitude = site_params.latitude * 1.0u"°",
    longitude = site_params.longitude * 1.0u"°",
    elevation = site_params.elev * 1.0u"m",
    slope = site_params.slope * 1.0u"°",
    aspect = site_params.aspect * 1.0u"°",
    horizon_angles = fill(0.0u"°", 24),
    sky_view_fraction = 1.0,
    albedo = site_params.REFL,
    roughness_height = site_params.RUF * 1.0u"m",
    atmospheric_pressure = atmospheric_pressure(site_params.elev * 1.0u"m"),
)
boundary_layer_model = MoninObukhov(; karman_constant = 0.4, dyer_constant = 16.0)

# soiltype=4 lookup in micro_global is homogeneous with depth (all 19 of its
# own nodes identical) -- one constant per property regardless of our depths.
soil_profile = SoilProfile(;
    bulk_density = fill(soilprops.BD[1] * 1.0u"Mg/m^3", n_depths),
    mineral_density = fill(soilprops.DD[1] * 1.0u"Mg/m^3", n_depths),
    mineral_conductivity = fill(2.5u"W/m/K", n_depths),   # micro_global's Thcond default (not soiltype-adjusted)
    mineral_heat_capacity = fill(870.0u"J/kg/K", n_depths), # micro_global's SpecHeat default (not soiltype-adjusted)
    hydraulics = CampbellHydraulicProfile(;
        air_entry_water_potential = fill(soilprops.PE[1] * 1.0u"J/kg", n_depths),
        saturated_hydraulic_conductivity = fill(soilprops.KS[1] * 1.0u"kg*s/m^3", n_depths),
        campbell_b_parameter = fill(soilprops.BB[1], n_depths),
        root_density = fill(0.0, n_depths)u"m/m^3",
    ),
)
soil_properties_model = CampbelldeVriesSoilProperties(;
    de_vries_shape_factor = 0.1, recirculation_power = 4.0, return_flow_threshold = 0.162,
)
soil_hydraulic_model = example_soil_hydraulic_model()

# TSKYC -> downwelling longwave directly, so Julia sees NicheMapR's own sky
# radiation exactly rather than re-deriving it from a cloud-cover fraction.
longwave_radiation = FluidProperties.σ .* u"K".(metout.TSKYC .* u"°C") .^ 4

environment_hourly = HourlyTimeseries(;
    pressure = fill(site.atmospheric_pressure, n_hours),
    reference_temperature = metout.TAREF .* u"°C",
    reference_humidity = metout.RH ./ 100.0,
    reference_wind_speed = metout.VREF .* u"m/s",
    global_radiation = metout.SOLR .* u"W/m^2",
    cloud_cover = zeros(n_hours),  # unused for global_radiation (supplied directly) or sky radiation (longwave_radiation takes priority); only affects the direct/diffuse split
    rainfall = fill(0.0u"kg/m^2", n_hours),  # runmoist=0 in the reference run too -- soil moisture held static either side
    zenith_angle = metout.ZEN .* u"°",
    longwave_radiation = longwave_radiation,
)

# PCTWET (surface wetness driving the G&S energy term) and soil moisture are
# both static in the reference run (runmoist=0): PCTWET==0 and soilmoist==0.5
# at every hour -- confirmed from the raw CSVs.
environment_daily = DailyTimeseries(;
    shade = fill(0.0, length(days)),
    soil_wetness = fill(0.0, length(days)),
    surface_emissivity = fill(site_params.SLE, length(days)),
    cloud_emissivity = fill(1.0, length(days)),  # unused: longwave_radiation takes priority
    rainfall = fill(0.0u"kg/m^2", length(days)),
    deep_soil_temperature = fill(mean(soil_nmr[:, "D200cm"]) * u"°C", length(days)),  # R's own 200cm node is exactly constant (its own tannul boundary condition)
    leaf_area_index = fill(0.0, length(days)),
)

config = MicroConfig(;
    convergence = FixedIterationConvergence(3),  # micro_global's own ndmax=3 default
    soil_moisture_strategy = PrescribedSoilMoisture(; precomputed_soil_moisture = fill(0.5, n_depths, length(days))),
)

model = MicroModel(;
    hours = collect(0.0:1:23.0),
    depths, heights,
    soil_properties_model, soil_hydraulic_model, boundary_layer_model,
    condensation_model = GarrattSegalCondensation(),
    config,
)
inputs = MicroInputs(;
    site, soil_profile,
    environment_minmax = nothing,
    environment_daily, environment_hourly,
    initial_soil_temperature = nothing,
    initial_soil_moisture = fill(0.5, n_depths),
)
problem = MicroProblem(model, inputs; days, time_mode = NonConsecutiveDayMode(; iterations_per_day = 3))

@time output = Microclimate.solve(problem)

# --- Phase 1: surface temperature and near-surface (1 cm) relative humidity ---
soil_surface_julia = uconvert.(u"°C", output.soil_temperature[:, 1])
soil_surface_nmr = soil_nmr[:, "D0cm"] .* u"°C"
rh_1cm_julia = output.profile.relative_humidity[:, 1] .* 100.0
rh_1cm_nmr = metout.RHLOC

t = 1:n_hours
p1 = plot(t, ustrip.(u"°C", soil_surface_julia); label = "Julia", color = :red, title = "Surface (0cm) soil temperature", ylabel = "°C")
plot!(p1, t, ustrip.(u"°C", soil_surface_nmr); label = "NicheMapR", color = :black)
p2 = plot(t, rh_1cm_julia; label = "Julia", color = :red, title = "1cm relative humidity", ylabel = "%")
plot!(p2, t, rh_1cm_nmr; label = "NicheMapR", color = :black)
display(plot(p1, p2; layout = (2, 1), size = (900, 600)))

st_diff = abs.(ustrip.(u"°C", soil_surface_julia) .- ustrip.(u"°C", soil_surface_nmr))
println("Surface temp max abs diff (°C): ", maximum(st_diff))
println("1cm RH max abs diff (%): ", maximum(abs.(rh_1cm_julia .- rh_1cm_nmr)))
imax = argmax(st_diff)
println("Max diff at index $imax (DOY=$(metout.DOY[imax]), TIME=$(metout.TIME[imax])):")
println("  Julia surface T = ", ustrip(u"°C", soil_surface_julia[imax]), "  NMR surface T = ", ustrip(u"°C", soil_surface_nmr[imax]))
println("  TAREF=", metout.TAREF[imax], " RH=", metout.RH[imax], " VREF=", metout.VREF[imax], " SOLR=", metout.SOLR[imax], " ZEN=", metout.ZEN[imax], " TSKYC=", metout.TSKYC[imax])
println("  NMR TALOC=", metout.TALOC[imax], " RHLOC=", metout.RHLOC[imax])
println("Mean surface T diff (°C): ", mean(st_diff), "  RMSE: ", sqrt(mean(st_diff.^2)))
# per-hour trace around the worst index
lo, hi = max(1, imax-3), min(n_hours, imax+3)
for i in lo:hi
    println("  i=$i DOY=$(metout.DOY[i]) TIME=$(metout.TIME[i]) Julia=$(round(ustrip(u"°C",soil_surface_julia[i]);digits=2)) NMR=$(round(ustrip(u"°C",soil_surface_nmr[i]);digits=2)) TAREF=$(metout.TAREF[i]) SOLR=$(metout.SOLR[i])")
end

# --- Phase 2: dew/frost ---
dew_julia = ustrip.(u"kg/m^2", output.ground_dew)
frost_julia = ustrip.(u"kg/m^2", output.ground_frost)
dew_nmr = metout.DEW
frost_nmr = metout.FROST

p3 = plot(t, dew_julia; label = "Julia", color = :red, title = "Dew formed (mm/h)", ylabel = "mm")
plot!(p3, t, dew_nmr; label = "NicheMapR", color = :black)
p4 = plot(t, frost_julia; label = "Julia", color = :red, title = "Frost formed (mm/h)", ylabel = "mm")
plot!(p4, t, frost_nmr; label = "NicheMapR", color = :black)
display(plot(p3, p4; layout = (2, 1), size = (900, 600)))

println("Julia dew events (>0): ", count(>(0), dew_julia), "  max: ", maximum(dew_julia))
println("NicheMapR dew events (>0): ", count(>(0), dew_nmr), "  max: ", maximum(dew_nmr))
println("Julia frost events (>0): ", count(>(0), frost_julia), "  max: ", maximum(frost_julia))
println("NicheMapR frost events (>0): ", count(>(0), frost_nmr), "  max: ", maximum(frost_nmr))

# Nighttime-only (SOLR==0) agreement -- the regime dew/frost actually forms in.
night = metout.SOLR .== 0.0
st_diff_night = abs.(ustrip.(u"°C", soil_surface_julia[night]) .- ustrip.(u"°C", soil_surface_nmr[night]))
rh_diff_night = abs.(rh_1cm_julia[night] .- rh_1cm_nmr[night])
println("\n--- Nighttime-only (SOLR==0, $(count(night)) hours) ---")
println("Surface temp: mean abs diff = ", mean(st_diff_night), "  max = ", maximum(st_diff_night))
println("1cm RH: mean abs diff = ", mean(rh_diff_night), "  max = ", maximum(rh_diff_night))
inight_max = findall(night)[argmax(st_diff_night)]
println("Worst nighttime index $inight_max (DOY=$(metout.DOY[inight_max]), TIME=$(metout.TIME[inight_max])): Julia=$(round(ustrip(u"°C",soil_surface_julia[inight_max]);digits=2)) NMR=$(round(ustrip(u"°C",soil_surface_nmr[inight_max]);digits=2)) TAREF=$(metout.TAREF[inight_max]) RH=$(metout.RH[inight_max]) VREF=$(metout.VREF[inight_max]) TSKYC=$(metout.TSKYC[inight_max])")

# --- Full depth-resolved comparison (now that grids match 1:1) ---
depth_cols = ["D0cm", "D2.5cm", "D5cm", "D10cm", "D15cm", "D20cm", "D30cm", "D50cm", "D100cm", "D200cm"]
println("\n--- Depth-resolved RMSE (°C) ---")
for (i, col) in enumerate(depth_cols)
    julia_col = ustrip.(u"°C", uconvert.(u"°C", output.soil_temperature[:, i]))
    nmr_col = soil_nmr[:, col]
    println("  $(col): RMSE=$(round(sqrt(mean((julia_col .- nmr_col).^2));digits=2))  max_abs_diff=$(round(maximum(abs.(julia_col .- nmr_col));digits=2))")
end
