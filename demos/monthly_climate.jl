using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile,
    example_soil_properties_model, example_soil_hydraulic_model
using Rasters, RasterDataSources
using Rasters.Extents: Extent
using Dates, Statistics, Unitful, Plots

ENV["RASTERDATASOURCES_PATH"] = "z:\\"#"c:/Spatial_Data/" #"/data/scratch/projects/punim0593/rasters" #

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"

points = [geocode("Alice Springs, Australia"),
          geocode("Sydney, Australia"),
          geocode("Melbourne, Australia"),
          geocode("Perth, Australia"),
          geocode("Adelaide, Australia"),
          geocode("Brisbane, Australia"),
          geocode("Darwin, Australia"),
          geocode("Hobart, Australia"),
          geocode("Canberra, Australia"),
          geocode("Cairns, Australia"),
          geocode("Broome, Australia"),
          geocode("Katherine, Australia"),]

# CRUCL2 is a 1961–1990 climatology — the year is ignored; any year works.
dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
    ),
    dem_source              = CRUCL2,
    weather_source          = CRUCL2, #WorldClim{Climate},
    #soil_moisture_source    = CPCSoil,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = false,
)
problem = MicroVectorProblem(;
    model,
    points,
    dates,
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths))),
)

@time output = solve(problem);

site = 1
soil_T = output.soil_temperature[point=site]           # (Ti × depth)
air_T  = output.air_temperature[point=site, height=2]  # 1.2 m
rad    = output.global_radiation[point=site]
rh     = output.relative_humidity[point=site, height=2]
ws     = output.wind_speed[point=site, height=2]

ti           = lookup(soil_T, Ti)
depth_labels = reshape(string.(depths), 1, :)

p1 = plot(collect(uconvert.(u"°C", soil_T)); label = depth_labels,
          title = points[site].display_name, legend = false)
p2 = plot(collect(uconvert.(u"°C", air_T));  title = "Air temperature (1.2 m)",   legend = false)
p3 = plot(collect(rad);    title = "Global radiation",           legend = false)
p4 = plot(collect(rh);     title = "Relative humidity (1.2 m)", legend = false)
p5 = plot(collect(ws);     title = "Wind speed (1.2 m)",        legend = false)

display(plot(p1, p2, p3, p4, p5; layout = (5, 1), size = (900, 1100), link = :x))

# =============================================================================
# Raster run — central Australia at CRUCL2 (10-minute) resolution
# =============================================================================

raster_area = Extent(X = (110.0, 154), Y = (-45.0, -10))
raster_area = Extent(X = (144.0, 149), Y = (-44.0, -40))
location = "Alice Springs, Australia"
raster_area = geocode(location; buffer = 2.0).extent

# Use the CRUCL2 elv layer as the template so the output is at 10-minute resolution.
crucl2_template = load_template(CRUCL2, raster_area)

display(plot(crucl2_template;
    seriestype = :heatmap,
    title = "CRUCL2 elevation — raster area",
    xlabel = "Longitude", ylabel = "Latitude",
    yflip = false))

raster_problem = MicroRasterProblem(;
    model,
    area         = raster_area,
    dates,
    template     = crucl2_template,
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
)

@time raster_output = solve(raster_problem)

# Annual max/min soil surface temperature (depth index 1 = 0 cm).
T_surface = raster_output.soil_temperature[depth=1]
max_Tsoil = uconvert.(u"°C", dropdims(maximum(T_surface; dims = Ti); dims = Ti))
min_Tsoil = uconvert.(u"°C", dropdims(minimum(T_surface; dims = Ti); dims = Ti))

T_palette = cgrad([:blue, :lightblue, :orange, :red, :purple])

display(plot(max_Tsoil;
    title = "Annual max soil surface temperature (CRUCL2)",
    seriestype = :heatmap, color = T_palette, yflip = false))

display(plot(min_Tsoil;
    title = "Annual min soil surface temperature (CRUCL2)",
    seriestype = :heatmap, color = T_palette, yflip = false))

display(plot(max_Tsoil .- min_Tsoil;
    title = "Annual soil surface temperature range (CRUCL2)",
    seriestype = :heatmap, color = T_palette, yflip = false))

# now SRTM DEM as the template

model_srtm = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
    ),
    dem_source              = SRTM,
    weather_source          = CRUCL2,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = true,
)

#srtm_area = Extent(X = (132.8, 133.0), Y = (-22.2, -22.0))
#srtm_area = Extent(X = (132.673, 132.773), Y = (-23.684, -23.584))
location = "Alice Springs, Australia"
srtm_area = geocode(location; buffer = 0.025).extent

srtm_template = load_template(SRTM, srtm_area)

display(plot(srtm_template;
    seriestype = :heatmap,
    title = "SRTM elevation",
    xlabel = "Longitude", ylabel = "Latitude",
    color = cgrad([:green, :tan, :brown, :white]),
    yflip = false))

srtm_problem = MicroRasterProblem(;
    model    = model_srtm,
    area     = srtm_area,
    dates,
    template = SRTM,          # framework loads + crops the tile automatically
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
)

@time srtm_output = solve(srtm_problem)

# Annual max/min soil surface temperature (depth index 1 = 0 cm).
# Strip units before min/max to avoid Julia's typemax(Float64)=Inf initialiser
# being rejected by Unitful's unit-compatibility check; reattach after.
T_surface = srtm_output.soil_temperature[depth=1]
max_Tsoil = uconvert.(u"°C", dropdims(maximum(T_surface; dims = Ti); dims = Ti))
min_Tsoil = uconvert.(u"°C", dropdims(minimum(T_surface; dims = Ti); dims = Ti))

T_palette = cgrad([:blue, :lightblue, :orange, :red, :purple])

display(plot(max_Tsoil;
    title = "Annual max soil surface temperature — Alice Springs SRTM",
    seriestype = :heatmap, color = T_palette, yflip = false))

display(plot(min_Tsoil;
    title = "Annual min soil surface temperature — Alice Springs SRTM",
    seriestype = :heatmap, color = T_palette, yflip = false))

display(plot(max_Tsoil .- min_Tsoil;
    title = "Annual soil surface temperature range — Alice Springs SRTM",
    seriestype = :heatmap, color = T_palette, yflip = false))

# =============================================================================
# Solar-only mode — clear-sky spectral output with SRTM terrain correction
# =============================================================================
# solar_only = true skips the microclimate ODE entirely and runs only the
# SolarRadiation.jl clear-sky spectral model (slope/aspect/horizon-corrected).
# CRUCL2 still provides elevation and atmospheric parameters for the solar calc.

#srtm_area = Extent(X = (132.673, 132.773), Y = (-23.684, -23.584))
#srtm_area = Extent(X = (133.774, 133.986), Y = (-23.807, -23.592))
location = "Alice Springs, Australia"
srtm_area = geocode(location; buffer = 0.1).extent

srtm_template = read(crop(Raster(SRTM; extent = srtm_area, lazy = true); to = srtm_area, touches = true))

display(plot(srtm_template;
    seriestype = :heatmap,
    title = "SRTM elevation",
    xlabel = "Longitude", ylabel = "Latitude",
    color = cgrad([:green, :tan, :brown, :white]),
    yflip = false))

model_solar_clearsky = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
    ),
    dem_source              = SRTM,
    weather_source          = CRUCL2,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = true,
    solar_only              = true,
    solar_output_layers     = (
        SOLAR_BROADBAND,    # :global_radiation  — broadband on horizontal surface
        SOLAR_PAR,          # :par               — 400–700 nm integral
        SOLAR_UVB,          # :uv_b              — 290–315 nm integral
        SOLAR_NIR,          # :nir               — 700–4000 nm integral
        SolarOutputLayer(; name = :global_terrain,     component = :global_terrain),
        SolarOutputLayer(; name = :direct_horizontal,  component = :direct_horizontal),
        SolarOutputLayer(; name = :diffuse_horizontal, component = :diffuse_horizontal),
    ),
)

solar_clearsky_problem = MicroRasterProblem(;
    model        = model_solar_clearsky,
    area         = srtm_area,
    dates,
    template     = SRTM,
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
)

# Use init/solve! explicitly so we can extract the terrain for reuse.
# The horizon-angle sweep over 14 000+ SRTM pixels is the expensive step;
# reusing it saves that cost on every subsequent run over the same area.
solar_clearsky_cache = init(solar_clearsky_problem)
srtm_terrain = terrain(solar_clearsky_cache)   # RasterStack — save for reuse
@time solar_clearsky = solve!(solar_clearsky_cache)

rad_palette = cgrad([:black, :blue, :lightblue, :orange, :yellow, :gold, :purple]);

# options for the solar output layers in the MicroRasterProblem:
# diffuse_horizontal
# direct_horizontal
# global_radiation
# global_terrain
# nir
# par
# uv_b
variable_to_plot = dropdims(maximum(solar_clearsky.global_terrain; dims = Ti); dims = Ti)

display(plot(variable_to_plot;
    title = "",
    seriestype = :heatmap, color = rad_palette, yflip = false))

# =============================================================================
# Cloud-corrected solar-only — Ångström-scaled clear-sky, no microclimate ODE
# =============================================================================
# solar_only=true + cloud_correct_solar=true: runs SolarRadiation.jl clear-sky
# model, then multiplies each hourly output by the Ångström sunshine fraction
# sf = a + b × (1 − cloud)^γ  (default a=0.36, b=0.64, γ=1).
# Weather is loaded for cloud cover; the microclimate ODE is skipped entirely.
# Terrain is reused from the clear-sky cache — no DEM re-download or horizon sweep.

model_solar_cloud = MicroMapModel(;
    micro_model = MicroModel(;
        depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = NoSnow(),
    ),
    dem_source              = SRTM,
    weather_source          = CRUCL2,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain         = true,
    solar_only              = true,
    cloud_correct_solar     = true,
    solar_output_layers     = (
        SolarOutputLayer(; name = :global_terrain,     component = :global_terrain),
        SolarOutputLayer(; name = :global_horizontal,  component = :global_horizontal),
        SolarOutputLayer(; name = :direct_horizontal,  component = :direct_horizontal),
        SolarOutputLayer(; name = :diffuse_horizontal, component = :diffuse_horizontal),
        SOLAR_PAR,
    ),
)

solar_cloud_problem = MicroRasterProblem(;
    model        = model_solar_cloud,
    area         = srtm_area,
    dates,
    template     = SRTM,
    soil_profile = example_soil_profile(depths),
    init         = (; soil_moisture = fill(0.2, length(depths))),
    data         = (; terrain = srtm_terrain),  # reuse terrain — no DEM re-download or horizon sweep
)

@time solar_cloud_output = solve(solar_cloud_problem)

# =============================================================================
# Solar output plots — snapshot, monthly daily total, annual total
# =============================================================================
# Ti layout for CRUCL2 monthly-resolution data:
#   12 months × 24 hours = 288 steps
#   month m (1-based), hour h (0-based) → Ti index (m-1)*24 + h + 1
#
# All energy plots use MJ/m²: W/m² × 1 h × 3600 s/h ÷ 1e6 J/MJ = × 0.0036

ti_of(month, hour) = (month - 1) * 24 + hour + 1

# Compare cloud-adjusted vs clear-sky: ratio is the effective monthly sunshine
# fraction as seen by the terrain-corrected surface.
sf_terrain = solar_cloud_output.global_terrain ./ solar_clearsky.global_terrain

display(plot(sf_terrain[Ti(ti_of(7, 12))];
    title  = "Cloud-adjusted / clear-sky ratio — July noon terrain radiation",
    seriestype = :heatmap, color = cgrad([:white, :blue]), clims = (0, 1),
    yflip = false))

# Variable to plot — swap to any layer in solar_clearsky, e.g. :par, :global_terrain
solar_var = solar_clearsky.global_terrain

# --- 1. Snapshot: irradiance at a chosen hour of a chosen month (W/m²) -------
snap_month = 7    # July
snap_hour  = 12   # local solar hour (0-based; 12 = noon)

snap = solar_var[Ti(ti_of(snap_month, snap_hour))]

month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

display(plot(snap;
    title  = "Clear-sky terrain radiation at $(snap_hour):00 — $(month_names[snap_month])",
    seriestype = :heatmap, color = rad_palette, yflip = false))

# --- 2. Daily total for each month, plot one chosen month --------------------
# Each hourly step represents 1 h; summing and multiplying by 1 hr gives W·hr/m²,
# which uconvert then expresses as MJ/m².

function solar_monthly_daily_totals(raster)
    [dropdims(
        uconvert.(u"MJ/m^2",
            sum(raster[Ti((m-1)*24+1 : m*24)]; dims = Ti) .* 1u"hr");
        dims = Ti)
     for m in 1:12]
end

monthly_daily = solar_monthly_daily_totals(solar_var)

plot_month = 7   # change to 1–12 to plot a different month

display(plot(monthly_daily[plot_month];
    title  = "Daily total clear-sky terrain radiation — $(month_names[plot_month])",
    seriestype = :heatmap, color = rad_palette, yflip = false))

# Convenience: show all 12 months in a small-multiples layout with shared z scale
clim_lo = minimum(minimum(m) for m in monthly_daily)
clim_hi = maximum(maximum(m) for m in monthly_daily)
p_months = [plot(monthly_daily[m];
                 title = month_names[m], seriestype = :heatmap,
                 color = rad_palette, yflip = false,
                 clims = (clim_lo, clim_hi),
                 colorbar = false, axis = false, ticks = false)
            for m in 1:12]
display(plot(p_months...; layout = (3, 4), size = (1200, 700),
             plot_title = "Monthly daily total clear-sky terrain radiation"))

# --- 3. Annual total: daily total × days in each month ----------------------
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

annual_solar = sum(monthly_daily[m] .* days_per_month[m] for m in 1:12)

display(plot(annual_solar;
    title  = "Annual total clear-sky terrain radiation",
    seriestype = :heatmap, color = rad_palette, yflip = false))

# --- 4. Animation: radiation over a day for a chosen month -------------------
anim_month = 1   # change to 1–12

clim_lo_anim = minimum(minimum(solar_var[Ti(ti_of(anim_month, h))]) for h in 0:23)
clim_hi_anim = maximum(maximum(solar_var[Ti(ti_of(anim_month, h))]) for h in 0:23)

anim = @animate for h in 0:23
    plot(solar_var[Ti(ti_of(anim_month, h))];
         title      = "$(month_names[anim_month]) — $(lpad(h, 2, '0')):00 local solar time",
         seriestype = :heatmap,
         color      = rad_palette,
         clims      = (clim_lo_anim, clim_hi_anim),
         yflip      = false)
end

gif(anim, "solar_$(month_names[anim_month]).gif"; fps = 1)
    