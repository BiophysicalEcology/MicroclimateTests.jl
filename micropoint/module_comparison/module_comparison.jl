# module_comparison.jl — compares Microclimate.jl's low-level canopy
# sub-model functions against the real R functions the below-canopy
# vignette (github.com/ilyamaclean/microclimlearn) uses: twostream,
# longwavebelow, windprofile_below, SolveEnergyBalance, temp_below/
# relhum_below. Run run_micropoint.R first to produce outputs/.
#
# Run with:
#   cd c:/git/BiophysicalEcologyEnv
#   julia --project=. c:/git/MicroclimateTests.jl/micropoint/module_comparison/module_comparison.jl

using Microclimate, Unitful, Statistics, Printf, CSV, DataFrames, Plots

gr()
outdir = joinpath(@__DIR__, "outputs")
struct Stat; r::Float64; rmse::Float64; bias::Float64; end
function compute_stats(a, b)
    m = isfinite.(a) .& isfinite.(b)
    count(m) < 2 && return Stat(NaN, NaN, NaN)
    Stat(cor(a[m], b[m]), sqrt(mean((b[m] .- a[m]) .^ 2)), mean(b[m] .- a[m]))
end
fmt(s::Stat) = @sprintf("r=%+.4f  RMSE=%8.3f  bias=%+8.3f", s.r, s.rmse, s.bias)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  1: two-stream shortwave                                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝
println("="^72); println("1: two-stream shortwave"); println("="^72)

vegp = DataFrame(CSV.File(joinpath(outdir, "vegp.csv")))[1, :]
groundp = DataFrame(CSV.File(joinpath(outdir, "groundp.csv")))[1, :]
d1 = DataFrame(CSV.File(joinpath(outdir, "01_paii.csv")))
inp1 = DataFrame(CSV.File(joinpath(outdir, "01_inputs.csv")))[1, :]
r1 = DataFrame(CSV.File(joinpath(outdir, "01_twostream.csv")))

n1 = nrow(d1)
canopy_height1 = vegp.h * 1.0u"m"
plant_area_index1 = reverse(d1.paii)  # R bottom-to-top -> Julia top-to-bottom
zenith_angle1 = inp1.zen * 1.0u"°"
direct1 = (inp1.swdown - inp1.difrad) * 1.0u"W/m^2"
diffuse1 = inp1.difrad * 1.0u"W/m^2"

sw_model = TwoStreamRadiation(; leaf_reflectance=vegp.lref, leaf_transmittance=vegp.ltra)
sw_buf = Microclimate.allocate_shortwave(sw_model, canopy_height1, plant_area_index1, n1, vegp.x)
Microclimate.canopy_shortwave!(sw_buf, sw_model, plant_area_index1;
    zenith_angle=zenith_angle1, direct_horizontal_irradiance=direct1, diffuse_horizontal_irradiance=diffuse1,
    ground_reflectance=groundp.gref)

jl_Rdirdown = reverse(ustrip.(u"W/m^2", sw_buf.downward_direct_shortwave))  # -> bottom-to-top, matching R
jl_Rdifdown = reverse(ustrip.(u"W/m^2", sw_buf.downward_diffuse_shortwave))
jl_Rswup    = reverse(ustrip.(u"W/m^2", sw_buf.upward_diffuse_shortwave))
# R's Rleafabs is per unit LEAF area; absorbed_shortwave is per unit GROUND
# area (documented on the buffer) -- convert via the same leaf-area factor
# canopy_energy_balance! itself uses (2*(1-layer_transmission)).
layer_transmission1 = exp.(-plant_area_index1)
jl_Rleafabs = reverse(ustrip.(u"W/m^2", sw_buf.absorbed_shortwave) ./ (2.0 .* (1.0 .- layer_transmission1)))

for (label, jl, r) in [("Rdirdown", jl_Rdirdown, r1.Rdirdown), ("Rdifdown", jl_Rdifdown, r1.Rdifdown),
    ("Rswup", jl_Rswup, r1.Rswup), ("Rleafabs", jl_Rleafabs, r1.Rleafabs)]
    println("  " * rpad(label, 12) * fmt(compute_stats(r, jl)))
end

fig1 = plot(layout=(1, 4), size=(1400, 500), dpi=120, left_margin=6Plots.mm, bottom_margin=6Plots.mm)
for (i, (label, jl, r)) in enumerate([("Rdirdown", jl_Rdirdown, r1.Rdirdown), ("Rdifdown", jl_Rdifdown, r1.Rdifdown),
    ("Rswup", jl_Rswup, r1.Rswup), ("Rleafabs", jl_Rleafabs, r1.Rleafabs)])
    plot!(fig1[i], jl, d1.z; label="Julia", color=:tomato, title=label, xlabel="W/m²", ylabel=i == 1 ? "Height (m)" : "")
    plot!(fig1[i], r, d1.z; label="R", color=:steelblue)
end
display(fig1)
savefig(fig1, joinpath(outdir, "01_twostream_comparison.png"))

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  2: below-canopy longwave -- Rlwdown should match closely (both models   ║
# ║  transmit the sky term the same way); Rlwup/RlwLabs are expected to      ║
# ║  diverge in profile shape, not just noise: R's longwavebelow() computes  ║
# ║  the upward/reflected stream via an all-pairs view-factor matrix         ║
# ║  (inter-layer weights from cumulative-PAI distance), while Microclimate.jl's║
# ║  LayeredLongwaveExchange uses a sequential top-down/bottom-up cascade -- ║
# ║  a real, already-documented structural-fidelity difference, not a bug.   ║
# ╚══════════════════════════════════════════════════════════════════════════╝
println("\n" * "="^72); println("2: below-canopy longwave"); println("="^72)

d2 = DataFrame(CSV.File(joinpath(outdir, "02_longwavebelow.csv")))
n2 = nrow(d2)
canopy_height2 = 10.0u"m"
plant_area_index2 = reverse(d2.paii)
leaf_temperature2 = u"K".(reverse(d2.tleaf) .* u"°C")
ground_temperature2 = u"K"(15.0u"°C")

lw_model = LayeredLongwaveExchange()
lw_buf = Microclimate.allocate_longwave(lw_model, plant_area_index2, n2)
site2 = Site(; latitude=0.0u"°", longitude=0.0u"°", elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
    horizon_angles=fill(0.0u"°", 24), sky_view_fraction=1.0, albedo=0.15, roughness_height=0.01u"m",
    atmospheric_pressure=101325.0u"Pa")
env2 = (; atmospheric_pressure=101325.0u"Pa", reference_humidity=0.7, reference_temperature=288.15u"K",
    cloud_emissivity=1.0, cloud_cover=0.5, shade=0.0, surface_emissivity=0.97, longwave_radiation=300.0u"W/m^2")
Microclimate.canopy_longwave!(lw_buf, lw_model, 0.97;
    leaf_temperature=leaf_temperature2, ground_temperature=ground_temperature2, ground_emissivity=0.97,
    site=site2, environment_instant=env2)

jl_Rlwdown = reverse(ustrip.(u"W/m^2", lw_buf.boundary_downward_longwave[1:n2]))
jl_Rlwup   = reverse(ustrip.(u"W/m^2", lw_buf.boundary_upward_longwave[1:n2]))
# R's RlwLabs = 0.5*vegem*(Rlwdown+Rlwup) (micropoint's longwavemodelCpp), gross not net.
jl_RlwLabs = 0.5 .* 0.97 .* (jl_Rlwdown .+ jl_Rlwup)

for (label, jl, r) in [("Rlwdown", jl_Rlwdown, d2.Rlwdown), ("Rlwup", jl_Rlwup, d2.Rlwup), ("RlwLabs", jl_RlwLabs, d2.RlwLabs)]
    println("  " * rpad(label, 14) * fmt(compute_stats(r, jl)))
end

fig2 = plot(layout=(1, 3), size=(1100, 500), dpi=120, left_margin=6Plots.mm, bottom_margin=6Plots.mm)
for (i, (label, jl, r)) in enumerate([("Rlwdown", jl_Rlwdown, d2.Rlwdown), ("Rlwup", jl_Rlwup, d2.Rlwup), ("RlwLabs", jl_RlwLabs, d2.RlwLabs)])
    plot!(fig2[i], jl, d2.z; label="Julia", color=:tomato, title=label, xlabel="W/m²", ylabel=i == 1 ? "Height (m)" : "")
    plot!(fig2[i], r, d2.z; label="R", color=:steelblue)
end
display(fig2)
savefig(fig2, joinpath(outdir, "02_longwave_comparison.png"))

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  3: below-canopy wind                                                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝
println("\n" * "="^72); println("3: below-canopy wind"); println("="^72)

d3 = DataFrame(CSV.File(joinpath(outdir, "03_windprofile.csv")))
n3 = nrow(d3)
canopy_height3 = 25.0u"m"
karman_constant = MoninObukhov().karman_constant
uh3 = 3.0u"m/s"

pai_skew3 = reverse(d3.paii_skew)
pai_unif3 = reverse(d3.paii_unif)
jl_wind_skew = reverse(Microclimate.wind_attenuation_profile(pai_skew3, canopy_height3, sum(pai_skew3), karman_constant)) .* uh3
jl_wind_unif = reverse(Microclimate.wind_attenuation_profile(pai_unif3, canopy_height3, sum(pai_unif3), karman_constant)) .* uh3

for (label, jl, r) in [("wind (skewed PAI)", ustrip.(u"m/s", jl_wind_skew), d3.wind_skew), ("wind (uniform PAI)", ustrip.(u"m/s", jl_wind_unif), d3.wind_unif)]
    println("  " * rpad(label, 20) * fmt(compute_stats(r, jl)))
end

fig3 = plot(layout=(1, 2), size=(1000, 500), dpi=120, left_margin=6Plots.mm, bottom_margin=6Plots.mm)
plot!(fig3[1], ustrip.(u"m/s", jl_wind_skew), d3.z; label="Julia", color=:tomato, title="Skewed PAI", xlabel="m/s", ylabel="Height (m)")
plot!(fig3[1], d3.wind_skew, d3.z; label="R", color=:steelblue)
plot!(fig3[2], ustrip.(u"m/s", jl_wind_unif), d3.z; label="Julia", color=:tomato, title="Uniform PAI", xlabel="m/s")
plot!(fig3[2], d3.wind_unif, d3.z; label="R", color=:steelblue)
display(fig3)
savefig(fig3, joinpath(outdir, "03_wind_comparison.png"))

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  4: leaf temperature -- fed R's own Rabs/conductance/wind directly, to   ║
# ║  isolate the energy-balance solve itself from stomatal-model differences ║
# ╚══════════════════════════════════════════════════════════════════════════╝
println("\n" * "="^72); println("4: leaf temperature (R's Rabs/gs/wind fed directly into Julia's solver)"); println("="^72)

d4 = DataFrame(CSV.File(joinpath(outdir, "04_leaftemperature.csv")))
n4 = nrow(d4)
vegp4_len, vegp4_wid = 0.15u"m", 0.07u"m"  # createplant_inputs("BDT")
body4 = leaf_body(vegp4_len, vegp4_wid, 1.6)

air_temperature4 = fill(u"K"(25.0u"°C"), n4)
relative_humidity4 = fill(0.80, n4)
atmospheric_pressure4 = 101.3u"kPa" |> p -> uconvert(u"Pa", p)
leaf_water_potential4 = -0.5e6u"J/kg"

jl_tleaf4 = zeros(n4)
for i in 1:n4
    absorbed_radiation_i = d4.Rabs[i] * 1.0u"W/m^2"
    wind_speed_i = d4.uz[i] * 1.0u"m/s"
    stomatal_conductance_i = LeafEvaporationParameters(;
        abaxial_vapour_conductance=d4.gs[i] * 1.0u"mol/m^2/s", adaxial_vapour_conductance=0.0u"mol/m^2/s",
        cuticular_conductance=0.0u"mol/m^2/s")
    solved = Microclimate.leaf_temperature(RootFindLeafTemperature(), absorbed_radiation_i, air_temperature4[i],
        relative_humidity4[i], wind_speed_i, atmospheric_pressure4, 0.97, stomatal_conductance_i,
        leaf_water_potential4, body4; leaf_area=1.0u"m^2", leaf_temperature_guess=air_temperature4[i])
    jl_tleaf4[i] = ustrip(u"°C", solved)
end

println("  " * rpad("tleaf", 12) * fmt(compute_stats(d4.tleaf, jl_tleaf4)))

fig4 = plot(size=(600, 500), dpi=120, xlabel="Leaf temperature (°C)", ylabel="Height (m)", title="Leaf temperature")
plot!(fig4, jl_tleaf4, d4.z; label="Julia", color=:tomato, lw=2)
plot!(fig4, d4.tleaf, d4.z; label="R", color=:steelblue, lw=2)
savefig(fig4, joinpath(outdir, "04_leaftemperature_comparison.png"))
display(fig4)
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  5: below-canopy air temperature and humidity -- full canopy_energy_     ║
# ║  balance! Picard solve vs. R's hand-rolled Aitken-relaxed loop. Different ║
# ║  air-profile theory/iteration on each side (KTheoryAirProfile vs.        ║
# ║  Raupach LNF) and different stomatal models, so qualitative agreement is ║
# ║  the bar here, not the tight match sections 1-4 aim for.                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝
println("\n" * "="^72); println("5: below-canopy air temperature/humidity (qualitative)"); println("="^72)

d5 = DataFrame(CSV.File(joinpath(outdir, "05_airtemphumidity.csv")))
can5 = DataFrame(CSV.File(joinpath(outdir, "05_wholecanopy.csv")))[1, :]
n5 = nrow(d5)
canopy_height5 = 20.0u"m"  # createplant_inputs("BDT")$h
heights5_below = collect(1:n5) .* (20.0 / n5) .* u"m"
heights5 = vcat(heights5_below, [22.0u"m"])
d5_paii = DataFrame(CSV.File(joinpath(outdir, "05_paii.csv")))
plant_area_index5 = reverse(d5_paii.paii)  # R bottom-to-top -> Julia top-to-bottom

boundary_layer_model5 = MoninObukhov()
model5 = MultilayerCanopy(;
    canopy_height=canopy_height5, plant_area_index=plant_area_index5,
    shortwave_model=TwoStreamRadiation(; leaf_reflectance=0.4, leaf_transmittance=0.2),
    leaf_parameters=LeafParameters(; leaf_length=0.15u"m", leaf_width=0.07u"m", leaf_emissivity=0.97, leaf_angle_distribution_parameter=1.6),
    leaf_temperature_solver=RootFindLeafTemperature(),
    air_profile_model=RaupachLTheoryAirProfile(; far_field_mode=Val(:bulk)),  # R's temp_below/relhum_below are Raupach LNF, not K-theory
    convergence_model=PicardCanopyConvergence(;
        convergence=IterationToleranceConvergence(; tolerance=0.01u"K", max_iterations_per_day=30)),
)
buffers5 = Microclimate.allocate_canopy(model5, heights5, boundary_layer_model5)

site5 = Site(; latitude=49.96807u"°", longitude=-5.215668u"°", elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
    horizon_angles=fill(0.0u"°", 24), sky_view_fraction=1.0, albedo=0.15, roughness_height=0.004u"m",
    atmospheric_pressure=101.399u"kPa" |> p -> uconvert(u"Pa", p))
environment_instant5 = (;
    atmospheric_pressure=site5.atmospheric_pressure, reference_temperature=u"K"(25.684u"°C"),
    reference_wind_speed=1.431u"m/s", reference_humidity=0.5372, zenith_angle=30.0u"°",
    cloud_emissivity=1.0, cloud_cover=0.2, surface_emissivity=0.97, shade=0.0,
    longwave_radiation=414.7725u"W/m^2",
)
inputs5 = Microclimate.CanopyEnergyBalanceInputs(model5;
    site=site5, environment_instant=environment_instant5, zenith_angle=30.0u"°",
    direct_horizontal_irradiance=(839.008 - 223.0) * 1.0u"W/m^2", diffuse_horizontal_irradiance=223.0u"W/m^2",
    ground_reflectance=0.15, ground_temperature=u"K"(22.0u"°C"), ground_emissivity=0.97,
    ground_relative_humidity=0.87, canopy_source_temperature=u"K"(25.684u"°C"),
)
result5 = Microclimate.canopy_energy_balance!(buffers5, model5, boundary_layer_model5, inputs5)
println("Julia Picard iterations: ", result5.iterations, "  (R Aitken iterations: see run_micropoint.R output)")

jl_tleaf5 = reverse(ustrip.(u"°C", buffers5.leaf.leaf_temperature))
jl_tair5 = reverse(ustrip.(u"°C", buffers5.air_profile.air_temperature))
jl_rh5 = reverse(buffers5.air_profile.relative_humidity .* 100.0)

for (label, jl, r) in [("tleaf", jl_tleaf5, d5.tleaf), ("tair", jl_tair5, d5.tair), ("rh", jl_rh5, d5.rh)]
    println("  " * rpad(label, 12) * fmt(compute_stats(r, jl)))
end

fig5 = plot(layout=(1, 3), size=(1300, 500), dpi=120, left_margin=6Plots.mm, bottom_margin=6Plots.mm)
plot!(fig5[1], jl_tleaf5, d5.z; label="Julia", color=:tomato, title="Leaf temp", xlabel="°C", ylabel="Height (m)")
plot!(fig5[1], d5.tleaf, d5.z; label="R", color=:steelblue)
plot!(fig5[2], jl_tair5, d5.z; label="Julia", color=:tomato, title="Air temp", xlabel="°C")
plot!(fig5[2], d5.tair, d5.z; label="R", color=:steelblue)
plot!(fig5[3], jl_rh5, d5.z; label="Julia", color=:tomato, title="Relative humidity", xlabel="%")
plot!(fig5[3], d5.rh, d5.z; label="R", color=:steelblue)
savefig(fig5, joinpath(outdir, "05_airtemphumidity_comparison.png"))
display(fig5)

println("\nDone -- comparison plots saved to $outdir")
