# write_ozflux_micropoint_inputs.jl — writes climdata.csv/params.csv/
# canopy_params.csv/canopy_layers.csv from a solved comparisons/ozflux
# result, for run_micropoint_ozflux_vegetated.R to read.
#
# Requires comparisons/ozflux/{config,utils,pipeline}.jl already included
# (for site_albedo/SITE_ALBEDO, emissivity, PAI_SHAPES, SITE_PAI_SHAPE,
# SITE_UTC_OFFSET_MINUTES, flatten_soil_profile).

# micropoint needs an evenly-spaced 0..canopy_height layer grid -- this
# evaluates the same density shape function pipeline.jl uses (PAI_SHAPES/
# SITE_PAI_SHAPE) on that even grid instead, so both models start from the
# same PAI *shape*, not just the same total. z_m/paii are bottom(1)-to-
# top(n), micropoint's own convention.
function micropoint_canopy_layers(site_name, canopy_height, target_pai; n=50)
    shape = PAI_SHAPES[get(SITE_PAI_SHAPE, site_name, :uniform)]
    layer_thickness = canopy_height / n
    z_m = [(i / n) * canopy_height for i in 1:n]
    density = shape(z_m, canopy_height)
    raw = density .* layer_thickness
    paii = raw .* (target_pai / sum(raw))
    return ustrip.(u"m", z_m), ustrip.(paii)
end

function write_ozflux_micropoint_inputs(result, outdir)
    mkpath(outdir)
    (; output, t_model, resolved, forcing_used) = result
    n = length(t_model)

    haskey(SITE_UTC_OFFSET_MINUTES, resolved.site_name) ||
        error("No UTC offset configured for site \"$(resolved.site_name)\" -- add it to SITE_UTC_OFFSET_MINUTES in comparisons/ozflux/config.jl")
    utc_offset = Minute(SITE_UTC_OFFSET_MINUTES[resolved.site_name])
    obs_time_str = Dates.format.(t_model .- utc_offset, dateformat"yyyy-mm-dd HH:MM:SS")

    # difrad reuses output.diffuse_fraction (not an independent recomputation
    # from real Fsd) so both models see the same direct/diffuse split
    # Microclimate.jl's own canopy solve actually used internally
    #
    # windspeed is forcing_used.Ws (gap-filled tower+SILO reference wind), not
    # output.canopy.wind_speed[:,1] -- that's Julia's own canopy-resolved
    # wind, which shifts with canopy/convergence settings and would make
    # micropoint's input (and so its output) implicitly track Julia's own
    # parameter choices instead of being an independent comparison baseline.
    # Floor guards the same near-zero-wind crash the old value happened to
    # avoid by construction.
    windspeed = forcing_used.Ws
    climdata = DataFrame(
        obs_time  = obs_time_str,
        temp      = forcing_used.Ta,
        relhum    = forcing_used.RH,
        pres      = ustrip.(u"kPa", output.pressure[1:n]),
        swdown    = forcing_used.Fsd,
        difrad    = forcing_used.Fsd .* output.diffuse_fraction[1:n],
        lwdown    = forcing_used.Fld,
        windspeed = max.(windspeed, 0.1),
        winddir   = zeros(n),
        precip    = forcing_used.Precip,
    )
    CSV.write(joinpath(outdir, "climdata.csv"), climdata)

    # micropoint takes one flat slab -- flatten_soil_profile (comparisons/
    # ozflux/utils.jl) depth-weight-averages whatever profile the Julia run
    # actually used (a no-op if soil_source was already :slga_uniform, a real
    # collapse if it was the full per-depth :slga -- either way this is the
    # SAME flattening Microclimate.jl itself applies for :slga_uniform, so a
    # :slga_uniform run gives both models the exact same numbers). Feeding
    # micropoint the full per-depth SLGA profile directly (not flattened at
    # all) previously made RunModelFull hang indefinitely (bisected to Smax;
    # see README).
    #
    # soil_source==:slga/:slga_uniform: createsoilc(soiltype="Clay loam") is
    # only a structural placeholder (Vq/Vm/Vo/Smin/n), overridden below with
    # the real SLGA-derived Smax/b/psi_e/Ksat/rho. soil_source==texture
    # class: that mix-and-match caused a total NaN cascade (Vq/Vm/Vo from
    # "Clay loam" paired with e.g. sandy loam's Smax/b/psi_e/Ksat) -- so
    # instead `soiltype` names the real texture and `soil_override=false`
    # tells R to trust createsoilc(soiltype=<that texture>)'s own fully
    # self-consistent table verbatim (still flat/non-layered; organic_cap can
    # still apply).
    source = soil_source(resolved.site_name)
    flat = flatten_soil_profile(result.problem.inputs.soil_profile, result.depths)
    smax = 1.0 - ustrip(u"Mg/m^3", flat.bulk_density[1]) / ustrip(u"Mg/m^3", flat.mineral_density[1])
    b = flat.hydraulics.campbell_b_parameter[1]
    psi_e = abs(ustrip(u"J/kg", flat.hydraulics.air_entry_water_potential[1])) / ustrip(u"m/s^2", Unitful.gn)
    ksat = ustrip(flat.hydraulics.saturated_hydraulic_conductivity[1])
    rho_kgm3 = ustrip(u"kg/m^3", flat.bulk_density[1])
    is_slga = source === :slga || source === :slga_uniform
    soiltype = is_slga ? "Clay loam" : CAMPBELL_NORMAN_MICROPOINT_NAME[source]
    soil_override = is_slga

    canopy_model = result.problem.model.canopy_model
    (; leaf_parameters, shortwave_model, interception_model) = canopy_model
    canopy_height_m = ustrip(u"m", canopy_model.canopy_height)
    pai = sum(canopy_model.plant_area_index)

    params = DataFrame(
        lat = [ustrip(u"°", resolved.latitude)], lon = [ustrip(u"°", resolved.longitude)],
        elev_m = [ustrip(u"m", result.problem.inputs.site.elevation)],
        zref = [ustrip(u"m", resolved.reference_height)],
        gref = [site_albedo(resolved.site_name)], groundem = [emissivity],
        nlayers = [15], totaldepth = [2.0],
        organic_cap = [organic_cap(resolved.site_name)],
        soiltype = [soiltype], soil_override = [soil_override],
        smax = [smax], b = [b], psi_e = [psi_e], ksat = [ksat], rho_kgm3 = [rho_kgm3],
    )
    CSV.write(joinpath(outdir, "params.csv"), params)
    canopy_params = DataFrame(
        canopy_height_m = [canopy_height_m], pai = [pai],
        leaf_length_m = [ustrip(u"m", leaf_parameters.leaf_length)],
        leaf_width_m = [ustrip(u"m", leaf_parameters.leaf_width)],
        leaf_emissivity = [leaf_parameters.leaf_emissivity],
        leaf_reflectance = [shortwave_model.leaf_reflectance],
        leaf_transmittance = [shortwave_model.leaf_transmittance],
        leaf_angle_x = [leaf_parameters.leaf_angle_distribution_parameter],
        leaf_water_storage_mm = [ustrip(u"kg/m^2", interception_model.leaf_water_storage_capacity)],
    )
    CSV.write(joinpath(outdir, "canopy_params.csv"), canopy_params)

    z_m, paii = micropoint_canopy_layers(resolved.site_name, canopy_model.canopy_height, pai)
    CSV.write(joinpath(outdir, "canopy_layers.csv"), DataFrame(z_m = z_m, paii = paii))

    println("Wrote climdata.csv ($(nrow(climdata)) rows), params.csv, canopy_params.csv, canopy_layers.csv to $outdir")
    return outdir
end
