# run_climaland.jl — build and solve ClimaLand's Soil.EnergyHydrology model,
# driven by the forcing from forcing_climaland.jl, with soil hydraulics/thermal
# parameters translated from the same Campbell/de Vries calibration
# comparisons/scan_snotel/config.jl uses for Microclimate.jl (see config.jl's
# comment block for the translation).
#

function run_climaland(cl_forcing, prep)
    (; t_sec, T_air_h, wind_h, q_h, SW_h, LW_h, P_h, rain_flux, snow_flux, nhours) = cl_forcing

    toml_dict = LP.create_toml_dict(FT_CL)

    soil_domain = Column(; zlim = (ZMIN_CL, ZMAX_CL), nelements = NELEM_CL)
    dz_cl = (ZMAX_CL - ZMIN_CL) / NELEM_CL
    depths_cl = [(i - 0.5) * dz_cl for i in 1:NELEM_CL]   # shallow-first, m

    hydrology_cm = vanGenuchten{FT_CL}(; α = VG_ALPHA, n = VG_N)

    # Uniform mean soil composition (organic-over-mineral layering approximated
    # by its depth-weighted mean) -- same simplification the original script
    # fell back to; a depth-varying SpaceVaryingInput would be more accurate
    # but wasn't reliably working there either.
    nu_ss_om     = FT_CL((NU_SS_OM_TOP + NU_SS_OM_DEEP) / 2)
    nu_ss_quartz = FT_CL((NU_SS_QUARTZ_TOP + NU_SS_QUARTZ_DEEP) / 2)

    # Per-site (SoilGrids) porosity/K_sat -- same g/rho unit conversion as
    # K_SAT_CL, applied to prep.site_sat_hydraulic_cond instead of the fixed
    # site-329 global. VG_ALPHA/VG_N (retention curve shape) stay fixed --
    # see config.jl's comment.
    site_K_sat_cl = FT_CL(ustrip(prep.site_sat_hydraulic_cond) * 9.81 / 1000)   # m/s

    params_cl = Soil.EnergyHydrologyParameters(
        toml_dict;
        ν            = FT_CL(prep.site_saturation_moisture),
        ν_ss_om      = nu_ss_om,
        ν_ss_quartz  = nu_ss_quartz,
        ν_ss_gravel  = NU_SS_GRAVEL,
        hydrology_cm,
        K_sat        = site_K_sat_cl,
        S_s          = S_S_CL,
        θ_r          = THETA_R,
        z_0m         = FT_CL(ustrip(u"m", roughness_height)),
        z_0b         = FT_CL(ustrip(u"m", roughness_height) * 0.1),
        emissivity   = FT_CL(emissivity),
    )

    # ── Bottom BC: deep temperature from MicroclimateMapper's own
    # deep_soil_temperature series (prep.f_tannul, daily) + free drainage.
    T_deep_daily  = ustrip.(u"K", prep.f_tannul)
    T_deep_hourly = repeat(T_deep_daily; inner = 24)
    T_deep_hourly = [T_deep_hourly[1:nhours]; T_deep_hourly[nhours]]
    T_bottom_itp  = LinearInterpolation(t_sec, T_deep_hourly; extrapolation_bc = Flat())

    ref_height = ustrip(u"m", profile_heights[end])
    start_dt   = DateTime(prep.sim_start)   # actual resolved start (includes spin-up)

    top_bc = AtmosDrivenFluxBC(
        PrescribedAtmosphere(
            TimeVaryingInput(t_sec, -rain_flux),
            TimeVaryingInput(t_sec, -snow_flux),
            TimeVaryingInput(t_sec, T_air_h),
            TimeVaryingInput(t_sec, wind_h),
            TimeVaryingInput(t_sec, q_h),
            TimeVaryingInput(t_sec, P_h),
            start_dt,
            FT_CL(ref_height),
            toml_dict,
        ),
        PrescribedRadiativeFluxes(
            FT_CL,
            TimeVaryingInput(t_sec, SW_h),
            TimeVaryingInput(t_sec, LW_h),
            start_dt,
        ),
    )
    bottom_bc = WaterHeatBC(;
        water = FreeDrainage(),
        heat  = TemperatureStateBC((_p, t) -> FT_CL(T_bottom_itp(Float64(t)))),
    )

    soil_cl = Soil.EnergyHydrology{FT_CL}(;
        parameters          = params_cl,
        domain               = soil_domain,
        boundary_conditions  = (; top = top_bc, bottom = bottom_bc),
        sources              = (PhaseChange{FT_CL}(),),
    )

    # ── Initial conditions: interpolate Microclimate's 19-node profile
    # (prep.initial_st / prep.initial_sm, already obs-informed) onto
    # ClimaLand's 40-node grid by depth.
    mc_depths_m = ustrip.(u"m", depths)   # 19-node depths, m, shallow-first
    initial_st  = ustrip.(u"K", prep.initial_st)
    initial_sm  = prep.initial_sm

    function set_ic_cl!(Y, p, t0, model)
        z_field = model.domain.fields.z
        pars    = model.parameters
        depth_cl_here = vec(parent(z_field)) .* FT_CL(-1)   # positive downward

        sm_itp = LinearInterpolation(mc_depths_m, FT_CL.(initial_sm); extrapolation_bc = Flat())
        T_itp  = LinearInterpolation(mc_depths_m, FT_CL.(initial_st); extrapolation_bc = Flat())

        ϑ_l_ic = [clamp(sm_itp(d), FT_CL(pars.θ_r) + FT_CL(1e-4), FT_CL(pars.ν) - FT_CL(1e-4)) for d in depth_cl_here]
        T_ic   = FT_CL.(T_itp.(depth_cl_here))
        θi_ic  = zeros(FT_CL, length(depth_cl_here))

        θ_l_ic  = Soil.volumetric_liquid_fraction.(ϑ_l_ic, FT_CL(pars.ν), FT_CL(pars.θ_r))
        ρc_s_ic = Soil.volumetric_heat_capacity.(θ_l_ic, θi_ic, pars.ρc_ds, pars.earth_param_set)
        ρe_ic   = Soil.volumetric_internal_energy.(θi_ic, ρc_s_ic, T_ic, pars.earth_param_set)

        parent(Y.soil.ϑ_l)    .= ϑ_l_ic
        parent(Y.soil.θ_i)    .= θi_ic
        parent(Y.soil.ρe_int) .= ρe_ic
    end

    ode_algo = CTS.IMEXAlgorithm(
        CTS.ARS111(),
        CTS.NewtonsMethod(max_iters = 6, update_j = CTS.UpdateEvery(CTS.NewNewtonIteration)),
    )

    stop_dt = start_dt + Day(prep.ndays)

    diag_writer = ClimaDiagnostics.Writers.DictWriter()
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        soil_cl, start_dt, mktempdir();
        output_writer     = diag_writer,
        output_vars        = ["tsoil", "swc", "si"],
        reduction_period   = :every_dt,
        reduction_type     = :instantaneous,
        dt                 = 3600.0,
    )

    simulation = LandSimulation(
        start_dt, stop_dt, 3600.0, soil_cl;
        set_ic!         = set_ic_cl!,
        timestepper     = ode_algo,
        solver_kwargs   = (; saveat = Second(3600)),
        updateat        = Second(3600),
        diagnostics,
        user_callbacks  = (),
    )

    println("Running ClimaLand EnergyHydrology ($(prep.ndays) days)...")
    cl_solve_time = @elapsed solve!(simulation)
    @printf("  ClimaLand solver: %.2f s\n", cl_solve_time)

    # ── Unpack the DictWriter into plain time-ordered arrays.
    # ClimaCore stores columns bottom-to-top (index 1 = deepest).
    #
    # DictWriter's top-level keys are `output_short_name(diagnostic)`, which
    # (confirmed against a real run) is NOT the bare short name passed to
    # `output_vars` -- ClimaDiagnostics appends a period/reduction suffix
    # (e.g. "tsoil_1h_inst" rather than "tsoil"), following the convention
    # used across the CliMA diagnostics stack.
    
    function unpack(varname)
        available = collect(keys(diag_writer))
        matches = filter(k -> startswith(string(k), varname), available)
        isempty(matches) && error(
            "No diagnostic key starting with \"$varname\" found. " *
            "Available keys: $available")
        length(matches) > 1 && @warn "Multiple diagnostic keys match \"$varname\": $matches — using the first."
        d = diag_writer[first(matches)]
        ts = sort(collect(keys(d)))
        mat = reduce(vcat, [reshape(parent(d[t]), 1, :) for t in ts])
        return Float64.(ts), mat[:, end:-1:1]
    end

    ts_tsoil, T_cl    = unpack("tsoil")
    _,        thetal_cl = unpack("swc")
    _,        thetai_cl = unpack("si")

    depths_cl_shallow_first = depths_cl   # already shallow-first (index 1 = shallow)

    return (; cl_solve_time, t = ts_tsoil, depths_cl = depths_cl_shallow_first,
              T = T_cl, theta_l = thetal_cl, theta_i = thetai_cl, start_dt)
end
