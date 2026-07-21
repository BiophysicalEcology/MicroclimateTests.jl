# write_micropoint_inputs.jl — writes climdata.csv + params.csv for
# run_micropoint.R to read (same "Julia writes CSVs, Rscript runs, Julia
# reads CSVs back" pattern as scan_snotel/pipeline.jl's NicheMapR
# integration).
#
# climdata is built from Microclimate.jl's own solved hourly output
# (`micro_out::MicroResult`) — both models see the same air temperature/
# humidity/wind/pressure/shortwave/sky-temperature series. `winddir` is a
# constant placeholder: bare-ground point mode has no wind-direction
# dependence, but micropoint's climdata format requires the column.

const SIGMA_SB_MP = 5.670374419e-8   # W/m^2/K^4

function write_micropoint_inputs(micro_out, prep, outdir)
    nhours = prep.ndays * 24
    mkpath(outdir)

    T_air_C  = ustrip.(u"°C", uconvert.(u"°C", micro_out.reference_temperature[1:nhours]))
    RH_pct   = clamp.(micro_out.reference_humidity[1:nhours], 0.0, 1.0) .* 100.0
    P_kPa    = ustrip.(u"kPa", uconvert.(u"kPa", micro_out.pressure[1:nhours]))
    SW_h     = ustrip.(u"W/m^2", micro_out.global_radiation[1:nhours])
    difrad_h = SW_h .* micro_out.diffuse_fraction[1:nhours]
    LW_h     = SIGMA_SB_MP .* ustrip.(u"K", micro_out.sky_temperature[1:nhours]).^4
    wind_ms  = ustrip.(u"m/s", micro_out.reference_wind_speed[1:nhours])
    winddir  = zeros(nhours)

    rain_mm_day = ustrip.(u"kg/m^2", prep.f_rain)
    precip_mm_hr = zeros(nhours)
    for d in 1:prep.ndays, h in 0:23
        precip_mm_hr[(d - 1) * 24 + h + 1] = rain_mm_day[d] / 24.0
    end

    t_model = [DateTime(d) + Hour(h) for d in prep.dates_vec for h in 0:23]

    # Microclimate.jl's `hour` is LOCAL SOLAR TIME (hour 12 = true local solar
    # noon) -- SolarRadiation.jl's hour_angle() defaults longitude_correction
    # to 0 and it's never overridden anywhere in Microclimate.jl/
    # MicroclimateMapper.jl (confirmed via source). micropoint's obs_time must
    # be genuine UTC (its own docs). Without this shift, micropoint computes
    # solar geometry for the wrong instant -- confirmed empirically: Rdirdown
    # was collapsing to exactly 0 at the labeled "solar noon" on the clearest
    # days, because a site ~111°W's ~7.4h offset shifted what micropoint saw
    # as "UTC noon" to before dawn. Relabelling only (same physical instants,
    # same row order/count) -- not a data reinterpolation.
    utc_offset_hours = round(Int, -prep.lon_dd / 15)
    t_model_utc = t_model .+ Hour(utc_offset_hours)
    obs_time_str = Dates.format.(t_model_utc, dateformat"yyyy-mm-dd HH:MM:SS")

    climdata = DataFrame(
        obs_time  = obs_time_str,
        temp      = T_air_C,
        relhum    = RH_pct,
        pres      = P_kPa,
        swdown    = SW_h,
        difrad    = difrad_h,
        lwdown    = LW_h,
        windspeed = wind_ms,
        winddir   = winddir,
        precip    = precip_mm_hr,
    )
    CSV.write(joinpath(outdir, "climdata.csv"), climdata)

    heights_m = ustrip.(u"m", profile_heights)

    # Per-site soil values (SoilGrids fetch, see prepare_site) -- same unit
    # conversions as config.jl's MP_PSI_E/MP_KSAT (psi_e: J/kg -> m via /g;
    # Ksat: kg*s/m^3, exact match to createsoilc's own units, no conversion).
    site_smax  = prep.site_saturation_moisture
    site_b     = prep.site_campbell_b
    site_psi_e = ustrip(u"J/kg", prep.site_air_entry_potential) / 9.81
    site_ksat  = ustrip(prep.site_sat_hydraulic_cond)
    site_rho_kgm3 = ustrip(u"kg/m^3", uconvert(u"kg/m^3", prep.site_bulk_density))

    params = DataFrame(
        lat = [prep.lat_dd], lon = [prep.lon_dd], elev_m = [prep.elev_m],
        ndays = [prep.ndays], sim_start = [string(sim_start)], sim_end = [string(sim_end)],
        depth1_cm = [MCF_DEPTHS_CM[1]], depth2_cm = [MCF_DEPTHS_CM[2]], depth3_cm = [MCF_DEPTHS_CM[3]],
        height1_m = [heights_m[1]], height2_m = [heights_m[end]],
        smax = [site_smax], b = [site_b], psi_e = [site_psi_e], ksat = [site_ksat],
        gref = [MCF_GROUNDR], groundem = [MP_GROUNDEM], rho_kgm3 = [site_rho_kgm3],
        nlayers = [MP_NLAYERS], totaldepth = [MP_TOTALDEPTH], zm = [MP_ZM],
    )
    CSV.write(joinpath(outdir, "params.csv"), params)

    println("Wrote climdata.csv ($(nrow(climdata)) rows) + params.csv to $outdir")
    return outdir
end
