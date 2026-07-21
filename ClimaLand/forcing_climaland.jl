# forcing_climaland.jl — build ClimaLand TimeVaryingInputs from Microclimate.jl's
# own *solved* hourly output (`MicroResult`), so both models are driven by
# literally the same air temperature/humidity/wind/shortwave/sky-temperature
# series — instead of independently re-deriving an approximate diurnal curve
# (what the original Dropbox scan329_comparison.jl did, and flagged as such).
#
# `MicroResult` (see Microclimate.jl's `outputs.jl`) already carries, at
# hourly resolution, exactly what Microclimate.jl used internally:
#   reference_temperature   (K)
#   reference_humidity      (fractional, 0-1)
#   reference_wind_speed    (m/s)
#   pressure                (Pa)
#   global_radiation        (W/m^2)      -- downward shortwave
#   sky_temperature         (K)          -- downward longwave = sigma*T_sky^4
#
# The one thing MicroResult does *not* expose is a realized hourly
# precipitation/snowfall water-equivalent series (evaporation/transpiration/
# drainage are computed inside Microclimate.jl's ODE right-hand side each
# step but discarded, and `snow_fall` is a depth-rate, not water-equivalent)
# — so rain/snow flux is still derived by hand here, from the same daily
# total (`prep.f_rain`) and the same rain/snow temperature threshold
# (`snow_temp_threshold`, from scan_snotel/config.jl) that Microclimate.jl's
# own SnowModel uses, spread uniformly across each day's 24 hours.

const SIGMA_SB = 5.670374419e-8   # W/m^2/K^4, Stefan-Boltzmann

function build_climaland_forcing(micro_out, prep)
    ndays  = prep.ndays
    nhours = ndays * 24

    T_air_h = ustrip.(u"K", micro_out.reference_temperature[1:nhours])
    RH_h    = clamp.(micro_out.reference_humidity[1:nhours], 0.0, 1.0)
    wind_h  = ustrip.(u"m/s", micro_out.reference_wind_speed[1:nhours])
    P_h     = ustrip.(u"Pa", micro_out.pressure[1:nhours])
    SW_h    = ustrip.(u"W/m^2", micro_out.global_radiation[1:nhours])
    LW_h    = SIGMA_SB .* ustrip.(u"K", micro_out.sky_temperature[1:nhours]).^4

    # Specific humidity from RH, using the same saturation-vapour-pressure
    # formulation Microclimate.jl defaults to (GoffGratch) for consistency.
    q_h = similar(T_air_h)
    for i in eachindex(T_air_h)
        e_sat = ustrip(u"Pa", FluidProperties.vapour_pressure(FluidProperties.GoffGratch(), T_air_h[i] * u"K"))
        e     = RH_h[i] * e_sat
        q_h[i] = 0.622 * e / (P_h[i] - 0.378 * e)
    end

    # ── Rain/snow split: daily total spread uniformly across 24h, split by
    # hourly air temperature vs the same snow threshold Microclimate.jl uses.
    rain_flux = zeros(nhours)
    snow_flux = zeros(nhours)
    snow_thresh_K = ustrip(u"K", snow_temp_threshold)
    rain_daily_m  = ustrip.(u"m", prep.f_rain ./ (1000.0u"kg/m^3"))   # kg/m^2/day -> m/day (rho_water = 1000 kg/m^3)
    for d in 1:ndays
        precip_rate = rain_daily_m[d] / 86400.0   # m/s, uniform over the day
        for h in 0:23
            idx = (d - 1) * 24 + h + 1
            if T_air_h[idx] < snow_thresh_K
                snow_flux[idx] = precip_rate
            else
                rain_flux[idx] = precip_rate
            end
        end
    end

    # Extend all arrays by one point so TimeVaryingInput covers the final
    # simulation time (nhours*3600 s) -- same pattern as the original script.
    T_air_h   = [T_air_h;   T_air_h[end]]
    wind_h    = [wind_h;    wind_h[end]]
    q_h       = [q_h;       q_h[end]]
    SW_h      = [SW_h;      SW_h[end]]
    LW_h      = [LW_h;      LW_h[end]]
    P_h       = [P_h;       P_h[end]]
    rain_flux = [rain_flux; rain_flux[end]]
    snow_flux = [snow_flux; snow_flux[end]]
    t_sec     = Float64.(0:3600:nhours*3600)

    return (; t_sec, T_air_h, wind_h, q_h, SW_h, LW_h, P_h, rain_flux, snow_flux, nhours)
end
