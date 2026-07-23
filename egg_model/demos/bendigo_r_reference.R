# R reference run of the egg-development state machine (ported verbatim from
# seasonal_hatch_forecast2.R's per-hour loop), driven by the same Bendigo
# soil_temperature/soil_humidity/soil_water_potential series the Julia model
# uses (bendigo_feb1_driver.csv, median-filtered, Feb 1 2020 lay date, 150 days).
# Isolates the comparison to the egg-development logic itself, not differences
# between the R and Julia microclimate models.

library(deSolve)

source("C:/git/NicheMapR/R/VAPPRS.R")
source("C:/git/NicheMapR/R/DRYAIR.R")
source("C:/git/NicheMapR/R/WETAIR.R")
source("C:/git/NicheMapR/R/egg_water.R")

driver <- read.csv("c:/git/MicroclimateTests.jl/egg_model/demos/bendigo_feb1_driver.csv")
n <- nrow(driver)

# real fitted parameters (seasonal_hatch_forecast2.R's egg_pars, Arrhenius_just_eggs.R's arr_pars)
desic_tol        <- 0.6
psi_e_init       <- -709.4682
K_e_base         <- 2.347802e-09
f_air            <- 0.5
pct_wet_egg_init <- 0.35
egg_mass_init    <- 0.0036
min_egg_mass     <- 0.0026
max_egg_mass     <- 2 * egg_mass_init  # placeholder turgid-mass ceiling (membrane tension counters uptake); needs a real estimate
shape_b          <- 0.69 / 1.82
P_e   <- -0.5052209
b     <- 1.41005
K_sat <- 0.003733307
elev  <- 0
vel   <- 0.001
spec_hyd <- 0.000304

T_A  <- 6641.6175; T_AL <- 33600.0; T_AH <- 48000.0
T_L  <- 289.15; T_H <- 314.65; T_REF <- 301.65; kdot_ref <- 1 / 17.4  # per day

ArrFunc5 <- function(x, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref) {
  exp(T_A / T_REF - T_A / x) *
    (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) /
    (1 + exp(T_AL / x - T_AL / T_L) + exp(T_AH / T_H - T_AH / x)) *
    kdot_ref
}

cold_thresh <- -15           # diapause disabled in this scenario (matches Julia's cold_hour_threshold=1000hr, effectively off)
cold_hours_thresh <- 1000
diapause_hours_thresh <- 240

# hourly interpolators from the exported driver
times_s <- driver$hour * 3600
Tsoilf       <- approxfun(times_s, driver$soil_T_C, rule = 2)
Tbsf         <- approxfun(times_s, driver$soil_T_C, rule = 2)   # SoilTemperatureEquals: egg temp == soil temp
RHsoilf      <- approxfun(times_s, driver$soil_RH, rule = 2)
PSIsoilf_raw <- approxfun(times_s, driver$soil_psi_Jkg, rule = 2)
# egg_water.R only clamps psi_s to <=0 (`psi_s <- min(0, PSIsoilf(t))`), so a
# near-saturated/rain spell (psi_s -> 0) sends k_s = K_sat*(P_e/psi_s)^n to
# infinity -- clamp here to <=P_e (air entry potential, "as wet as saturated
# soil gets") instead, matching the defensive floor already added on the
# Julia side (hydric.jl's soil_hydraulic_conductivity).
PSIsoilf <- function(t) min(P_e, PSIsoilf_raw(t))
Tbs      <- driver$soil_T_C

dev <- cold_hours <- diapause_hrs <- numeric(n)
diapause_pot <- rep(TRUE, n)
in_diapause <- Q1 <- Q2 <- rep(FALSE, n)
Ww_e <- rep(egg_mass_init, n)
psi_e <- rep(psi_e_init, n)
egg_desic <- rep(NA, n)
cold <- 0; m_init <- egg_mass_init; m_max <- egg_mass_init
hatchdate <- NA

for (i in 2:n) {
  if (Tbs[i] < cold_thresh) cold <- cold + 1
  cold_hours[i] <- cold
  if (cold_hours[i] > cold_hours_thresh || diapause_hrs[i] > diapause_hours_thresh)
    diapause_pot[i] <- FALSE
  if (dev[i-1] > 0.25 && dev[i-1] < 0.30 && diapause_pot[i]) {
    diapause_hrs <- diapause_hrs + 1
    in_diapause[i] <- TRUE
  }

  K_e     <- if (dev[i-1] > 0.25) K_e_base * 3 else 0
  pct_wet <- if (dev[i-1] > 0.45) pct_wet_egg_init * 100 else pct_wet_egg_init
  if (min(-0.001, PSIsoilf(i * 3600)) < psi_e[i-1]) K_e <- 0

  indata <- list(m_init = m_init / 1000, psi_e_init = psi_e[i-1],
                 shape = shape_b, f_air = f_air, K_e = K_e,
                 spec_hyd = spec_hyd, pct_wet = pct_wet,
                 P_e = P_e, b = b, K_sat = K_sat, elev = elev, vel = vel,
                 Tsoilf = Tsoilf, Tbf = Tbsf, PSIsoilf = PSIsoilf, RHsoilf = RHsoilf)
  egg.water.out <- as.data.frame(ode(
    y = c(m_init / 1000, psi_e[i-1]),
    times = c(i * 3600, i * 3600 + 3600),
    func = egg_water, parms = indata))
  Ww_e[i]      <- min(max_egg_mass, max(min_egg_mass, egg.water.out[2, 2] * 1000))
  psi_e[i]     <- min(-0.001, egg.water.out[2, 3])
  egg_desic[i] <- min(1, (Ww_e[i] - egg_mass_init) / (m_max - egg_mass_init))
  if (is.na(egg_desic[i])) egg_desic[i] <- 1
  m_init <- Ww_e[i]
  if (m_init > m_max) m_max <- m_init

  Q1[i] <- egg_desic[i] < desic_tol && dev[i-1] > 0.25 && dev[i-1] < 0.30
  Q2[i] <- !in_diapause[i] && egg_desic[i] < desic_tol && dev[i-1] > 0.45 && dev[i-1] < 0.50

  dev[i] <- if (Q1[i] || Q2[i] || in_diapause[i]) {
    dev[i-1]
  } else {
    dev[i-1] + ArrFunc5(Tbs[i] + 273.15, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref) / 24
  }
  if (dev[i] >= 1 && is.na(hatchdate)) hatchdate <- i

  if (i %% 200 == 0) cat(sprintf("hour=%d dev=%.4f mass=%.4fmg psi_e=%.1f m_max=%.4fmg Q1=%s Q2=%s\n",
                                   i, dev[i], Ww_e[i], psi_e[i], m_max, Q1[i], Q2[i]))
}

cat("\nhatch hour:", hatchdate, ifelse(is.na(hatchdate), paste("(did not hatch in", n, "hours)"), paste("=", round(hatchdate/24,1), "days")), "\n")
cat("final dev:", tail(dev,1), " final mass:", tail(Ww_e,1), "mg  m_max:", m_max, "mg\n")

out <- data.frame(hour = 0:(n-1), dev = dev, mass_mg = Ww_e, psi_e = psi_e, egg_desic = egg_desic, Q1 = Q1, Q2 = Q2, in_diapause = in_diapause)
write.csv(out, "c:/git/MicroclimateTests.jl/egg_model/demos/bendigo_r_reference_output.csv", row.names = FALSE)
