# run_below_canopy_r.R — runs the below-canopy vignette examples (two-stream
# shortwave, longwave, wind, leaf temperature, air temperature/humidity)
# using the real microclimlearn functions + micropoint's forestparams/
# groundparams/climdata, and writes every input/output to CSV for
# below_canopy_comparison.jl to reproduce with Microclimate.jl's own
# equivalent functions.

options(error = function() { traceback(2); quit(status = 1) })

suppressMessages(library(micropoint))
suppressMessages(library(microclimlearn))

outdir <- "C:/git/MicroclimateTests.jl/micropoint/module_comparison/outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

lat <- 49.96807
long <- -5.215668

vegp <- micropoint::forestparams
groundp <- micropoint::groundparams
write.csv(as.data.frame(unclass(vegp)), file.path(outdir, "vegp.csv"), row.names = FALSE)
write.csv(as.data.frame(unclass(groundp)), file.path(outdir, "groundp.csv"), row.names = FALSE)

# ── 1: two-stream shortwave ─────────────────────────────────────────────────
n1 <- 1000
paii1 <- PAIgeometry(PAI = vegp$pai, skew = vegp$skew, spread = vegp$spread, n = n1)
z1 <- ((1:n1) / n1) * vegp$h
write.csv(data.frame(z = z1, paii = paii1), file.path(outdir, "01_paii.csv"), row.names = FALSE)

weather <- micropoint::climdata
hr <- 4094
swdown <- weather$swdown[hr]
difrad <- weather$difrad[hr]
tme <- as.POSIXlt(weather$obs_time[hr], tz = "UTC")
solp <- sunposition(tme, lat = lat, long = long)
write.csv(data.frame(swdown = swdown, difrad = difrad, zen = solp$zen, azi = solp$azi),
    file.path(outdir, "01_inputs.csv"), row.names = FALSE)

tsmod <- twostream(vegp, groundp, paii1, swdown, difrad, lat, long, solp)
write.csv(data.frame(z = z1, Rdirdown = tsmod$Rdirdown, Rdifdown = tsmod$Rdifdown,
    Rswup = tsmod$Rswup, Rleafabs = tsmod$Rleafabs), file.path(outdir, "01_twostream.csv"), row.names = FALSE)

# ── 2: below-canopy longwave ────────────────────────────────────────────────
n2 <- 100
paii2 <- PAIgeometry(PAI = 4, skew = 7, spread = 70, n = n2)
z2 <- (1:n2) / 10
tleaf2 <- 3 * cos(1 - z2 / 10) + 13.5
lwb <- longwavebelow(paii2, lwdown = 300, tleaf2, tground = 15, vegem = 0.97, groundem = 0.97)
write.csv(data.frame(z = z2, paii = paii2, tleaf = tleaf2, Rlwdown = lwb$Rlwdown,
    Rlwup = lwb$Rlwup, RlwLabs = lwb$RlwLabs), file.path(outdir, "02_longwavebelow.csv"), row.names = FALSE)

# ── 3: below-canopy wind ────────────────────────────────────────────────────
n3 <- 100
paii_skew <- PAIgeometry(PAI = 4, skew = 3, spread = 70, n = n3)
paii_unif <- rep(4 / n3, n3)
uh <- 3
wind_skew <- windprofile_below(hgt = 25, paii_skew, uh)
wind_unif <- windprofile_below(hgt = 25, paii_unif, uh)
z3 <- ((1:n3) / n3) * 25
write.csv(data.frame(z = z3, paii_skew = paii_skew, paii_unif = paii_unif,
    wind_skew = wind_skew, wind_unif = wind_unif), file.path(outdir, "03_windprofile.csv"), row.names = FALSE)

# ── 4: leaf temperature ─────────────────────────────────────────────────────
n4 <- 100
vegp4 <- createplant_inputs("BDT")
paii4 <- PAIgeometry(PAI = vegp4$pai, skew = 7, spread = 70, n = n4)
z4 <- ((1:n4) / n4) * vegp4$h

tleaf4 <- rep(26, n4)
tair4 <- rep(25, n4)
tground4 <- 25
rh4 <- rep(80, n4)

tme4 <- as.POSIXlt(0, origin = "2026-05-21 12:00", tz = "UTC")
solp4 <- sunposition(tme4, lat = lat, long = long)
tsmod4 <- twostream(vegp4, groundp, paii4, swdown = 839, difrad = 223, lat = lat, long = long, solp4)
lwb4 <- longwavebelow(paii4, lwdown = 300, tleaf4, tground4, vegem = 0.97, groundem = 0.97)
Rabs4 <- lwb4$RlwLabs + tsmod4$Rleafabs

uh4 <- 2
uz4 <- windprofile_below(vegp4$h, paii4, uh4)
dT4 <- abs(tleaf4 - tair4)
rHa4 <- leafresistance(tair4, dT4, uz4, vegp4$len, vegp4$wid, vegp4$x)

tsmod4_par <- twostream(vegp4, groundp, paii4, swdown = 839, difrad = 223, lat = lat, long = long, solp4, PAR = TRUE)
gs_sun4 <- stomatalcond_calc(430, tsmod4_par$Rsun, tair4, tleaf4, rh4, 101.3, -0.5, vegp4, z4)
gs_shade4 <- stomatalcond_calc(430, tsmod4_par$Rshade, tair4, tleaf4, rh4, 101.3, -0.5, vegp4, z4)
gs4 <- gs_sun4 * tsmod4_par$sunlitfrac + (1 - tsmod4_par$sunlitfrac) * gs_shade4
ph4 <- rhohair(tair4, 101.3)
rS4 <- ph4 / gs4

tleaf4_solved <- SolveEnergyBalance(tair4, 101.3, rh4, Rabs4, rHa4, rS4, G = 0, em = 0.97,
    method = "Penman", iters = 4)
write.csv(data.frame(z = z4, Rabs = Rabs4, rHa = rHa4, gs = gs4, rS = rS4, uz = uz4, tleaf = tleaf4_solved),
    file.path(outdir, "04_leaftemperature.csv"), row.names = FALSE)

# ── 5: below-canopy air temperature and humidity ────────────────────────────
n5 <- 100
vegp5 <- createplant_inputs("BDT")
paii5 <- PAIgeometry(PAI = vegp5$pai, skew = 7, spread = 70, n = n5)
write.csv(data.frame(paii = paii5), file.path(outdir, "05_paii.csv"), row.names = FALSE)
Ca5 <- Cafromyear(2017)

can5 <- solve_wholecanopy(weather, vegp5, groundp, hr = hr, lat = lat, long = long,
    zref = 22, G = 20, theta = 0.2, soilrh = 87, Ca = Ca5)
write.csv(data.frame(uf = can5$uf, uh = can5$uh, th = can5$th, TL = can5$TL),
    file.path(outdir, "05_wholecanopy.csv"), row.names = FALSE)

lwdown5 <- weather$lwdown[hr]
swdown5 <- weather$swdown[hr]
difrad5 <- weather$difrad[hr]
pk5 <- weather$pres[hr]
uf5 <- can5$uf
uh5 <- can5$uh
th5 <- can5$th
TL5 <- can5$TL
tground5 <- 22
psi_r5 <- -0.5
soilrh5 <- 87

tleaf5 <- rep(th5 + 1, n5)
tair5 <- rep(th5, n5)
rh5 <- rep(weather$relhum[hr], n5)

tme5 <- as.POSIXlt(weather$obs_time[hr], tz = "UTC")
solp5 <- sunposition(tme5, lat, long)
tsmod5 <- twostream(vegp5, groundp, paii5, swdown5, difrad5, lat, long, solp5)

uz5 <- windprofile_below(vegp5$h, paii5, uh5)
z5 <- ((1:n5) / n5) * vegp5$h

error5 <- 1e99
aitkt <- list(oldv = tair5, newv = tair5, z = z5, hgt = vegp5$h)
ea5 <- satvap(tair5) * (rh5 / 100)
aitkr <- list(oldv = ea5, newv = ea5, z = z5, hgt = vegp5$h)
itr5 <- 0

while (error5 > 1e-2 && itr5 < 200) {
  lwb5 <- longwavebelow(paii5, lwdown5, tleaf5, tground5, 0.97, 0.97)
  Rabs5 <- lwb5$RlwLabs + tsmod5$Rleafabs

  dT5 <- abs(tleaf5 - tair5)
  rHa5 <- leafresistance(tair5, dT5, uz5, vegp5$len, vegp5$wid, vegp5$x)

  tsmod5_par <- twostream(vegp5, groundp, paii5, swdown5, difrad5, lat, long, solp5, PAR = TRUE)
  gs_sun5 <- stomatalcond_calc(Ca5, tsmod5_par$Rsun, tair5, tleaf5, rh5, pk5, psi_r5, vegp5, z5)
  gs_shade5 <- stomatalcond_calc(Ca5, tsmod5_par$Rshade, tair5, tleaf5, rh5, pk5, psi_r5, vegp5, z5)
  gs5 <- gs_sun5 * tsmod5_par$sunlitfrac + (1 - tsmod5_par$sunlitfrac) * gs_shade5
  rS5 <- rhohair(tair5, pk5) / gs5

  tleaf5 <- SolveEnergyBalance(tair5, pk5, rh5, Rabs5, rHa5, rS5, G = 0, em = 0.97,
      method = "Penman", iters = 4)

  ph5 <- rhohair(tair5, pk5)
  cp5 <- cpair(tair5)
  Hz5 <- (ph5 * cp5 / rHa5) * (tleaf5 - tair5)
  tairn5 <- temp_below(vegp5$h, paii5, TL5, uf5, pk5, Hz5, tair5, tleaf5, tground5)

  rV5 <- rHa5 + rS5
  es5 <- satvap(tleaf5)
  ea5 <- satvap(tair5) * (rh5 / 100)
  Lz5 <- ((latvap(tleaf5) * ph5) / (rV5 * pk5)) * (es5 - ea5)
  rhn5 <- relhum_below(vegp5$h, paii5, TL5, uf5, pk5, Lz5, rh5, tair5, tleaf5, tground5, soilrh5)
  ean5 <- satvap(tairn5) * (rhn5 / 100)

  error5 <- max(abs(tair5 - tairn5))
  aitkt$oldv <- tair5; aitkr$oldv <- ea5
  aitkt$newv <- tairn5; aitkr$newv <- ean5
  aitkt <- aitken_weightdif(aitkt)
  aitkr <- aitken_weightdif(aitkr)
  tair5 <- aitkt$newv
  rh5 <- (aitkr$newv / satvap(tair5)) * 100

  itr5 <- itr5 + 1
}
cat(sprintf("Section 5: converged in %d iterations, error=%.5f\n", itr5, error5))

es5 <- satvap(tleaf5)
ea5 <- satvap(tair5) * (rh5 / 100)
write.csv(data.frame(z = z5, tleaf = tleaf5, tair = tair5, rh = rh5, es = es5, ea = ea5),
    file.path(outdir, "05_airtemphumidity.csv"), row.names = FALSE)

# ── 5c: zero-stomatal-conductance, ground-temperature sweep ─────────────────
# Isolates the heat-only (sensible-flux-driven) below-canopy transport from
# stomatal/latent-heat coupling -- gs forced near-zero (not exactly 0, to
# avoid rS = rhohair/gs blowing up to literal Inf) so rS is astronomically
# large and Lz5 (latent) collapses to ~0, leaving temp_below purely driven
# by Hz5/tground5. Reuses section 5's already-computed uf5/uh5/th5/TL5/
# tsmod5/uz5/z5 (wind & radiation don't depend on tground or gs).
gs5_zero <- 1e-9
tground_sweep <- c(5, 10, 15, 20, 25, 30, 35)

for (tg in tground_sweep) {
  tleaf5c <- rep(th5 + 1, n5)
  tair5c <- rep(th5, n5)
  rh5c <- rep(weather$relhum[hr], n5)
  ea5c <- satvap(tair5c) * (rh5c / 100)
  aitkt5c <- list(oldv = tair5c, newv = tair5c, z = z5, hgt = vegp5$h)
  aitkr5c <- list(oldv = ea5c, newv = ea5c, z = z5, hgt = vegp5$h)
  error5c <- 1e99
  itr5c <- 0

  while (error5c > 1e-2 && itr5c < 200) {
    lwb5c <- longwavebelow(paii5, lwdown5, tleaf5c, tg, 0.97, 0.97)
    Rabs5c <- lwb5c$RlwLabs + tsmod5$Rleafabs

    dT5c <- abs(tleaf5c - tair5c)
    rHa5c <- leafresistance(tair5c, dT5c, uz5, vegp5$len, vegp5$wid, vegp5$x)
    rS5c <- rhohair(tair5c, pk5) / gs5_zero

    tleaf5c <- SolveEnergyBalance(tair5c, pk5, rh5c, Rabs5c, rHa5c, rS5c, G = 0, em = 0.97,
        method = "Penman", iters = 4)

    ph5c <- rhohair(tair5c, pk5)
    cp5c <- cpair(tair5c)
    Hz5c <- (ph5c * cp5c / rHa5c) * (tleaf5c - tair5c)
    tairn5c <- temp_below(vegp5$h, paii5, TL5, uf5, pk5, Hz5c, tair5c, tleaf5c, tg)

    rV5c <- rHa5c + rS5c
    es5c <- satvap(tleaf5c)
    ea5c <- satvap(tair5c) * (rh5c / 100)
    Lz5c <- ((latvap(tleaf5c) * ph5c) / (rV5c * pk5)) * (es5c - ea5c)
    rhn5c <- relhum_below(vegp5$h, paii5, TL5, uf5, pk5, Lz5c, rh5c, tair5c, tleaf5c, tg, soilrh5)
    ean5c <- satvap(tairn5c) * (rhn5c / 100)

    error5c <- max(abs(tair5c - tairn5c))
    aitkt5c$oldv <- tair5c; aitkr5c$oldv <- ea5c
    aitkt5c$newv <- tairn5c; aitkr5c$newv <- ean5c
    aitkt5c <- aitken_weightdif(aitkt5c)
    aitkr5c <- aitken_weightdif(aitkr5c)
    tair5c <- aitkt5c$newv
    rh5c <- (aitkr5c$newv / satvap(tair5c)) * 100

    itr5c <- itr5c + 1
  }
  cat(sprintf("Section 5c (tground=%g): converged in %d iterations, error=%.5f\n", tg, itr5c, error5c))

  # Diagnostic: ground-most resistance rHa at i=1 -- matches temp_below's own
  # internal computation exactly (single term, since sumRH only has RH[1] at
  # i=1). Direct analog of Julia's buffers.air_profile.resistance_to_ground[n]
  # (n=ground-most in Julia's opposite indexing convention), for the
  # ground-coupling-gain comparison.
  ow_diag <- uf5 * (0.75 + 0.5 * cos(pi * (1 - z5[1] / vegp5$h)))
  KH_diag <- TL5 * ow_diag^2
  rHa_ground <- max(2, (1 / KH_diag) * (vegp5$h / n5))
  cat(sprintf("  ground_resistance(R rHa_ground) = %.4f s/m\n", rHa_ground))

  # Diagnostic: the actual GT value fed into H at i=1 (ground-most), using
  # the converged tair5c -- includes rHa_ground AND the dz factor together,
  # so this is the real operational weight R gives the ground term, not
  # just the resistance in isolation.
  dz_diag <- vegp5$h / n5
  GT_ground <- (rhohair(tair5c[1], pk5) * cpair(tair5c[1]) / rHa_ground) * (tg - tair5c[1]) * dz_diag
  cat(sprintf("  ground_flux(R GT, i=1) = %.4f  (dz=%.4f, tair[1]=%.3f)\n", GT_ground, dz_diag, tair5c[1]))

  # Converged-state Hz (raw, pre-paii-weighted, matching temp_below's own SS
  # = paii*Hz split) -- exported so Julia can inject the *exact same* SS
  # source terms into canopy_air_profile! directly, single-pass, isolating
  # the transport formula itself from any leaf-energy-balance formula
  # differences between the two models.
  dT5c_final <- abs(tleaf5c - tair5c)
  rHa5c_final <- leafresistance(tair5c, dT5c_final, uz5, vegp5$len, vegp5$wid, vegp5$x)
  ph5c_final <- rhohair(tair5c, pk5)
  cp5c_final <- cpair(tair5c)
  Hz5c_final <- (ph5c_final * cp5c_final / rHa5c_final) * (tleaf5c - tair5c)

  write.csv(data.frame(z = z5, tleaf = tleaf5c, tair = tair5c, rh = rh5c, Hz = Hz5c_final, paii = paii5, uf = uf5, TL = TL5, th = th5),
      file.path(outdir, sprintf("05c_gs0_tg%g_micropoint.csv", tg)), row.names = FALSE)
}

# ── 5d: real (non-zero) stomatal conductance, ground-temperature sweep ──────
# Same sweep as 5c, but with actual stomatalcond_calc-derived gs (light/CO2/
# VPD-responsive, as section 5 itself uses) instead of forced near-zero --
# shows how much the stomatal/photosynthesis-model difference (Julia's flat
# PrescribedStomatalConductance vs R's stomatalcond_calc) adds on top of the
# already-verified transport formula. Exports both Hz and Lz (latent flux,
# W/m^2) so raupach_formula_isolation.jl can inject the full heat+vapor
# source terms.
for (tg in tground_sweep) {
  tleaf5d <- rep(th5 + 1, n5)
  tair5d <- rep(th5, n5)
  rh5d <- rep(weather$relhum[hr], n5)
  ea5d <- satvap(tair5d) * (rh5d / 100)
  aitkt5d <- list(oldv = tair5d, newv = tair5d, z = z5, hgt = vegp5$h)
  aitkr5d <- list(oldv = ea5d, newv = ea5d, z = z5, hgt = vegp5$h)
  error5d <- 1e99
  itr5d <- 0

  while (error5d > 1e-2 && itr5d < 200) {
    lwb5d <- longwavebelow(paii5, lwdown5, tleaf5d, tg, 0.97, 0.97)
    Rabs5d <- lwb5d$RlwLabs + tsmod5$Rleafabs

    dT5d <- abs(tleaf5d - tair5d)
    rHa5d <- leafresistance(tair5d, dT5d, uz5, vegp5$len, vegp5$wid, vegp5$x)

    tsmod5d_par <- twostream(vegp5, groundp, paii5, swdown5, difrad5, lat, long, solp5, PAR = TRUE)
    gs_sun5d <- stomatalcond_calc(Ca5, tsmod5d_par$Rsun, tair5d, tleaf5d, rh5d, pk5, psi_r5, vegp5, z5)
    gs_shade5d <- stomatalcond_calc(Ca5, tsmod5d_par$Rshade, tair5d, tleaf5d, rh5d, pk5, psi_r5, vegp5, z5)
    gs5d <- gs_sun5d * tsmod5d_par$sunlitfrac + (1 - tsmod5d_par$sunlitfrac) * gs_shade5d
    rS5d <- rhohair(tair5d, pk5) / gs5d

    tleaf5d <- SolveEnergyBalance(tair5d, pk5, rh5d, Rabs5d, rHa5d, rS5d, G = 0, em = 0.97,
        method = "Penman", iters = 4)

    ph5d <- rhohair(tair5d, pk5)
    cp5d <- cpair(tair5d)
    Hz5d <- (ph5d * cp5d / rHa5d) * (tleaf5d - tair5d)
    tairn5d <- temp_below(vegp5$h, paii5, TL5, uf5, pk5, Hz5d, tair5d, tleaf5d, tg)

    rV5d <- rHa5d + rS5d
    es5d <- satvap(tleaf5d)
    ea5d <- satvap(tair5d) * (rh5d / 100)
    Lz5d <- ((latvap(tleaf5d) * ph5d) / (rV5d * pk5)) * (es5d - ea5d)
    rhn5d <- relhum_below(vegp5$h, paii5, TL5, uf5, pk5, Lz5d, rh5d, tair5d, tleaf5d, tg, soilrh5)
    ean5d <- satvap(tairn5d) * (rhn5d / 100)

    error5d <- max(abs(tair5d - tairn5d))
    aitkt5d$oldv <- tair5d; aitkr5d$oldv <- ea5d
    aitkt5d$newv <- tairn5d; aitkr5d$newv <- ean5d
    aitkt5d <- aitken_weightdif(aitkt5d)
    aitkr5d <- aitken_weightdif(aitkr5d)
    tair5d <- aitkt5d$newv
    rh5d <- (aitkr5d$newv / satvap(tair5d)) * 100

    itr5d <- itr5d + 1
  }
  cat(sprintf("Section 5d (tground=%g): converged in %d iterations, error=%.5f, mean gs=%.4f\n",
      tg, itr5d, error5d, mean(gs5d)))

  # Converged-state Hz and Lz (raw, pre-paii-weighted), recomputed once more
  # from the converged state for export, same pattern as section 5c.
  dT5d_final <- abs(tleaf5d - tair5d)
  rHa5d_final <- leafresistance(tair5d, dT5d_final, uz5, vegp5$len, vegp5$wid, vegp5$x)
  ph5d_final <- rhohair(tair5d, pk5)
  cp5d_final <- cpair(tair5d)
  Hz5d_final <- (ph5d_final * cp5d_final / rHa5d_final) * (tleaf5d - tair5d)
  rS5d_final <- rhohair(tair5d, pk5) / gs5d
  rV5d_final <- rHa5d_final + rS5d_final
  es5d_final <- satvap(tleaf5d)
  ea5d_final <- satvap(tair5d) * (rh5d / 100)
  Lz5d_final <- ((latvap(tleaf5d) * ph5d_final) / (rV5d_final * pk5)) * (es5d_final - ea5d_final)

  write.csv(data.frame(z = z5, tleaf = tleaf5d, tair = tair5d, rh = rh5d, Hz = Hz5d_final, Lz = Lz5d_final,
      gs = gs5d, paii = paii5, uf = uf5, TL = TL5, th = th5),
      file.path(outdir, sprintf("05d_gsreal_tg%g_micropoint.csv", tg)), row.names = FALSE)
}

cat("Done -- wrote all below_canopy CSVs to", outdir, "\n")
