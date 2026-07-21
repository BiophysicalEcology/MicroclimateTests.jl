# run_below_canopy_r.R — runs the below-canopy vignette examples (two-stream
# shortwave, longwave, wind, leaf temperature, air temperature/humidity)
# using the real microclimlearn functions + micropoint's forestparams/
# groundparams/climdata, and writes every input/output to CSV for
# below_canopy_comparison.jl to reproduce with Microclimate.jl's own
# equivalent functions.

options(error = function() { traceback(2); quit(status = 1) })

suppressMessages(library(micropoint))
suppressMessages(library(microclimlearn))

outdir <- "C:/git/MicroclimateTests.jl/micropoint/below_canopy/outputs"
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

cat("Done -- wrote all below_canopy CSVs to", outdir, "\n")
