# run_micropoint.R — reads climdata.csv/params.csv, runs micropoint::RunMicro
# per height/depth in bare-ground mode, writes micropoint_out.csv +
# micropoint_timing.csv.
#
# Usage: Rscript run_micropoint.R <outdir>

options(error = function() {
  traceback(2)
  quit(status = 1)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
if (is.na(outdir)) stop("Usage: Rscript run_micropoint.R <outdir>")

suppressMessages(library(micropoint))

climdata <- read.csv(file.path(outdir, "climdata.csv"), stringsAsFactors = FALSE)
climdata$obs_time <- as.POSIXlt(climdata$obs_time, tz = "UTC")

p <- as.list(read.csv(file.path(outdir, "params.csv"))[1, ])

soilc <- createsoilc(soiltype = "Clay loam", nlayers = p$nlayers, totalDepth = p$totaldepth)
# Smax/b/psi_e/Ksat/rho are per-soil-layer vectors (length nlayers+1, see
# ?createsoilc) -- rep() to replace every layer with our calibrated scalar
# rather than overwrite the vector with a length-1 value.
nl <- length(soilc$Smax)
soilc$Smax  <- rep(p$smax, nl)
soilc$b     <- rep(p$b, nl)
soilc$psi_e <- rep(p$psi_e, nl)
soilc$Ksat  <- rep(p$ksat, nl)
soilc$rho   <- rep(p$rho_kgm3, nl)
# gref/groundem are single scalars (not per-layer), unlike the above.
soilc$gref <- p$gref
soilc$groundem <- p$groundem

# Organic top layer, matching Microclimate.jl's depth-varying
# mineral_conductivity/mineral_heat_capacity (scan_snotel/config.jl): an
# organic litter layer in the top ~2.5cm transitioning to mineral soil by
# ~5cm. Node depths here (from micropoint's own geometricCpp(nlayers,
# totalDepth), verified against the C++ source) are approximately
# [0, 0.24, 1.21, 3.39, 7.26, ...] cm for nlayers=15/totalDepth=2m -- nodes
# 1:3 fall in the organic zone, node 4 near the transition, nodes 5+ mineral.
# thermalConductivityCpp derives conductivity from Vq/Vm/Vo composition as
# kSolid = 2.5^(Vq/Vs) * 8.8^(Vm/Vs) * 0.25^(Vo/Vs) (from the C++ source) --
# pure organic (Vo/Vs ~ 1) gives kSolid ~ 0.25 W/m/K, close to
# Microclimate.jl's 0.2 organic value (the closest reachable under this
# formula, since 0.2 is below the pure-organic floor); a 50/50 mix gives
# ~1.35, matching the transition node. heatCapacityCpp uses the same Vq/Vm/Vo
# fractions, so this one change moves both properties together, same as
# Microclimate.jl's paired conductivity/heat_capacity vectors. Deeper nodes
# are left at createsoilc's own texture default -- not forced to match
# Microclimate.jl's 2.5 W/m/K mineral value exactly, since Clay loam's raw
# composition would need an implausibly large organic fraction throughout
# the whole 2m profile to hit that under this formula (its Vq/Vm split gives
# a much higher intrinsic mineral conductivity than 2.5); this is a
# structural organic-over-mineral match, not a literal number-for-number
# conductivity calibration.
Vsolid_total <- soilc$Vq[1] + soilc$Vm[1] + soilc$Vo[1]   # preserved across createsoilc's own taper
vq_ratio <- soilc$Vq[1] / (soilc$Vq[1] + soilc$Vm[1])
vm_ratio <- soilc$Vm[1] / (soilc$Vq[1] + soilc$Vm[1])
set_organic_node <- function(i, f) {
  soilc$Vo[i] <<- f * Vsolid_total
  soilc$Vq[i] <<- (1 - f) * Vsolid_total * vq_ratio
  soilc$Vm[i] <<- (1 - f) * Vsolid_total * vm_ratio
}
set_organic_node(1, 0.98)   # ~0 cm
set_organic_node(2, 0.98)   # ~0.24 cm
set_organic_node(3, 0.98)   # ~1.21 cm
set_organic_node(4, 0.50)   # ~3.39 cm, transition

# Bare ground -- see config.jl for why (no calibration data for micropoint's
# fuller vegetated parameter set).
vegp <- NA
paii <- NA
Lfrac <- NA

# RunMicro needs ~a full year of climdata to estimate boundaryT (mean annual
# air temp) itself; passing it explicitly avoids relying on that heuristic.
boundaryT <- mean(climdata$temp)

# Picks whichever of `candidates` RunMicro actually returned -- confirmed via
# a live test run that bare-ground output uses tair/tground/tsoil/soilwater,
# kept as a safety net rather than a single hardcoded name in case vegp/NA
# handling changes upstream.
getcol <- function(mout, candidates) {
  for (nm in candidates) if (nm %in% names(mout)) return(mout[[nm]])
  stop("None of [", paste(candidates, collapse = ", "), "] found in RunMicro output. Available: ",
       paste(names(mout), collapse = ", "))
}

t0 <- Sys.time()

ThetaIni <- InitailizeWater(climdata, vegp = NA, soilc, lat = p$lat, long = p$lon,
                             zmr = p$zm, boundaryT = boundaryT)

mout_h1 <- RunMicro(climdata, reqhgt = p$height1_m, vegp, soilc, paii, Lfrac,
                     lat = p$lat, long = p$lon, zref = p$height2_m, zm = p$zm, boundaryT = boundaryT)
mout_h2 <- RunMicro(climdata, reqhgt = p$height2_m, vegp, soilc, paii, Lfrac,
                     lat = p$lat, long = p$lon, zref = p$height2_m, zm = p$zm, boundaryT = boundaryT)

mout_d1 <- RunMicro(climdata, reqhgt = -p$depth1_cm / 100, vegp, soilc, paii, Lfrac,
                     lat = p$lat, long = p$lon, zref = p$height2_m, zm = p$zm, boundaryT = boundaryT, ThetaIni = ThetaIni)
mout_d2 <- RunMicro(climdata, reqhgt = -p$depth2_cm / 100, vegp, soilc, paii, Lfrac,
                     lat = p$lat, long = p$lon, zref = p$height2_m, zm = p$zm, boundaryT = boundaryT, ThetaIni = ThetaIni)
mout_d3 <- RunMicro(climdata, reqhgt = -p$depth3_cm / 100, vegp, soilc, paii, Lfrac,
                     lat = p$lat, long = p$lon, zref = p$height2_m, zm = p$zm, boundaryT = boundaryT, ThetaIni = ThetaIni)

elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("micropoint (2 heights + 3 depths, bare ground): %.2f s\n", elapsed_s))

obs_time_out <- sprintf("%04d-%02d-%02d %02d:00:00", mout_h1$year, mout_h1$month, mout_h1$day, mout_h1$hour)

# RunMicro's Rdirdown is beam-normal (perpendicular to the solar beam), not
# horizontal (confirmed via ?RunMicro) -- convert to horizontal by
# multiplying by cos(zenith), matching the two-stream vignette's own inverse
# conversion (Rbeam <- (swdown - difrad) / cos(zen)). Without this, summing
# Rdirdown + Rdifdown blows up at low sun angles (cos(zen) -> 0) -- the same
# failure mode as microclimf's own documented radiation-spike bug.
solar <- solarposition(p$lat, p$lon, year = climdata$obs_time)
coszen <- pmax(cos(solar$zen * pi / 180), 0)
Rdirdown_horiz <- getcol(mout_h1, c("Rdirdown")) * coszen

out <- data.frame(
  obs_time  = obs_time_out,
  Tz_h1     = getcol(mout_h1, c("tair")),
  Tz_h2     = getcol(mout_h2, c("tair")),
  Tsurf     = getcol(mout_h1, c("tground")),
  T_5cm     = getcol(mout_d1, c("tsoil")),
  T_20cm    = getcol(mout_d2, c("tsoil")),
  T_50cm    = getcol(mout_d3, c("tsoil")),
  soilm_bulk = getcol(mout_d1, c("soilwater")),
  Rdirdown  = Rdirdown_horiz,
  Rdifdown  = getcol(mout_h1, c("Rdifdown")),
  Rlwdown   = getcol(mout_h1, c("Rlwdown"))
)
write.csv(out, file.path(outdir, "micropoint_out.csv"), row.names = FALSE)
write.csv(data.frame(elapsed_s = elapsed_s), file.path(outdir, "micropoint_timing.csv"), row.names = FALSE)
cat(sprintf("micropoint run complete: %.2f s, %d rows\n", elapsed_s, nrow(out)))
