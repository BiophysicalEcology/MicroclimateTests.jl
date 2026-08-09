# run_microclimf_topo.R — reads the coarse CRUCL2 forcing grid + dtmc
# (coarse elevation reference) + dtm (fine SRTM elevation) written by
# topo_scaling.jl's Stage 2, applies microclimf's own topographic correction
# (`altcorrect` elevation lapse-rate + the grid model's own slope/aspect/
# horizon-angle solar and wind-shelter correction), and writes back
# snapshots for the same variables Julia's own side produces (soil surface
# T, near-surface air T, soil moisture) for a side-by-side comparison.
#
# Soil surface T + soil moisture and near-surface air T each need their own
# `runpointmodela` call (matching the above/below-ground split already
# established in run_microclimf.R/run_microclimf_grid.R) -- only the
# near-surface (above-ground) call is included in the timed portion, mirroring
# how the other R scripts in this folder separate the "headline" timed call
# from additional untimed snapshot-gathering calls.
#
# Snow is NOT included in this first pass -- array-mode `runsnowmodel` is a
# whole further untested code path, and the site was already snow-free in
# Stage 1's July snapshot, so it's low value to chase here.
#
# Usage: Rscript run_microclimf_topo.R <topo_dir>

options(error = function() {
  traceback(2)
  quit(status = 1)
})

args <- commandArgs(trailingOnly = TRUE)
topo_dir <- args[1]
if (is.na(topo_dir)) stop("Usage: Rscript run_microclimf_topo.R <topo_dir>")

suppressMessages({
  library(microclimf)
  library(terra)
})

# Padded to 2 identical bands: runpointmodela's internal processing (e.g.
# winddir's `apply(..., MARGIN=3, .getmode)`) expects a genuine time
# dimension to summarise over; a true single-band raster collapses to 2D
# somewhere internally (R silently drops singleton dims), leaving no 3rd
# dimension for MARGIN=3 to act on. Since we only want one hour's snapshot,
# duplicating it into a trivial length-2 series doesn't change the result
# (both bands are identical) -- it just gives the array-mode code the shape
# it expects.
rd <- function(name) { r <- rast(file.path(topo_dir, name)); c(r, r) }

dtm  <- rast(file.path(topo_dir, "dtm.tif"))     # fine SRTM elevation -- the actual target grid (single-band, not padded)
dtmc <- rast(file.path(topo_dir, "dtmc.tif"))    # coarse CRUCL2 elevation -- altcorrect's lapse-rate reference (single-band)

climarrayr <- list(
  temp      = rd("coarse_temp.tif"),
  relhum    = rd("coarse_relhum.tif"),
  pres      = rd("coarse_pres.tif"),
  swdown    = rd("coarse_swdown.tif"),
  difrad    = rd("coarse_difrad.tif"),
  lwdown    = rd("coarse_lwdown.tif"),
  windspeed = rd("coarse_windspeed.tif"),
  winddir   = rd("coarse_winddir.tif"),
  precip    = rd("coarse_precip.tif")
)

# Two identical timestamps, one hour apart, matching the padded climarrayr
# above -- matches topo_scaling.jl's "day 7 (July), hour 12 (noon)" snapshot
# choice for the first band. The exact date only matters for solar geometry,
# which dtm/dtmc's real CRS/lat-long already anchors correctly.
tme <- as.POSIXlt(c("2000-07-15 12:00:00", "2000-07-15 13:00:00"), tz = "UTC")

p <- as.list(read.csv(file.path(topo_dir, "topo_params.csv"))[1, ])

mk <- function(v) { r <- dtm; values(r) <- v; r }

vegp <- list(
  pai = mk(p$pai), hgt = mk(p$hgt), x = mk(p$x), gsmax = mk(p$gsmax),
  clump = mk(p$clump), leafr = mk(p$leafr), leafd = mk(p$leafd), leaft = mk(p$leaft)
)
class(vegp) <- "vegparams"

soilc <- list(
  soiltype = mk(7), groundr = mk(p$groundr),
  Smax = mk(p$smax), Smin = mk(0.01), Ksat = mk(p$ksat), b = mk(p$b), psi_e = mk(p$psi_e),
  Vq = mk(0.3), Vm = mk(0.55), Vo = mk(0.01), Mc = mk(0.2), rho = mk(p$rho_Mgm3)
)
class(soilc) <- "soilcharac"

# dtm (the small fine SRTM patch) sits within a single CRUCL2 pixel, so
# runpointmodela sees only one overlapping climate cell and refuses to run
# in array mode ("Only one grid cell of the climate data overlaps with dtm").
# resampleclimdata upsamples climarrayr onto a finer intermediate grid
# matching dtm's own extent -- doesn't add real information (the source
# climate genuinely is uniform at sub-CRUCL2-pixel scale), just satisfies
# runpointmodela's structural requirement for more than one cell. Its
# element order matches ours (temp, relhum, pres, swdown, difrad, lwdown,
# windspeed, winddir, precip): confirmed from source, it indexes
# climarrayr[[7]]/[[8]] as windspeed/winddir internally.
climarrayr <- resampleclimdata(climarrayr, dtm)

# ── Above-ground (near-surface air T): the timed, headline call ─────────────
t0 <- Sys.time()
micropointa_air <- runpointmodela(climarrayr, tme, p$height1_m, dtm, vegp, soilc)
mout_air <- runmicro(micropointa_air, p$height1_m, vegp, soilc, dtm, dtmc, altcorrect = 1)
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("microclimf topo (near-surface air T, altcorrect=1): %.2f s\n", elapsed_s))

# ── Below-ground/surface (soil surface T + soil moisture): untimed ──────────
micropointa_srf <- runpointmodela(climarrayr, tme, 0, dtm, vegp, soilc)
mout_srf <- runmicro(micropointa_srf, 0, vegp, soilc, dtm, dtmc, altcorrect = 1)

# ── Extract the single exported hour's 2D grid for each variable ────────────
airT  <- mout_air$Tz[, , 1]
soilT <- mout_srf$Tz[, , 1]
soilm <- mout_srf$soilm[, , 1]

# Matches microclimf's own internal .rast helper (build plain, then assign
# ext/crs) rather than a matrix+extent+crs keyword constructor, which isn't
# independently confirmed to work the same way.
.torast <- function(m, template) {
  r <- rast(m)
  ext(r) <- ext(template)
  crs(r) <- crs(template)
  r
}
writeRaster(.torast(airT,  dtm), file.path(topo_dir, "r_airT.tif"),  overwrite = TRUE)
writeRaster(.torast(soilT, dtm), file.path(topo_dir, "r_soilT.tif"), overwrite = TRUE)
writeRaster(.torast(soilm, dtm), file.path(topo_dir, "r_soilm.tif"), overwrite = TRUE)

write.csv(data.frame(elapsed_s = elapsed_s), file.path(topo_dir, "r_timing.csv"), row.names = FALSE)
cat(sprintf("microclimf topo snapshots written to %s\n", topo_dir))
