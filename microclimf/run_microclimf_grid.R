# run_microclimf_grid.R — grid-mode timing helper for grid_scaling.jl.
# Builds a synthetic uniform grid_n x grid_n DTM/vegp/soilc (reusing the
# already-translated soil/veg parameters from a prior comparison.jl run's
# params.csv) and times only the `runmicro()` grid solve (not the constant-
# cost point-model prep step), for a fair "how does grid solve time scale
# with N" comparison against MicroclimateMapper.jl's raster solve.
#
# <n_hours> takes the first N rows of climdata.csv (already written by
# comparison.jl) as forcing -- for a pure grid-solve-time comparison, the
# actual forcing values/source don't need to match what the Julia raster run
# used, only the row count (simulated duration) is what matters for a fair
# timing basis.
#
# Usage: Rscript run_microclimf_grid.R <outdir> <grid_n> <lat> <lon> [n_hours]

args    <- commandArgs(trailingOnly = TRUE)
outdir  <- args[1]
grid_n  <- as.integer(args[2])
lat     <- as.numeric(args[3])
lon     <- as.numeric(args[4])
n_hours <- if (length(args) >= 5) as.integer(args[5]) else NA
if (any(is.na(c(grid_n, lat, lon)))) stop("Usage: Rscript run_microclimf_grid.R <outdir> <grid_n> <lat> <lon> [n_hours]")

suppressMessages({
  library(microclimf)
  library(terra)
})

climdata <- read.csv(file.path(outdir, "climdata.csv"), stringsAsFactors = FALSE)
if (!is.na(n_hours)) climdata <- climdata[1:min(n_hours, nrow(climdata)), ]
climdata$obs_time <- as.POSIXlt(climdata$obs_time, tz = "UTC")

params <- read.csv(file.path(outdir, "params.csv"), stringsAsFactors = FALSE)
p <- as.list(params[1, ])

# ── Synthetic uniform DTM at the requested size (same CRS-construction
# approach as run_microclimf.R — see that file's header for caveats) ────────
zone <- floor((lon + 180) / 6) + 1
epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
pt <- vect(cbind(lon, lat), type = "points", crs = "EPSG:4326")
pt_utm <- project(pt, paste0("EPSG:", epsg))
xy <- crds(pt_utm)
px <- p$pixel_m
half <- (grid_n * px) / 2
e <- ext(xy[1, 1] - half, xy[1, 1] + half, xy[1, 2] - half, xy[1, 2] + half)
dtm <- rast(e, nrows = grid_n, ncols = grid_n, crs = paste0("EPSG:", epsg))
values(dtm) <- p$elev_m

mk <- function(v) { r <- dtm; values(r) <- v; r }

vegp <- list(
  pai = mk(p$pai), hgt = mk(p$hgt), x = mk(p$x), gsmax = mk(p$gsmax),
  clump = mk(p$clump), leafr = mk(p$leafr), leafd = mk(p$leafd), leaft = mk(p$leaft)
)
class(vegp) <- "vegparams"

soilc <- list(
  soiltype = mk(7), groundr = mk(p$groundr),
  Smax = mk(p$smax), Smin = mk(0.01), Ksat = mk(p$ksat), b = mk(p$b), psi_e = mk(p$psi_e),
  Vq = mk(0.3), Vm = mk(0.55), Vo = mk(0.02), Mc = mk(0.2), rho = mk(p$rho_Mgm3)
)
class(soilc) <- "soilcharac"

# ── Point model prep: constant cost regardless of grid size, so NOT
# included in the timed portion below ───────────────────────────────────────
micropoint <- runpointmodel(climdata, reqhgt = p$height1_m, dtm, vegp, soilc)

t0 <- Sys.time()
mout <- runmicro(micropoint, reqhgt = p$height1_m, vegp, soilc, dtm)
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

write.csv(data.frame(elapsed_s = elapsed_s, grid_n = grid_n),
          file.path(outdir, "microclimf_grid_timing.csv"), row.names = FALSE)
cat(sprintf("microclimf grid run: %dx%d, %.2f s (runmicro only)\n", grid_n, grid_n, elapsed_s))
