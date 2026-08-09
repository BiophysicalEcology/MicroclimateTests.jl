# run_micropoint_vegetated.R — runs micropoint's own vegetated RunModelFull
# on a real site's actual GRIDMET-derived forcing/soil (climdata.csv +
# params.csv, already written by canopy_site_check.jl via
# write_micropoint_inputs.jl) plus canopy_params.csv (this run's own
# MultilayerCanopy spec), so canopy_site_check.jl can compare real leaf/air
# profiles against Microclimate.jl's MultilayerCanopy on the *same* real
# weather -- unlike micropoint/vignette/, which used micropoint's own canned
# demo data. vegp's Eller-et-al. photosynthesis/hydraulic fields have no
# Microclimate.jl equivalent (see micropoint/README.md) and are left at
# createvegp("BET.Te")'s own defaults; only the fields that do map are
# overridden.
#
# Usage: Rscript run_micropoint_vegetated.R <outdir>

options(error = function() {
  traceback(2)
  quit(status = 1)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
if (is.na(outdir)) stop("Usage: Rscript run_micropoint_vegetated.R <outdir>")

suppressMessages(library(micropoint))

climdata <- read.csv(file.path(outdir, "climdata.csv"), stringsAsFactors = FALSE)
climdata$obs_time <- as.POSIXlt(climdata$obs_time, tz = "UTC")

p <- as.list(read.csv(file.path(outdir, "params.csv"))[1, ])

# Near-zero-but-nonzero swdown/difrad right at sunrise/sunset (zenith ~90deg,
# a real floating-point residue from Julia's own solar-geometry calc, not a
# data error) NaNs out RunModelFull's vegetated solver from that hour on --
# confirmed by bisection: the NaN onset lines up exactly with the first tiny
# positive swdown after a run of exact zeros, at completely different wind
# speeds each time, and vanishes entirely once these are zeroed. Bare-ground
# RunMicro doesn't hit this -- it's specific to the vegetated/in-canopy
# shortwave code. zenith>=89deg is a generous margin past where real
# irradiance is meaningful anyway.
solar_pre <- solarposition(p$lat, p$lon, year = climdata$obs_time)
lowsun <- solar_pre$zen >= 89
climdata$swdown[lowsun] <- 0
climdata$difrad[lowsun] <- 0
cp <- as.list(read.csv(file.path(outdir, "canopy_params.csv"))[1, ])

# ── Soil: same per-site override pattern as run_micropoint.R (bare ground) ──
soilc <- createsoilc(soiltype = "Clay loam", nlayers = p$nlayers, totalDepth = p$totaldepth)
nl <- length(soilc$Smax)
soilc$Smax  <- rep(p$smax, nl)
soilc$b     <- rep(p$b, nl)
soilc$psi_e <- rep(p$psi_e, nl)
soilc$Ksat  <- rep(p$ksat, nl)
soilc$rho   <- rep(p$rho_kgm3, nl)
soilc$gref  <- p$gref
soilc$groundem <- p$groundem

# ── Vegetation: createvegp("BET.Te") as a structural base (only field names
# that map onto MultilayerCanopy's own parameters are overridden -- Eller
# photosynthesis/hydraulic fields (Vcmx25, Tup, Tlow, Dcrit, alpha, Kxmx, hv,
# f0, fd, psi50, apsi, pTAW, rootskew) have no Microclimate.jl equivalent and
# are left at BET.Te's own defaults, same documented mismatch as vignette/) ──
vegp <- createvegp(vegtype = "BET.Te")
vegp$h    <- cp$canopy_height_m
vegp$pai  <- cp$pai
vegp$x    <- cp$leaf_angle_x
vegp$lref <- cp$leaf_reflectance
vegp$ltra <- cp$leaf_transmittance
vegp$vegem<- cp$leaf_emissivity
vegp$len  <- cp$leaf_length_m
vegp$wid  <- cp$leaf_width_m
# mwft (interception storage capacity) left at BET.Te's own default -- the
# Julia side uses NoInterception, but there's no clean way to express that
# here: mwft=0 with canopy water also starting at 0 makes micropoint's own
# wetted-fraction calc a 0/0 at hour 1, NaN-ing the whole run (confirmed).
# vignette/ leaves mwft untouched for the same reason -- same known mismatch.

n <- cp$n_layers
# Uniform PAI across layers -- Julia side uses a scalar plant_area_index
# (split evenly), not a tapered profile, so a flat vector is the faithful
# match (not PAIgeometry, which tapers -- that's vignette/'s own choice for
# its own reasons, not applicable here).
paii  <- rep(cp$pai / n, n)
# No Julia equivalent -- mid-range, matching vignette/'s own seq(0.5, 0.85)
# rather than the boundary value 1.0 (which risks hitting an edge case in
# Eller's respiration term the same way mwft=0 did above).
Lfrac <- rep(0.75, n)

# zref must be above canopy_height -- canopy_site_check.jl's own reference
# height (last entry of its `canopy_heights`, 2 m) is exactly the height
# write_micropoint_inputs.jl's climdata.csv forcing represents.
zref <- 2.0

t0 <- Sys.time()
mout <- RunModelFull(climdata, soilc, vegp, paii, Lfrac, lat = p$lat, long = p$lon, zref = zref)
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("micropoint RunModelFull (vegetated, %d layers, %d hours): %.2f s\n",
            n, nrow(climdata), elapsed_s))

# Rdirdown is beam-normal (perpendicular to the solar beam), not horizontal
# (confirmed via ?RunModelFull, same as run_micropoint.R's bare-ground
# conversion) -- convert every layer's column to horizontal so it's directly
# comparable to Microclimate.jl's own horizontal-irradiance convention.
coszen <- pmax(cos(solar_pre$zen * pi / 180), 0)
mout$Rdirdown <- mout$Rdirdown * coszen

for (nm in names(mout)) {
  write.csv(mout[[nm]], file.path(outdir, paste0("veg_", nm, ".csv")), row.names = FALSE)
}
write.csv(data.frame(elapsed_s = elapsed_s, n_layers = n, n_hours = nrow(climdata)),
          file.path(outdir, "veg_timing.csv"), row.names = FALSE)

cat(sprintf("Wrote %d vegetated output matrices + veg_timing.csv to %s\n", length(mout), outdir))
