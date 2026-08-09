# run_micropoint_ozflux_vegetated.R — runs micropoint's own vegetated
# RunModelFull on OzFlux tower forcing (climdata.csv/params.csv,
# already written by write_ozflux_micropoint_inputs.jl from the same
# comparisons/ozflux `result` Microclimate.jl solved) plus canopy_params.csv
# (that run's own MultilayerCanopy spec) and canopy_layers.csv (micropoint's
# own evenly-spaced PAI grid, shaped to match Julia's PAI profile). Adapted
# from micropoint/run_micropoint_vegetated.R (same vegp mapping) and
# micropoint/run_micropoint.R (organic-top-layer soil override), but reads
# zref/PAI shape from files instead of hardcoding.
#
# Usage: Rscript run_micropoint_ozflux_vegetated.R <outdir>

options(error = function() {
  traceback(2)
  quit(status = 1)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
if (is.na(outdir)) stop("Usage: Rscript run_micropoint_ozflux_vegetated.R <outdir>")

suppressMessages(library(micropoint))

climdata <- read.csv(file.path(outdir, "climdata.csv"), stringsAsFactors = FALSE)
climdata$obs_time <- as.POSIXlt(climdata$obs_time, tz = "UTC")

p  <- as.list(read.csv(file.path(outdir, "params.csv"))[1, ])
cp <- as.list(read.csv(file.path(outdir, "canopy_params.csv"))[1, ])
layers <- read.csv(file.path(outdir, "canopy_layers.csv"))

# Near-zero-but-nonzero swdown/difrad right at sunrise/sunset.
solar_pre <- solarposition(p$lat, p$lon, year = climdata$obs_time)
lowsun <- solar_pre$zen >= 89
climdata$swdown[lowsun] <- 0
climdata$difrad[lowsun] <- 0

# ── Soil: createsoilc(soiltype=p$soiltype) -- "Clay loam" placeholder for
# soil_source=:slga (overridden below with real SLGA Smax/b/psi_e/Ksat/rho),
# or the real texture name otherwise, trusted verbatim (soil_override=false)
# so Vq/Vm/Vo/Smin/n stay self-consistent with Smax/b/psi_e/Ksat instead of
# mixing a placeholder's structure with a different texture's numbers (that
# mix caused a total NaN cascade). Flat/non-layered either way.
soilc <- createsoilc(soiltype = p$soiltype, nlayers = p$nlayers, totalDepth = p$totaldepth)
if (tolower(as.character(p$soil_override)) == "true") {
  soilc$Smax  <- p$smax
  soilc$b     <- p$b
  soilc$psi_e <- p$psi_e
  soilc$Ksat  <- p$ksat
  soilc$rho   <- p$rho_kgm3
}
soilc$gref  <- p$gref
soilc$groundem <- p$groundem

# Organic litter-layer override -- per-site toggle (SITE_ORGANIC_CAP in
# comparisons/ozflux/config.jl)
if (tolower(as.character(p$organic_cap)) == "true") {
  Vsolid_total <- soilc$Vq[1] + soilc$Vm[1] + soilc$Vo[1]
  vq_ratio <- soilc$Vq[1] / (soilc$Vq[1] + soilc$Vm[1])
  vm_ratio <- soilc$Vm[1] / (soilc$Vq[1] + soilc$Vm[1])
  set_organic_node <- function(i, f) {
    soilc$Vo[i] <<- f * Vsolid_total
    soilc$Vq[i] <<- (1 - f) * Vsolid_total * vq_ratio
    soilc$Vm[i] <<- (1 - f) * Vsolid_total * vm_ratio
  }
  set_organic_node(1, 0.98)
  set_organic_node(2, 0.98)
  set_organic_node(3, 0.98)
  set_organic_node(4, 0.50)
}

# ── Vegetation: createvegp("BET.Te") as a structural base -- Eller
# photosynthesis/hydraulic fields have no Microclimate.jl equivalent and
# stay at BET.Te's own defaults (see micropoint/README.md). ──
vegp <- createvegp(vegtype = "BET.Te")
vegp$h     <- cp$canopy_height_m
vegp$pai   <- cp$pai
vegp$x     <- cp$leaf_angle_x
vegp$lref  <- cp$leaf_reflectance
vegp$ltra  <- cp$leaf_transmittance
vegp$vegem <- cp$leaf_emissivity
vegp$len   <- cp$leaf_length_m
vegp$wid   <- cp$leaf_width_m
vegp$mwft  <- cp$leaf_water_storage_mm

# paii from canopy_layers.csv: micropoint's own evenly-spaced bottom-to-top
# grid, shaped to match Julia's PAI profile.
paii <- layers$paii
n <- length(paii)
Lfrac <- rep(0.75, n)  # no Julia equivalent; mid-range, avoids the psi50=1.0 edge case

zref <- p$zref  # true tower-top reference height, not a hardcoded 2m

# RunModelFull estimates mean annual temperature from a full year of
# climdata by default -- pass it explicitly (same fallback run_micropoint.R
# uses for RunMicro/InitailizeWater) so a truncated run doesn't fail.
boundaryT <- mean(climdata$temp)

# ── Chunked run: RunModelFull can hit a permanent, non-recovering NaN
# cascade under difficult real forcing. Soil state (weeks to months of memory)
# carries forward between chunks via SoilTempIni/ThetaIni, using the last
# non-NaN row.
#
# Chunked by fixed row count, NOT by calendar month parsed from obs_time to
# avoid a ~10-hour leading fragment whenever the UTC offset doesn't land on 
# a local midnight.
chunk_size <- 24 * 30  # ~1 month; exact boundary doesn't matter, just needs to be small enough to contain a cascade
chunk_id <- ceiling(seq_len(nrow(climdata)) / chunk_size)
chunks <- split(seq_len(nrow(climdata)), chunk_id)
chunk_names <- as.character(seq_along(chunks))

mout_list <- vector("list", length(chunks))
SoilTempIni <- NA
ThetaIni <- NA
t0 <- Sys.time()
for (i in seq_along(chunks)) {
  cd_chunk <- climdata[chunks[[i]], ]
  m <- RunModelFull(cd_chunk, soilc, vegp, paii, Lfrac, lat = p$lat, long = p$lon, zref = zref,
                     boundaryT = boundaryT, SoilTempIni = SoilTempIni, ThetaIni = ThetaIni)
  n_na <- sum(is.na(m$tair[, 1]))
  if (n_na > 0) {
    cat(sprintf("  chunk %s: %d/%d hours NaN (cascade)\n", chunk_names[i], n_na, nrow(cd_chunk)))
  }
  good <- which(!is.na(m$tsoil[, 1]))
  if (length(good) > 0 && length(good) / nrow(cd_chunk) >= 0.9) {
    last_good <- max(good)
    SoilTempIni <- as.numeric(m$tsoil[last_good, ])
    ThetaIni <- as.numeric(m$theta[last_good, ])
  } else {
    SoilTempIni <- NA
    ThetaIni <- NA
  }
  mout_list[[i]] <- m
}
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("micropoint RunModelFull (vegetated, %d layers, %d hours, %d monthly chunks): %.2f s\n",
            n, nrow(climdata), length(chunks), elapsed_s))

# Concatenate every chunk's matrices back into one continuous full-period result.
var_names <- names(mout_list[[1]])
mout <- setNames(lapply(var_names, function(nm) do.call(rbind, lapply(mout_list, `[[`, nm))), var_names)

# Rdirdown is beam-normal, not horizontal (confirmed via ?RunModelFull) --
# convert to horizontal so it's directly comparable to Microclimate.jl's
# own horizontal-irradiance convention.
coszen <- pmax(cos(solar_pre$zen * pi / 180), 0)
mout$Rdirdown <- mout$Rdirdown * coszen

for (nm in names(mout)) {
  write.csv(mout[[nm]], file.path(outdir, paste0("veg_", nm, ".csv")), row.names = FALSE)
}
write.csv(data.frame(elapsed_s = elapsed_s, n_layers = n, n_hours = nrow(climdata)),
          file.path(outdir, "veg_timing.csv"), row.names = FALSE)

# tsoil/theta column depths -- micropoint's own geometric node spacing
# (geometricCpp, the same internal spacing createsoilc uses), not a
# hand-derived approximation. Written out so Julia can match columns to
# real depths instead of guessing the spacing.
r_depth_m <- getFromNamespace("geometricCpp", "micropoint")(p$nlayers, p$totaldepth)
write.csv(data.frame(depth_m = r_depth_m), file.path(outdir, "soil_depths.csv"), row.names = FALSE)

cat(sprintf("Wrote %d vegetated output matrices + veg_timing.csv + soil_depths.csv to %s\n", length(mout), outdir))
