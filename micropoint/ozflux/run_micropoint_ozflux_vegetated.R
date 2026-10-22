# run_micropoint_ozflux_vegetated.R — runs micropoint's own vegetated
# RunModelFull on OzFlux tower forcing (climdata.csv/params.csv/
# soil_profile.csv, already written by write_ozflux_micropoint_inputs.jl
# from the same comparisons/ozflux `result` Microclimate.jl solved) plus
# canopy_params.csv (that run's own MultilayerCanopy spec) and
# canopy_layers.csv (micropoint's own evenly-spaced PAI grid, shaped to
# match Julia's PAI profile). Adapted from micropoint/run_micropoint_
# vegetated.R (same vegp mapping) and micropoint/run_micropoint.R
# (organic-top-layer soil override), but reads zref/PAI shape from files
# instead of hardcoding.
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

lowsun <- climdata$zenith_deg >= 89
climdata$swdown[lowsun] <- 0
climdata$difrad[lowsun] <- 0

# Soil node depths as RunModelFull builds them (node nlayers at totalDepth)
r_depth_m <- head(getFromNamespace("geometricCpp", "micropoint")(p$nlayers, p$totaldepth), p$nlayers + 1)

# ── Soil: createsoilc(soiltype=p$soiltype) -- "Clay loam" placeholder for
# soil_source=:slga (overridden below with soil_profile.csv's real SLGA
# Smax/b/psi_e/Ksat and porosity), or the real texture name otherwise, trusted
# verbatim (soil_override=false) so Vq/Vm/Vo/Smin/n stay self-consistent
# with Smax/b/psi_e/Ksat instead of mixing a placeholder's structure with a
# different texture's numbers (that mix caused a total NaN cascade).
soilc <- createsoilc(soiltype = p$soiltype, nlayers = p$nlayers, totalDepth = p$totaldepth)
if (tolower(as.character(p$soil_override)) == "true") {
  # soil_profile.csv is Julia's own depth grid -- interpolate onto r_depth_m
  # so every vector handed to RunModelFull is the correct nlayers length
  # regardless of how the two depth grids differ (a length-1 or mismatched
  # vector here is undefined behavior in RunModelFull's per-layer indexing,
  # not a real numerical failure).
  soil_profile <- read.csv(file.path(outdir, "soil_profile.csv"))
  interp <- function(col) approx(soil_profile$depth_m, soil_profile[[col]], r_depth_m, rule = 2)$y
  soilc$Smax  <- interp("smax")
  soilc$b     <- interp("b")
  soilc$psi_e <- interp("psi_e")
  soilc$Ksat  <- interp("ksat")
  # createsoilc has no bulk density: porosity is 1 - (Vq + Vm + Vo), so rescale
  # the solid fractions to Julia's porosity, keeping their proportions
  solids <- soilc$Vq + soilc$Vm + soilc$Vo
  scale <- (1 - soilc$Smax) / solids
  soilc$Vq <- soilc$Vq * scale
  soilc$Vm <- soilc$Vm * scale
  soilc$Vo <- soilc$Vo * scale
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
# micropoint needs x > 0; Microclimate.jl's x = 0 (vertical leaves) becomes nearly vertical
vegp$x     <- max(cp$leaf_angle_x, 0.01)
if (cp$leaf_angle_x < 0.01) cat(sprintf("  leaf angle x = %g raised to 0.01 for micropoint\n", cp$leaf_angle_x))
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

# Julia's own initial soil state, interpolated onto r_depth_m -- seeds the
# run instead of RunModelFull's generic auto-generated guess.
initial_state <- read.csv(file.path(outdir, "initial_soil_state.csv"))
SoilTempIni <- approx(initial_state$depth_m, initial_state$temp_C, r_depth_m, rule = 2)$y
ThetaIni <- approx(initial_state$depth_m, initial_state$moisture_frac, r_depth_m, rule = 2)$y

# One continuous run over the whole period
t0 <- Sys.time()
mout <- RunModelFull(climdata, soilc, vegp, paii, Lfrac, lat = p$lat, long = p$lon, zref = zref,
                     boundaryT = boundaryT, SoilTempIni = SoilTempIni, ThetaIni = ThetaIni)
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
n_na <- sum(is.na(mout$tair[, 1]))
if (n_na > 0) {
  cat(sprintf("  %d/%d hours NaN\n", n_na, nrow(climdata)))
}
cat(sprintf("micropoint RunModelFull (vegetated, %d layers, %d hours): %.2f s\n",
            n, nrow(climdata), elapsed_s))

# Rdirdown is beam-normal, not horizontal (confirmed via ?RunModelFull) --
# convert to horizontal so it's directly comparable to Microclimate.jl's
# own horizontal-irradiance convention.
coszen <- pmax(cos(climdata$zenith_deg * pi / 180), 0)
mout$Rdirdown <- mout$Rdirdown * coszen

for (nm in names(mout)) {
  write.csv(mout[[nm]], file.path(outdir, paste0("veg_", nm, ".csv")), row.names = FALSE)
}
write.csv(data.frame(elapsed_s = elapsed_s, n_layers = n, n_hours = nrow(climdata)),
          file.path(outdir, "veg_timing.csv"), row.names = FALSE)

# tsoil/theta column depths, as RunModelFull reports them. Written out so
# Julia can match columns to real depths instead of guessing the spacing.
if (!isTRUE(all.equal(as.numeric(mout$zsoil), r_depth_m, tolerance = 1e-6))) {
  warning("RunModelFull's soil depths differ from r_depth_m used to interpolate the inputs")
}
write.csv(data.frame(depth_m = as.numeric(mout$zsoil)), file.path(outdir, "soil_depths.csv"), row.names = FALSE)

cat(sprintf("Wrote %d vegetated output matrices + veg_timing.csv + soil_depths.csv to %s\n", length(mout), outdir))
