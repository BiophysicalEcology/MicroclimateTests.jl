# run_micropoint_vignette.R — reproduces micropoint's own published vignette
# forest example (createvegp("BET.Te"), PAIgeometry, createsoilc("Clay loam"),
# RunModelFull) and writes every input/output as CSV for
# vignette_comparison.jl to read and reproduce with Microclimate.jl.
#
# Unlike the sibling micropoint/ comparison (real SNOTEL sites, bare ground,
# Julia writes inputs / R reads them), this one runs the OTHER direction: R's
# own canned climdata is the source of truth, so R writes climdata.csv and
# Julia reads it back to drive an identical-forcing MultilayerCanopy run. See
# README.md for the full parameter-mapping/known-mismatch writeup.
#
# Usage: Rscript run_micropoint_vignette.R <outdir>

options(error = function() {
  traceback(2)
  quit(status = 1)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
if (is.na(outdir)) stop("Usage: Rscript run_micropoint_vignette.R <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(micropoint))

lat  <- 49.96807
long <- -5.215668

# ── Vegetation: forest example (BET.Te — Broadleaf Evergreen Tree, Temperate) ─
vegp <- createvegp(vegtype = "BET.Te")
n <- 20
paii  <- PAIgeometry(vegp$pai, 7, 70, n = n)   # bottom-to-top, sums to vegp$pai
Lfrac <- seq(0.5, 0.85, length.out = length(paii))

# ── Soil: Clay loam (vignette default texture) ────────────────────────────────
soilc <- createsoilc(soiltype = "Clay loam")

# ── Forcing: micropoint's own built-in climdata, height-adjusted to above
# canopy (h=25m) exactly as the vignette's forest example does ───────────────
climdata_adj <- weatherhgt_adjust(climdata, zin = 2, zout = 27, lat = lat, long = long)

t0 <- Sys.time()
mout <- RunModelFull(climdata_adj, soilc, vegp, paii, Lfrac, lat = lat, long = long, zref = 27)
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("micropoint RunModelFull (forest, %d layers, %d hours): %.2f s\n",
            n, nrow(climdata_adj), elapsed_s))

# ── Write forcing (so Julia reads the exact series micropoint itself used) ───
write.csv(climdata_adj, file.path(outdir, "climdata.csv"), row.names = FALSE)

# ── Write vegetation/soil parameters (every mappable field; unmapped fields
# are still included for completeness/documentation, ignored by Julia side) ──
write.csv(as.data.frame(vegp), file.path(outdir, "vegp.csv"), row.names = FALSE)
write.csv(data.frame(paii = paii, Lfrac = Lfrac), file.path(outdir, "paii.csv"), row.names = FALSE)
write.csv(as.data.frame(soilc[c("Smax","Smin","n","Ksat","b","psi_e","rho")]),
          file.path(outdir, "soilc.csv"), row.names = FALSE)
write.csv(data.frame(gref = soilc$gref, groundem = soilc$groundem, grefPAR = soilc$grefPAR,
                      nLayers = soilc$nLayers, totalDepth = soilc$totalDepth,
                      lat = lat, long = long, zref = 27, canopy_height = vegp$h),
          file.path(outdir, "params.csv"), row.names = FALSE)

# Exact per-node soil depths (m) micropoint solves at -- geometrically spaced
# (millimetres near the surface, tens of cm at depth), not uniform. Export
# rather than approximate on the Julia side: the diurnal thermal wave needs
# fine near-surface resolution, and a uniform grid badly damps/lags it.
soil_depths <- micropoint:::geometricCpp(soilc$nLayers, soilc$totalDepth)[1:(soilc$nLayers + 1)]
write.csv(data.frame(depth_m = soil_depths), file.path(outdir, "soil_depths.csv"), row.names = FALSE)

# ── Write output matrices ─────────────────────────────────────────────────────
for (nm in names(mout)) {
  write.csv(mout[[nm]], file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
}

write.csv(data.frame(elapsed_s = elapsed_s, n_layers = n, n_hours = nrow(climdata_adj)),
          file.path(outdir, "timing.csv"), row.names = FALSE)

cat(sprintf("Wrote climdata/vegp/paii/soilc/params + %d output matrices + timing to %s\n",
            length(mout), outdir))
