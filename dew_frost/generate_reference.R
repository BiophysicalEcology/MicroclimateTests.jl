# generate_reference.R
#
# Runs NicheMapR's micro_global for Kimba, SA (cap=0, all other args default)
# and dumps the pieces needed to drive an equivalent Microclimate.jl run with
# byte-identical hourly forcing -- see compare_kimba.jl. Regenerate with:
#   Rscript generate_reference.R

library(NicheMapR)

out_dir <- file.path(dirname(sys.frame(1)$ofile %||% "."), "data")
if (!dir.exists(out_dir)) out_dir <- "data"  # fallback when sourced interactively
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

micro <- micro_global(loc = c(136.4, -33.1), cap = 0)

write.csv(micro$metout, file.path(out_dir, "metout.csv"), row.names = FALSE)
write.csv(micro$soil, file.path(out_dir, "soil.csv"), row.names = FALSE)
write.csv(micro$soilmoist, file.path(out_dir, "soilmoist.csv"), row.names = FALSE)
write.csv(data.frame(DEP = micro$DEP), file.path(out_dir, "DEP.csv"), row.names = FALSE)
write.csv(data.frame(PE = micro$PE, BD = micro$BD, DD = micro$DD, BB = micro$BB, KS = micro$KS),
          file.path(out_dir, "soilprops.csv"), row.names = FALSE)
write.csv(data.frame(RAINFALL = micro$RAINFALL), file.path(out_dir, "RAINFALL.csv"), row.names = FALSE)
# RUF/Usrhyt/Refhyt/SLE are micro_global's own defaults (not returned in the
# list), reproduced here explicitly since this call didn't override them.
write.csv(data.frame(longitude = micro$longlat[1], latitude = micro$longlat[2],
                      elev = micro$elev, REFL = micro$REFL,
                      RUF = 0.004, Usrhyt = 0.01, Refhyt = 1.2, SLE = 0.95,
                      slope = 0, aspect = 0, timeinterval = micro$timeinterval,
                      ndays = micro$ndays),
          file.path(out_dir, "site_params.csv"), row.names = FALSE)
write.csv(data.frame(hori = rep(0, 24)), file.path(out_dir, "hori.csv"), row.names = FALSE)

cat("Wrote reference CSVs to", out_dir, "\n")
cat("n dew events (DEW>0):", sum(micro$metout[,"DEW"] > 0), "\n")
cat("n frost events (FROST>0):", sum(micro$metout[,"FROST"] > 0), "\n")
