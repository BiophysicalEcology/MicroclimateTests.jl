# single_site.R — interactive single-site NMR run.
#
# comparison.jl launches run_nmr.R as a background Rscript subprocess (up to
# 16 concurrent), so a failing run's real R error/traceback gets lost --
# only the exit code survives. This sources run_nmr.R directly against a
# site's already-written CSVs (nmr_outputs/<site>/*.csv, from Julia's
# write_nmr_inputs) so RStudio shows the full error and leaves every
# intermediate object (micro, microinput, forcing, ...) in the environment
# for inspection.
#
# Usage: set `site` below, then run/source this file -- works from any
# working directory (paths below are absolute, not relative to getwd()).

site <- "a1"

oznet_dir <- "c:/git/MicroclimateTests.jl/comparisons/oznet"
outdir <- file.path(oznet_dir, "nmr_outputs", site)
if (!dir.exists(outdir)) stop("no nmr_outputs for site ", site, " -- run comparison.jl/single_site.jl first")

source(file.path(oznet_dir, "run_nmr.R"))

