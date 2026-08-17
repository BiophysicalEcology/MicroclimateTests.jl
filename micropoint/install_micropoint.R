# install_micropoint.R — one-time setup. Run once with:
#   Rscript install_micropoint.R
#
# micropoint (https://github.com/ilyamaclean/micropoint) isn't on CRAN.
# Few dependencies beyond a recent Rcpp.

cran_deps <- c("remotes", "Rcpp")
installed <- rownames(installed.packages())
to_install <- setdiff(cran_deps, installed)
if (length(to_install) > 0) {
    install.packages(to_install)
}

#if (!("micropoint" %in% installed.packages()[, "Package"])) {
    remotes::install_github("ilyamaclean/micropoint")
#} else {
#    cat("micropoint already installed.\n")
#}

cat("\nDone. Run the comparison with:\n")
cat("  cd c:/git/BiophysicalEcologyEnv && julia --project=. c:/git/MicroclimateTests.jl/micropoint/comparison.jl\n")
