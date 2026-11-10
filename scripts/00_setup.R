
# 00_setup.R - one-time environment setup. Run from the project root.
# renv.lock records exact package versions; others rebuild with renv::restore().

# install.packages("renv"); renv::init()   # already done

# CRAN: dplyr (wrangling), PRROC (PR/ROC curves), ggplot2 (plots),
#       doParallel + doRNG (let GENIE3 use several cores reproducibly)
install.packages(c("BiocManager", "dplyr", "PRROC",
                   "ggplot2", "doParallel", "doRNG"))

# Bioconductor: GENIE3 (the method), GEOquery (downloads GEO data)
BiocManager::install(c("GENIE3", "GEOquery"), update = FALSE, ask = FALSE)

renv::snapshot()   # writes renv.lock = the environment file for submission

