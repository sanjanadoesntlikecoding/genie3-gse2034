if (!file.exists("renv.lock")) {
  stop("Run this from the project root (open 'systems bio.Rproj' in RStudio).")
}

CFG <- list(
  gse_id      = "GSE2034",
  trrust_url  = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
  data_dir    = "data",
  raw_dir     = "data/raw",
  results_dir = "results",
  n_genes     = 1000,
  tree_method = "RF",
  K           = "sqrt",
  n_trees     = 500,
  seed        = 123,
  n_cores     = max(1, parallel::detectCores() - 1)
)

dir.create(CFG$results_dir, showWarnings = FALSE, recursive = TRUE)