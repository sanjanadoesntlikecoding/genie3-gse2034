
# 01_download_data.R - downloads all data (raw data is NOT committed to GitHub)
#   GSE2034 expression -> data/gse2034.rds
#   TRRUST v2 human    -> data/trrust_human.tsv  (gold standard)
# If GEO times out: download GSE2034_series_matrix.txt.gz and GPL96.annot.gz
# from GEO in a browser into data/raw/ and rerun this script.

source("scripts/config.R")
suppressPackageStartupMessages(library(GEOquery))
dir.create(CFG$raw_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1200)   # default 60 s is too short for the ~40 MB files

# 1. Expression data + probe annotation
gse_rds <- file.path(CFG$data_dir, "gse2034.rds")
if (!file.exists(gse_rds)) {
  # GSEMatrix = processed series matrix; getGPL = probe -> gene annotation
  gse <- getGEO(CFG$gse_id, GSEMatrix = TRUE, getGPL = TRUE,
                destdir = CFG$raw_dir)[[1]]
  saveRDS(gse, gse_rds)
}
gse <- readRDS(gse_rds)
cat("Probes x samples:", dim(Biobase::exprs(gse)), "\n")   # expect 22283 x 286

# 2. TRRUST gold standard
trrust_file <- file.path(CFG$data_dir, "trrust_human.tsv")
if (!file.exists(trrust_file)) download.file(CFG$trrust_url, trrust_file, mode = "wb")
tr <- read.delim(trrust_file, header = FALSE)
cat("TRRUST rows:", nrow(tr), "| unique TF-target pairs:", nrow(unique(tr[, 1:2])), "\n")
# 9396 rows but ~8427 unique pairs: some links are listed once per PubMed ID

