
# 03_run_genie3.R - main GENIE3 run
# For each target gene, a Random Forest predicts its expression across the
# 286 tumours from the 223 TFs. Each TF importance = weight of edge TF -> target.

source("scripts/config.R")
suppressPackageStartupMessages({ library(GENIE3); library(doParallel); library(doRNG) })  # doParallel + doRNG: multi-core, reproducible
d <- readRDS(file.path(CFG$data_dir, "processed.rds"))

set.seed(CFG$seed)
t <- system.time(
  W <- GENIE3(d$expr, regulators = d$regulators, treeMethod = CFG$tree_method,
              K = CFG$K, nTrees = CFG$n_trees, nCores = CFG$n_cores, verbose = TRUE)
)
cat("Done in", round(t["elapsed"] / 60, 1), "min\n")

saveRDS(W, file.path(CFG$results_dir, "genie3_weights.rds"))
links <- getLinkList(W, reportMax = 1000)
write.csv(links, file.path(CFG$results_dir, "genie3_top1000.csv"), row.names = FALSE)
print(head(links, 20))

