source("scripts/config.R")
suppressPackageStartupMessages({ library(Biobase); library(dplyr) })

gse <- readRDS(file.path(CFG$data_dir, "gse2034.rds"))
X   <- exprs(gse)
fd  <- fData(gse)

# 1. log2 transform if needed
if (quantile(X, 0.99, na.rm = TRUE) > 100) {
  X[X < 1] <- 1
  X <- log2(X)
  message("Applied log2 transform.")
}

# 2. probes -> gene symbols
sym <- trimws(as.character(fd[["Gene Symbol"]]))
keep <- !is.na(sym) & sym != "" & !grepl("///", sym) & !grepl("^AFFX", rownames(X))
X <- X[keep, ]; sym <- sym[keep]

# several probes per gene -> keep the highest-expressed one
ord <- order(sym, -rowMeans(X))
X <- X[ord, ]; sym <- sym[ord]
first <- !duplicated(sym)
G <- X[first, ]; rownames(G) <- sym[first]
cat("Genes after collapsing:", nrow(G), "\n")

# 3. keep only genes in TRRUST
tr <- read.delim(file.path(CFG$data_dir, "trrust_human.tsv"), header = FALSE,
                 col.names = c("TF", "target", "mode", "pmid"), stringsAsFactors = FALSE)
G <- G[rownames(G) %in% union(tr$TF, tr$target), ]

# 4. drop lowest-expressed 25%, keep the 1000 most variable
G <- G[rowMeans(G) > quantile(rowMeans(G), 0.25), ]
v <- apply(G, 1, var)
G <- G[names(sort(v, decreasing = TRUE))[seq_len(min(CFG$n_genes, length(v)))], ]

# 5. scale each gene to mean 0, sd 1
expr <- t(scale(t(G)))

# 6. regulators + gold standard
regulators <- intersect(unique(tr$TF), rownames(expr))
gold <- tr %>%
  filter(TF %in% regulators, target %in% rownames(expr), TF != target) %>%
  group_by(TF, target) %>%
  summarise(mode = if (n_distinct(mode) == 1) first(mode) else "Mixed",
            n_pmids = n(), .groups = "drop")

cat("Matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Regulators:", length(regulators), "| gold edges:", nrow(gold), "\n")
print(table(gold$mode))

saveRDS(list(expr = expr, regulators = regulators, gold = gold, pheno = pData(gse)),
        file.path(CFG$data_dir, "processed.rds"))
