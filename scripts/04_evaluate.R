
# 04_evaluate.R - GENIE3 vs correlation network vs random, scored on TRRUST

source("scripts/config.R"); source("scripts/utils.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

d <- readRDS(file.path(CFG$data_dir, "processed.rds"))
W <- readRDS(file.path(CFG$results_dir, "genie3_weights.rds"))
W <- W[d$regulators, rownames(d$expr)]

set.seed(CFG$seed)
methods <- list(
  GENIE3_RF    = W,
  Pearson_abs  = cor_scores(d$expr, d$regulators, "pearson"),
  Spearman_abs = cor_scores(d$expr, d$regulators, "spearman"),
  Random       = matrix(runif(length(W)), nrow(W), dimnames = dimnames(W))
)

# 1. metrics for each method
res <- lapply(names(methods), function(m) {
  e <- evaluate_scores(methods[[m]], d$gold, curve = TRUE)
  list(tab = cbind(method = m, e), curve = data.frame(method = m, attr(e, "curve")))
})
tab <- bind_rows(lapply(res, function(r) r$tab))
print(tab, digits = 3)
write.csv(tab, file.path(CFG$results_dir, "eval_methods.csv"), row.names = FALSE)

# 2. precision-recall curves
curves <- bind_rows(lapply(res, function(r) r$curve))
colnames(curves)[2:3] <- c("recall", "precision")
curves <- curves %>% group_by(method) %>%
  slice(unique(round(seq(1, n(), length.out = 2000))))
p <- ggplot(curves, aes(recall, precision, colour = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = tab$base_rate[1], linetype = 2) +
  labs(title = "GENIE3 vs co-expression on GSE2034 (gold standard: TRRUST)",
       subtitle = "Dashed = random guessing") + theme_bw()
ggsave(file.path(CFG$results_dir, "pr_curves.png"), p, width = 7, height = 5)

# 3. sign check: GENIE3 weights are always >= 0, so they cannot say
#    activation vs repression. Does the correlation sign match TRRUST?
rho <- cor(t(d$expr[d$regulators, ]), t(d$expr), method = "spearman")
sign_check <- d$gold %>%
  filter(mode %in% c("Activation", "Repression")) %>%
  mutate(rho = rho[cbind(TF, target)],
         genie3_weight = W[cbind(TF, target)],
         agrees = ifelse(mode == "Activation", rho > 0, rho < 0))
print(sign_check %>% group_by(mode) %>%
        summarise(n = n(), frac_agree = mean(agrees), median_rho = median(rho)))
write.csv(sign_check, file.path(CFG$results_dir, "sign_check.csv"), row.names = FALSE)

# 4. hub TFs among the top 1000 GENIE3 edges
top <- score_table(W) %>% arrange(desc(score)) %>% head(1000) %>%
  mutate(in_trrust = paste(regulator, target) %in% paste(d$gold$TF, d$gold$target))
hubs <- top %>% group_by(regulator) %>%
  summarise(n_edges = n(), n_confirmed = sum(in_trrust),
            total_weight = sum(score), .groups = "drop") %>%
  arrange(desc(n_edges))
print(head(hubs, 15))
write.csv(hubs, file.path(CFG$results_dir, "top_hubs.csv"), row.names = FALSE)
write.csv(top, file.path(CFG$results_dir, "genie3_top1000_annotated.csv"), row.names = FALSE)

