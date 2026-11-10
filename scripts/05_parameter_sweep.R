
# 05_parameter_sweep.R - effect of K, nTrees, RF vs ET, and sample size.
# Uses the same 300 targets for every run so results are comparable.
# Each run is cached in results/sweep/, so it is safe to stop and restart.

source("scripts/config.R"); source("scripts/utils.R")
suppressPackageStartupMessages({ library(GENIE3); library(dplyr); library(ggplot2) })

d <- readRDS(file.path(CFG$data_dir, "processed.rds"))
sweep_dir <- file.path(CFG$results_dir, "sweep")
dir.create(sweep_dir, showWarnings = FALSE)

set.seed(CFG$seed)
targets  <- sample(unique(d$gold$target), min(300, length(unique(d$gold$target))))
gold_sub <- d$gold %>% filter(target %in% targets)

grid <- bind_rows(
  data.frame(experiment = "K", K = c("1", "sqrt", "50", "all"), nTrees = 500, method = "RF", n = 286, rep = 1),
  data.frame(experiment = "nTrees", K = "sqrt", nTrees = c(10, 50, 250, 1000), method = "RF", n = 286, rep = 1),
  data.frame(experiment = "method", K = "sqrt", nTrees = 500, method = c("RF", "ET"), n = 286, rep = 1),
  expand.grid(experiment = "samples", K = "sqrt", nTrees = 500, method = "RF",
              n = c(30, 60, 120, 200, 286), rep = 1:3, stringsAsFactors = FALSE)
)
grid$id <- with(grid, paste(experiment, K, nTrees, method, n, rep, sep = "_"))

run_one <- function(g) {
  f <- file.path(sweep_dir, paste0(g$id, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  set.seed(CFG$seed + g$rep)
  cols <- if (g$n < ncol(d$expr)) sample(ncol(d$expr), g$n) else seq_len(ncol(d$expr))
  e <- t(scale(t(d$expr[, cols])))
  K <- if (grepl("^[0-9]+$", g$K)) as.numeric(g$K) else g$K
  t <- system.time(
    W <- GENIE3(e, regulators = d$regulators, targets = targets,
                treeMethod = g$method, K = K, nTrees = g$nTrees, nCores = CFG$n_cores)
  )
  W <- W[d$regulators, targets]
  pear <- cor_scores(e, d$regulators, "pearson")[, targets]
  out <- bind_rows(
    cbind(g, scorer = "GENIE3",  evaluate_scores(W, gold_sub), minutes = unname(t["elapsed"]) / 60),
    cbind(g, scorer = "Pearson", evaluate_scores(pear, gold_sub), minutes = NA)
  )
  saveRDS(out, f)
  message(g$id, ": AUPR = ", signif(out$AUPR[1], 3), " (", round(t["elapsed"] / 60, 1), " min)")
  out
}

sweep <- bind_rows(lapply(split(grid, seq_len(nrow(grid))), run_one))
write.csv(sweep, file.path(CFG$results_dir, "sweep_results.csv"), row.names = FALSE)

# plot 1: K / nTrees / method
pear_ref <- sweep %>% filter(scorer == "Pearson", n == 286) %>% pull(AUPR) %>% mean()
p1 <- sweep %>% filter(scorer == "GENIE3", experiment != "samples") %>%
  mutate(setting = ifelse(experiment == "K", K,
                   ifelse(experiment == "nTrees", as.character(nTrees), method))) %>%
  ggplot(aes(factor(setting, levels = unique(setting)), AUPR)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = pear_ref, linetype = 2, colour = "red") +
  facet_wrap(~experiment, scales = "free_x") +
  labs(x = NULL, subtitle = "Red dashed = |Pearson| baseline") + theme_bw()
ggsave(file.path(CFG$results_dir, "sweep_params.png"), p1, width = 9, height = 4)

# plot 2: sample size, GENIE3 vs Pearson
p2 <- sweep %>% filter(experiment == "samples") %>%
  group_by(scorer, n) %>%
  summarise(m = mean(AUPR), lo = min(AUPR), hi = max(AUPR), .groups = "drop") %>%
  ggplot(aes(n, m, colour = scorer)) +
  geom_line() + geom_pointrange(aes(ymin = lo, ymax = hi)) +
  labs(x = "Number of tumours", y = "AUPR") + theme_bw()
ggsave(file.path(CFG$results_dir, "sweep_samples.png"), p2, width = 6, height = 4)

print(sweep %>% filter(scorer == "GENIE3") %>%
        select(experiment, K, nTrees, method, n, rep, AUPR, prec_100, minutes), digits = 3)



# plot 1 (redrawn): settings in logical order
p1 <- sweep %>% filter(scorer == "GENIE3", experiment != "samples") %>%
  mutate(setting = case_when(experiment == "K" ~ K,
                             experiment == "nTrees" ~ as.character(nTrees),
                             TRUE ~ method),
         setting = factor(setting, levels = c("1","sqrt","50","all","10","250","500","1000","RF","ET"))) %>%
  ggplot(aes(setting, AUPR)) + geom_col(fill = "steelblue") +
  geom_hline(yintercept = pear_ref, linetype = 2, colour = "red") +
  facet_wrap(~experiment, scales = "free_x") +
  labs(x = NULL, subtitle = "Red dashed = |Pearson| baseline (same 300 targets)") + theme_bw()
ggsave(file.path(CFG$results_dir, "sweep_params.png"), p1, width = 9, height = 4)


# plot 1 (redrawn again): distinct labels so nTrees sorts correctly
p1 <- sweep %>% filter(scorer == "GENIE3", experiment != "samples") %>%
  mutate(setting = case_when(experiment == "K" ~ paste("K =", K),
                             experiment == "nTrees" ~ paste(nTrees, "trees"),
                             TRUE ~ method),
         setting = factor(setting, levels = c("K = 1","K = sqrt","K = 50","K = all",
                                              "RF","ET",
                                              "10 trees","50 trees","250 trees","1000 trees"))) %>%
  ggplot(aes(setting, AUPR)) + geom_col(fill = "steelblue") +
  geom_hline(yintercept = pear_ref, linetype = 2, colour = "red") +
  facet_wrap(~experiment, scales = "free_x") +
  labs(x = NULL, subtitle = "Red dashed = |Pearson| baseline (same 300 targets)") + theme_bw()
ggsave(file.path(CFG$results_dir, "sweep_params.png"), p1, width = 9, height = 4)
