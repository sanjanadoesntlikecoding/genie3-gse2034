
# utils.R - helpers for scoring a predicted network against TRRUST

# regulators x targets matrix -> long edge table (as.vector reads column-wise)
score_table <- function(S) {
  df <- data.frame(regulator = rep(rownames(S), times = ncol(S)),
                   target    = rep(colnames(S), each  = nrow(S)),
                   score     = as.vector(S), stringsAsFactors = FALSE)
  df[df$regulator != df$target, ]
}

# Pairs not in TRRUST count as negatives. TRRUST is incomplete, so some
# "false positives" may be real but undocumented -> AUPR is a lower bound.
evaluate_scores <- function(S, gold, curve = FALSE) {
  df <- score_table(S)
  df$label <- paste(df$regulator, df$target) %in% paste(gold$TF, gold$target)
  fg <- df$score[df$label]; bg <- df$score[!df$label]
  pr  <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = curve)
  roc <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg)
  ranked <- df$label[order(df$score, decreasing = TRUE)]
  out <- data.frame(AUPR = pr$auc.integral, AUROC = roc$auc,
                    base_rate = mean(df$label),
                    prec_100 = mean(ranked[1:100]),
                    prec_1000 = mean(ranked[1:1000]),
                    n_pos = sum(df$label), n_pairs = nrow(df))
  out$AUPR_fold <- out$AUPR / out$base_rate
  if (curve) attr(out, "curve") <- pr$curve
  out
}

# Co-expression baseline (course method): |correlation| regulator vs gene
cor_scores <- function(expr, regulators, method = "pearson") {
  abs(cor(t(expr[regulators, , drop = FALSE]), t(expr), method = method))
}

