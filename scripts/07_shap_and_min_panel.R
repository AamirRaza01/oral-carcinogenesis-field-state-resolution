# ==============================================================================
# Script: 07_shap_and_min_panel.R
# Purpose: Rank stable panel features by mean absolute SHAP values (XGBoost)
#          and perform forward selection to identify the minimal parsimonious
#          biomarker panel within AUC_TOLERANCE of the full signature.
# Reference: Section III-E, Section IV-D (Three-State Classification)
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

library(xgboost); library(pROC); library(ggplot2)

N_ITERS       <- 200    # Resampled iterations per panel size
AUC_TOLERANCE <- 0.02   # Maximum allowable AUC reduction vs. full panel
BALACC_FLOOR  <- 0.90   # Minimum balanced accuracy threshold
MIN_PANEL     <- 2

manifest <- load_manifest()
panels <- readRDS(file.path(res_dir, "stability_panels.rds"))
grid <- run_grid()
n_cores <- max(1, parallel::detectCores() - 1)

xgb_params <- function(nthread = n_cores)
  list(objective = "binary:logistic", eval_metric = "logloss", max_depth = 6,
       eta = 0.1, subsample = 0.8, colsample_bytree = 0.8, nthread = nthread)

evaluate_panel <- function(X, y, genes, n_iters = N_ITERS) {
  Xs <- X[, genes, drop = FALSE]; yb <- as.integer(y == "cancer")
  auc <- bacc <- rep(NA_real_, n_iters)
  for (it in seq_len(n_iters)) {
    sp <- stratified_split(y)
    dtr <- xgb.DMatrix(Xs[sp$train, , drop = FALSE], label = yb[sp$train])
    dvl <- xgb.DMatrix(Xs[sp$val,   , drop = FALSE], label = yb[sp$val])
    fit <- tryCatch(xgb.train(xgb_params(), dtr, nrounds = 500,
                              evals = list(val = dvl), early_stopping_rounds = 20,
                              verbose = 0), error = function(e) NULL)
    if (is.null(fit)) next
    p_te <- predict(fit, xgb.DMatrix(Xs[sp$test, , drop = FALSE]))
    p_vl <- predict(fit, dvl)
    r <- tryCatch(roc(yb[sp$val], p_vl, quiet = TRUE), error = function(e) NULL)
    th <- if (is.null(r)) 0.5 else r$thresholds[which.max(r$sensitivities + r$specificities - 1)]
    if (!is.finite(th)) th <- 0.5
    cls <- as.integer(p_te > th); yt <- yb[sp$test]
    se <- if (sum(yt == 1)) mean(cls[yt == 1] == 1) else NA_real_
    sn <- if (sum(yt == 0)) mean(cls[yt == 0] == 0) else NA_real_
    bacc[it] <- mean(c(se, sn), na.rm = TRUE)
    auc[it]  <- tryCatch(as.numeric(roc(yt, p_te, quiet = TRUE)$auc), error = function(e) NA_real_)
  }
  data.frame(auc = auc, bacc = bacc)
}

# ------------------------------------------------------------------------------
# SHAP ranking and forward panel selection
# ------------------------------------------------------------------------------
rank_rows <- list(); fwd_rows <- list(); min_rows <- list()

for (i in seq_len(nrow(grid))) {
  arm <- grid$arm[i]; set <- grid$set[i]; run <- grid$run[i]
  panel <- panels[[run]]
  if (is.null(panel) || length(panel) < MIN_PANEL) next

  expr <- load_arm(arm)
  samp <- intersect(set_samples(manifest, set), colnames(expr))
  genes <- intersect(panel, rownames(expr))
  X <- t(expr[genes, samp, drop = FALSE]); y <- set_labels(manifest, samp)

  # Tree SHAP feature ranking
  sp <- stratified_split(y); yb <- as.integer(y == "cancer")
  dtr <- xgb.DMatrix(X[sp$train, , drop = FALSE], label = yb[sp$train])
  dvl <- xgb.DMatrix(X[sp$val,   , drop = FALSE], label = yb[sp$val])
  fit <- xgb.train(xgb_params(), dtr, nrounds = 500, evals = list(val = dvl),
                   early_stopping_rounds = 20, verbose = 0)

  X_test <- X[sp$test, , drop = FALSE]
  shap_raw <- predict(fit, xgb.DMatrix(X_test), predcontrib = TRUE)
  shap <- shap_raw[, seq_along(genes), drop = FALSE]
  colnames(shap) <- genes
  mean_abs <- colMeans(abs(shap))
  ord <- order(mean_abs, decreasing = TRUE)
  gene_order <- genes[ord]

  rank_rows[[run]] <- data.frame(run = run, arm = arm, set = set,
                                 gene = gene_order, mean_abs_shap = mean_abs[ord],
                                 rank = seq_along(ord))

  # Forward selection along SHAP ranking
  fwd <- vector("list", length(gene_order))
  for (k in seq_along(gene_order)) {
    m <- evaluate_panel(X, y, gene_order[1:k])
    fwd[[k]] <- data.frame(
      run = run, arm = arm, set = set, panel_size = k,
      genes = paste(gene_order[1:k], collapse = ";"),
      last_added = gene_order[k],
      median_auc = median(m$auc, na.rm = TRUE),
      auc_lo = unname(quantile(m$auc, 0.025, na.rm = TRUE)),
      auc_hi = unname(quantile(m$auc, 0.975, na.rm = TRUE)),
      median_bacc = median(m$bacc, na.rm = TRUE))
  }
  fwd_df <- do.call(rbind, fwd)
  fwd_rows[[run]] <- fwd_df

  full_auc <- fwd_df$median_auc[nrow(fwd_df)]
  cand <- fwd_df[fwd_df$panel_size >= MIN_PANEL &
                 fwd_df$median_auc >= (full_auc - AUC_TOLERANCE) &
                 fwd_df$median_bacc >= BALACC_FLOOR, ]
  best <- if (nrow(cand)) cand[which.min(cand$panel_size), ] else fwd_df[nrow(fwd_df), ]
  best$full_panel_auc <- full_auc; best$delta_auc <- full_auc - best$median_auc
  best$met_criteria <- nrow(cand) > 0
  min_rows[[run]] <- best

  p <- ggplot(fwd_df, aes(panel_size, median_auc)) +
    geom_ribbon(aes(ymin = auc_lo, ymax = auc_hi), fill = "#3B4CC0", alpha = 0.15) +
    geom_line(color = "#3B4CC0", linewidth = 0.8) +
    geom_point(color = "#3B4CC0", size = 2) +
    geom_hline(yintercept = full_auc - AUC_TOLERANCE, linetype = "dashed", color = "red") +
    geom_vline(xintercept = best$panel_size, linetype = "dotted", color = "darkgreen") +
    scale_x_continuous(breaks = seq_len(nrow(fwd_df))) +
    labs(title = paste0("Forward selection, ", run),
         subtitle = sprintf("%s, set %s | full panel AUC %.4f", ARM_LABEL[[arm]], set, full_auc),
         x = "Panel size (genes)", y = "Median test AUC") +
    theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "plain"))
  ggsave(file.path(fig_dir, sprintf("elbow_%s.pdf", run)), p, width = 8, height = 5)
  ggsave(file.path(fig_dir, sprintf("elbow_%s.png", run)), p, width = 8, height = 5, dpi = 300)
}

write.csv(do.call(rbind, rank_rows), file.path(res_dir, "shap_rankings.csv"), row.names = FALSE)
write.csv(do.call(rbind, fwd_rows),  file.path(res_dir, "forward_selection.csv"), row.names = FALSE)
min_df <- do.call(rbind, min_rows)
write.csv(min_df, file.path(res_dir, "minimum_panels.csv"), row.names = FALSE)

# Export minimal panel derived from DONOR cohort for progression scoring in Step 08
donor_runs <- min_df[min_df$set == "DONOR", ]
if (nrow(donor_runs)) {
  writeLines(strsplit(donor_runs$genes[which.max(donor_runs$median_auc)], ";")[[1]],
             file.path(res_dir, "panel_for_ladder.txt"))
}