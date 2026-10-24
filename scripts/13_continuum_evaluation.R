# ==============================================================================
# Script: 13_continuum_evaluation.R
# Purpose: Multiclass evaluation of the three-state continuum panel across 200
#          stratified splits, Shapley-guided forward selection of the minimal panel,
#          and external projection of held-out dysplasia and field mucosa.
# Reference: Paper Section III-E, Section IV-D, Section IV-E, Figure 2, Table III
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

suppressPackageStartupMessages({
  library(xgboost); library(pROC); library(ggplot2)
})

N_ITER        <- 200
AUC_TOLERANCE <- 0.01
CONTINUUM     <- c("donor", "field", "cancer")

manifest <- load_manifest()
manifest$continuum <- local({
  lab <- rep(NA_character_, nrow(manifest))
  lab[manifest$class == "cancer"] <- "cancer"
  lab[manifest$class == "normal" & manifest$normal_provenance == "donor"] <- "donor"
  lab[manifest$class == "normal" & manifest$normal_provenance != "donor"] <- "field"
  factor(lab, levels = CONTINUUM)
})

dysplasia_samples <- manifest$sample[manifest$class == "dysplasia"]
panels <- readRDS(file.path(res_dir, "continuum_panels.rds"))

# ------------------------------------------------------------------------------
# 1. Multiclass performance metrics
# ------------------------------------------------------------------------------
macro_auc <- function(prob, truth) {
  aucs <- vapply(seq_along(CONTINUUM), function(k) {
    yk <- as.integer(truth == CONTINUUM[k])
    if (length(unique(yk)) < 2) return(NA_real_)
    as.numeric(pROC::roc(yk, prob[, k], quiet = TRUE)$auc)
  }, numeric(1))
  mean(aucs, na.rm = TRUE)
}

bal_acc <- function(pred, truth) {
  mean(vapply(CONTINUUM, function(k) {
    i <- truth == k
    if (!any(i)) return(NA_real_)
    mean(pred[i] == k)
  }, numeric(1)), na.rm = TRUE)
}

softprob_matrix <- function(raw, n, K) {
  p <- if (is.matrix(raw)) raw else matrix(raw, nrow = n, ncol = K, byrow = TRUE)
  stopifnot(nrow(p) == n, ncol(p) == K)
  if (max(abs(rowSums(p) - 1)) > 1e-4) {
    stop("multi:softprob outputs do not sum to 1.")
  }
  p
}

fit_eval <- function(X, y, genes, n_iter = N_ITER, collect_cm = FALSE) {
  yi <- as.integer(factor(y, levels = CONTINUUM)) - 1L
  aucs <- numeric(0); baccs <- numeric(0)
  K <- length(CONTINUUM)
  cm <- matrix(0, K, K, dimnames = list(CONTINUUM, CONTINUUM))
  for (it in seq_len(n_iter)) {
    sp <- stratified_split(factor(y, levels = CONTINUUM))
    dtr <- xgb.DMatrix(X[sp$train, genes, drop = FALSE], label = yi[sp$train])
    dvl <- xgb.DMatrix(X[sp$val,   genes, drop = FALSE], label = yi[sp$val])
    fit <- tryCatch(xgb.train(list(objective = "multi:softprob", num_class = K,
                                   eval_metric = "mlogloss", max_depth = 4,
                                   eta = 0.1, subsample = 0.8,
                                   colsample_bytree = 0.8, nthread = 1),
                              dtr, nrounds = 400, evals = list(val = dvl),
                              early_stopping_rounds = 20, verbose = 0),
                    error = function(e) NULL)
    if (is.null(fit)) next
    p <- softprob_matrix(predict(fit, xgb.DMatrix(X[sp$test, genes, drop = FALSE])),
                         length(sp$test), K)
    truth <- as.character(y[sp$test])
    pred <- CONTINUUM[max.col(p)]
    aucs <- c(aucs, macro_auc(p, truth)); baccs <- c(baccs, bal_acc(pred, truth))
    if (collect_cm)
      for (i in seq_along(truth)) cm[truth[i], pred[i]] <- cm[truth[i], pred[i]] + 1
  }
  list(auc = aucs, bacc = baccs, cm = cm)
}

# ------------------------------------------------------------------------------
# 2. Evaluation and forward selection per preprocessing arm
# ------------------------------------------------------------------------------
eval_rows <- list(); cm_rows <- list(); min_rows <- list()
field_rows <- list(); dys_rows <- list()

for (arm in names(ARMS)) {
  run <- paste("multinomial", arm, sep = "_")
  if (!run %in% names(panels)) next
  expr <- load_arm(arm)
  samp <- intersect(manifest$sample[!is.na(manifest$continuum)], colnames(expr))
  y <- manifest$continuum[match(samp, manifest$sample)]
  X <- t(expr[, samp, drop = FALSE])
  X_dys <- t(expr[, intersect(dysplasia_samples, colnames(expr)), drop = FALSE])
  
  genes <- intersect(panels[[run]], colnames(X))
  if (length(genes) < 2) next

  r <- fit_eval(X, y, genes, collect_cm = TRUE)
  q <- function(v) quantile(v, c(0.5, 0.025, 0.975), na.rm = TRUE)
  qa <- q(r$auc); qb <- q(r$bacc)
  
  eval_rows[[arm]] <- data.frame(arm = arm, arm_label = ARM_LABEL[[arm]],
                                 n_genes = length(genes), n_iter = length(r$auc),
                                 auc = round(qa[1], 4), auc_lo = round(qa[2], 4),
                                 auc_hi = round(qa[3], 4), bacc = round(qb[1], 4),
                                 bacc_lo = round(qb[2], 4), bacc_hi = round(qb[3], 4))
  cm_rows[[arm]] <- data.frame(arm = arm, as.data.frame(as.table(r$cm)))
  
  # Forward selection in Shapley importance order
  sp <- stratified_split(y)
  yi <- as.integer(y) - 1L
  dtr <- xgb.DMatrix(X[sp$train, genes, drop = FALSE], label = yi[sp$train])
  dvl <- xgb.DMatrix(X[sp$val, genes, drop = FALSE], label = yi[sp$val])
  fit <- xgb.train(list(objective = "multi:softprob", num_class = length(CONTINUUM),
                        eval_metric = "mlogloss", max_depth = 4, eta = 0.1,
                        subsample = 0.8, colsample_bytree = 0.8, nthread = 1),
                   dtr, nrounds = 400, evals = list(val = dvl),
                   early_stopping_rounds = 20, verbose = 0)
  contrib <- predict(fit, xgb.DMatrix(X[sp$test, genes, drop = FALSE]), predcontrib = TRUE)
  imp <- if (is.list(contrib)) {
    Reduce(`+`, lapply(contrib, function(m) colMeans(abs(m[, seq_along(genes), drop = FALSE]))))
  } else if (length(dim(contrib)) == 3L) {
    if (dim(contrib)[2] == length(genes) + 1L) {
      apply(abs(contrib[, seq_along(genes), , drop = FALSE]), 2, mean) * dim(contrib)[3]
    } else {
      apply(abs(contrib[, , seq_along(genes), drop = FALSE]), 3, mean) * dim(contrib)[2]
    }
  } else {
    colMeans(abs(contrib[, seq_along(genes), drop = FALSE]))
  }
  
  gene_order <- genes[order(imp, decreasing = TRUE)]
  full_auc <- median(r$auc, na.rm = TRUE)
  chosen <- NULL
  for (k in seq_along(gene_order)) {
    rk <- fit_eval(X, y, gene_order[seq_len(k)], n_iter = 60)
    ak <- median(rk$auc, na.rm = TRUE); bk <- median(rk$bacc, na.rm = TRUE)
    min_rows[[paste(arm, k)]] <- data.frame(arm = arm, size = k,
                                            gene_added = gene_order[k],
                                            auc = round(ak, 4), bacc = round(bk, 4))
    if (is.null(chosen) && ak >= full_auc - AUC_TOLERANCE) chosen <- k
  }
  if (!is.null(chosen)) {
    writeLines(gene_order[seq_len(chosen)],
               file.path(res_dir, sprintf("continuum_minimum_%s.txt", arm)))
  }
  
  # External projection of held-out dysplasia cohort
  if (nrow(X_dys)) {
    dp <- matrix(0, nrow(X_dys), length(CONTINUUM))
    n_ok <- 0
    for (m in seq_len(50)) {
      spd <- stratified_split(y)
      dtr2 <- xgb.DMatrix(X[spd$train, genes, drop = FALSE], label = yi[spd$train])
      dvl2 <- xgb.DMatrix(X[spd$val, genes, drop = FALSE], label = yi[spd$val])
      f3 <- tryCatch(xgb.train(list(objective = "multi:softprob", num_class = length(CONTINUUM),
                                    eval_metric = "mlogloss", max_depth = 4, eta = 0.1,
                                    subsample = 0.8, colsample_bytree = 0.8, nthread = 1),
                               dtr2, nrounds = 400, evals = list(val = dvl2),
                               early_stopping_rounds = 20, verbose = 0),
                     error = function(e) NULL)
      if (is.null(f3)) next
      dp <- dp + softprob_matrix(predict(f3, xgb.DMatrix(X_dys[, genes, drop = FALSE])),
                                 nrow(X_dys), length(CONTINUUM))
      n_ok <- n_ok + 1
    }
    dp <- dp / max(n_ok, 1)
    assign_to <- CONTINUUM[max.col(dp)]
    tb <- table(factor(assign_to, levels = CONTINUUM))
    dys_rows[[arm]] <- data.frame(arm = arm, n_dysplasia = nrow(X_dys),
                                  as.list(setNames(as.integer(tb), paste0("n_", CONTINUUM))),
                                  as.list(setNames(round(colMeans(dp), 4),
                                                   paste0("mean_p_", CONTINUUM))))
  }
  
  # Binary classifier evaluation on held-out field mucosa
  bin <- samp[y %in% c("donor", "cancer")]
  yb <- as.integer(manifest$continuum[match(bin, manifest$sample)] == "cancer")
  fld <- samp[y == "field"]
  probs <- matrix(NA_real_, length(fld), 50)
  for (m in seq_len(50)) {
    spb <- stratified_split(factor(yb))
    d1 <- xgb.DMatrix(X[bin[spb$train], genes, drop = FALSE], label = yb[spb$train])
    d2 <- xgb.DMatrix(X[bin[spb$val], genes, drop = FALSE], label = yb[spb$val])
    f2 <- tryCatch(xgb.train(list(objective = "binary:logistic", eval_metric = "logloss",
                                  max_depth = 4, eta = 0.1, subsample = 0.8,
                                  colsample_bytree = 0.8, nthread = 1),
                             d1, nrounds = 400, evals = list(val = d2),
                             early_stopping_rounds = 20, verbose = 0),
                   error = function(e) NULL)
    if (is.null(f2)) next
    probs[, m] <- predict(f2, xgb.DMatrix(X[fld, genes, drop = FALSE]))
  }
  s <- rowMeans(probs, na.rm = TRUE)
  field_rows[[arm]] <- data.frame(arm = arm, n_field = length(fld),
                                  above_half = sum(s > 0.5),
                                  above_nine = sum(s > 0.9),
                                  median_score = round(median(s), 4))
  
  # Row-normalized confusion matrix figure
  cmdf <- as.data.frame(as.table(r$cm / rowSums(r$cm)))
  names(cmdf) <- c("truth", "predicted", "fraction")
  cmdf$truth <- factor(cmdf$truth, levels = rev(CONTINUUM))
  cmdf$predicted <- factor(cmdf$predicted, levels = CONTINUUM)
  p <- ggplot(cmdf, aes(predicted, truth, fill = fraction)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", fraction)), size = 3.4, colour = "#0b0b0b") +
    scale_fill_gradient(low = "white", high = "#2a78d6", limits = c(0, 1)) +
    labs(title = sprintf("Arm %s, %s", arm, ARM_LABEL[[arm]])) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
  ggsave(file.path(fig_dir, sprintf("continuum_confusion_%s.pdf", arm)), p, width = 4, height = 3.6)
  ggsave(file.path(fig_dir, sprintf("continuum_confusion_%s.png", arm)), p, width = 4, height = 3.6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 3. Export validation summaries
# ------------------------------------------------------------------------------
for (nm in c("continuum_eval", "continuum_confusion", "continuum_minimum",
             "field_assignment", "dysplasia_assignment")) {
  rows <- switch(nm, continuum_eval = eval_rows, continuum_confusion = cm_rows,
                 continuum_minimum = min_rows, field_assignment = field_rows,
                 dysplasia_assignment = dys_rows)
  if (length(rows)) {
    write.csv(do.call(rbind, rows), file.path(res_dir, paste0(nm, ".csv")), row.names = FALSE)
  }
}