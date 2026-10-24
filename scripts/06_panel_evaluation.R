# ==============================================================================
# Script: 06_panel_evaluation.R
# Purpose: Multi-model validation (Elastic Net, XGBoost, Random Forest, Logistic
#          Regression) over 1,000 stratified resampled splits. Optimal decision
#          thresholds are determined on validation splits via Youden's J statistic.
# Reference: Section III-E (Evaluation) and Table III
# ==============================================================================

source("00_oral_common.R")
set.seed(123)

library(glmnet); library(xgboost); library(randomForest); library(pROC)

N_ITER <- 1000

manifest <- load_manifest()
panels <- readRDS(file.path(res_dir, "stability_panels.rds"))
grid <- run_grid()

# ------------------------------------------------------------------------------
# 1. Performance metrics computation
# ------------------------------------------------------------------------------
met_bin <- function(y_true, y_pred_class, y_prob, pos = "cancer") {
  tp <- sum(y_pred_class == pos & y_true == pos); fp <- sum(y_pred_class == pos & y_true != pos)
  tn <- sum(y_pred_class != pos & y_true != pos); fn <- sum(y_pred_class != pos & y_true == pos)
  se <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  sp <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  pc <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1 <- if (!is.na(pc) && !is.na(se) && (pc + se) > 0) 2 * pc * se / (pc + se) else NA_real_
  auc <- tryCatch(as.numeric(roc(as.numeric(y_true == pos), y_prob, quiet = TRUE)$auc),
                  error = function(e) NA_real_)
  c(acc = mean(y_pred_class == y_true), bacc = mean(c(se, sp), na.rm = TRUE),
    auc = auc, sens = se, spec = sp, prec = pc, f1 = f1)
}

youden <- function(y_val, p_val, pos = "cancer") {
  r <- tryCatch(roc(as.numeric(y_val == pos), p_val, quiet = TRUE), error = function(e) NULL)
  if (is.null(r)) return(0.5)
  th <- r$thresholds[which.max(r$sensitivities + r$specificities - 1)]
  if (!is.finite(th)) 0.5 else th
}

# ------------------------------------------------------------------------------
# 2. Resampled iteration workflow across model families
# ------------------------------------------------------------------------------
run_one_iter <- function(X, y) {
  sp <- stratified_split(y)
  Xtr <- X[sp$train, , drop = FALSE]; ytr <- y[sp$train]
  Xvl <- X[sp$val,   , drop = FALSE]; yvl <- y[sp$val]
  Xte <- X[sp$test,  , drop = FALSE]; yte <- y[sp$test]
  pos <- "cancer"; res <- list()

  fold <- rep(NA_integer_, length(ytr))
  for (lev in levels(ytr)) {
    li <- which(ytr == lev); fold[li] <- sample(rep(1:5, length.out = length(li)))
  }

  tryCatch({
    fit <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0.5,
                     type.measure = "class", foldid = fold)
    pb <- predict(fit, Xte, s = "lambda.1se", type = "response")[, 1]
    pv <- predict(fit, Xvl, s = "lambda.1se", type = "response")[, 1]
    th <- youden(yvl, pv)
    res$EN <- met_bin(yte, ifelse(pb >= th, pos, "normal"), pb)
  }, error = function(e) {})

  tryCatch({
    dtr <- xgb.DMatrix(Xtr, label = as.integer(ytr == pos))
    dvl <- xgb.DMatrix(Xvl, label = as.integer(yvl == pos))
    par <- list(objective = "binary:logistic", eval_metric = "logloss",
                max_depth = 6, eta = 0.1, subsample = 0.8,
                colsample_bytree = 0.8, nthread = 1)
    fit <- xgb.train(par, dtr, nrounds = 500, evals = list(val = dvl),
                     early_stopping_rounds = 20, verbose = 0)
    pb <- predict(fit, xgb.DMatrix(Xte)); pv <- predict(fit, dvl)
    th <- youden(yvl, pv)
    res$XGB <- met_bin(yte, ifelse(pb >= th, pos, "normal"), pb)
  }, error = function(e) {})

  tryCatch({
    fit <- randomForest(Xtr, ytr, ntree = 500)
    pb <- predict(fit, Xte, type = "prob")[, pos]
    pv <- predict(fit, Xvl, type = "prob")[, pos]
    th <- youden(yvl, pv)
    res$RF <- met_bin(yte, ifelse(pb >= th, pos, "normal"), pb)
  }, error = function(e) {})

  tryCatch({
    df_tr <- data.frame(Xtr, .y = ytr, check.names = TRUE)
    fit <- suppressWarnings(glm(.y ~ ., data = df_tr, family = binomial))
    pb <- predict(fit, data.frame(Xte, check.names = TRUE), type = "response")
    pv <- predict(fit, data.frame(Xvl, check.names = TRUE), type = "response")
    th <- youden(yvl, pv)
    res$LR <- met_bin(yte, ifelse(pb >= th, pos, "normal"), pb)
  }, error = function(e) {})

  res
}

# ------------------------------------------------------------------------------
# 3. Model evaluation across all configurations
# ------------------------------------------------------------------------------
per_iter <- list(); summ_rows <- list()
metric_names <- c("acc", "bacc", "auc", "sens", "spec", "prec", "f1")

for (i in seq_len(nrow(grid))) {
  arm <- grid$arm[i]; set <- grid$set[i]; run <- grid$run[i]
  panel <- panels[[run]]
  if (is.null(panel) || length(panel) < 2) next

  expr <- load_arm(arm)
  samp <- intersect(set_samples(manifest, set), colnames(expr))
  genes <- intersect(panel, rownames(expr))
  X <- t(expr[genes, samp, drop = FALSE]); y <- set_labels(manifest, samp)

  message("\nEvaluating run: ", run, " (", length(genes), " genes)")
  all_res <- lapply(seq_len(N_ITER), function(k) run_one_iter(X, y))

  for (k in seq_along(all_res)) for (m in names(all_res[[k]])) {
    v <- all_res[[k]][[m]]
    per_iter[[length(per_iter) + 1]] <- data.frame(
      run = run, arm = arm, set = set, iter = k, model = m,
      as.list(v), check.names = FALSE)
  }

  for (m in c("EN", "XGB", "RF", "LR")) {
    vals <- do.call(rbind, lapply(all_res, function(r) if (!is.null(r[[m]])) r[[m]] else NULL))
    if (is.null(vals)) next
    for (metric in metric_names) {
      v <- vals[, metric]; v <- v[!is.na(v)]
      if (!length(v)) next
      summ_rows[[length(summ_rows) + 1]] <- data.frame(
        run = run, arm = arm, set = set, model = m, metric = metric, n = length(v),
        median = round(median(v), 4),
        ci_lo = round(unname(quantile(v, 0.025)), 4),
        ci_hi = round(unname(quantile(v, 0.975)), 4))
    }
  }
}

per_iter_df <- do.call(rbind, per_iter)
summ_df <- do.call(rbind, summ_rows)
write.csv(per_iter_df, file.path(res_dir, "panel_eval_per_iter.csv"), row.names = FALSE)
write.csv(summ_df,     file.path(res_dir, "panel_eval_summary.csv"), row.names = FALSE)