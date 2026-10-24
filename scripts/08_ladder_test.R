# ==============================================================================
# Script: 08_ladder_test.R
# Purpose: Score the oral carcinogenesis progression ladder (donor < contralateral < 
#          adjacent < matched < dysplasia < cancer) with a locked DONOR-derived panel.
#          Evaluates progression via permutation-based Jonckheere-Terpstra tests and
#          performs within-series paired contrasts to control for study batch confounds.
# Reference: Section III-A, Section IV-F (Ordered Scores and Within-Deposit Contrasts)
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

library(xgboost); library(pROC); library(ggplot2)

TRAIN_SET  <- "DONOR"   # Held-out testing requires DONOR cohort training
N_PERM     <- 10000     # Permutations for Jonckheere-Terpstra p-value
N_MODELS   <- 50        # Number of bootstrap models for score ensembling
ALLOW_FULL <- FALSE
MIN_GROUP  <- 3         # Minimum sample size threshold per tier

if (TRAIN_SET != "DONOR" && !ALLOW_FULL)
  stop("TRAIN_SET must be DONOR to ensure patient-derived normal tissues remain held out.")

manifest <- load_manifest()
panel_file <- file.path(res_dir, "panel_for_ladder.txt")
if (!file.exists(panel_file)) stop("Missing panel file: ", panel_file)
panel <- readLines(panel_file)
panel <- panel[nzchar(panel)]

# ------------------------------------------------------------------------------
# 1. Non-parametric Jonckheere-Terpstra test
# ------------------------------------------------------------------------------
jt_statistic <- function(x, g) {
  lv <- levels(g); total <- 0
  for (i in seq_along(lv)) for (j in seq_along(lv)) {
    if (j <= i) next
    a <- x[g == lv[i]]; b <- x[g == lv[j]]
    na <- length(a); nb <- length(b)
    if (!na || !nb) next
    r <- rank(c(a, b))
    total <- total + sum(r[(na + 1):(na + nb)]) - nb * (nb + 1) / 2
  }
  total
}

jt_test <- function(x, g, n_perm = N_PERM) {
  g <- droplevels(g)
  obs <- jt_statistic(x, g)
  perm <- vapply(seq_len(n_perm), function(i) jt_statistic(x, factor(sample(as.character(g)), levels = levels(g))),
                 numeric(1))
  list(statistic = obs, p_value = (1 + sum(perm >= obs)) / (n_perm + 1),
       z = (obs - mean(perm)) / sd(perm), n_perm = n_perm)
}

# ------------------------------------------------------------------------------
# 2. Progression scoring across preprocessing arms
# ------------------------------------------------------------------------------
score_rows <- list(); test_rows <- list(); ws_rows <- list(); ts_rows <- list()

for (arm in names(ARMS)) {
  expr <- load_arm(arm)
  genes <- intersect(panel, rownames(expr))
  if (length(genes) < 2) next

  train_samp <- intersect(set_samples(manifest, TRAIN_SET), colnames(expr))
  y_train <- set_labels(manifest, train_samp)
  X_train <- t(expr[genes, train_samp, drop = FALSE])
  X_all   <- t(expr[genes, manifest$sample, drop = FALSE])

  message("\nScoring progression across Arm ", arm, " (", ARM_LABEL[[arm]], ")")

  # Ensemble bagging: average predicted probabilities over N_MODELS bootstrap iterations
  yb <- as.integer(y_train == "cancer")
  preds <- matrix(NA_real_, nrow = nrow(X_all), ncol = N_MODELS,
                  dimnames = list(rownames(X_all), NULL))
  for (m in seq_len(N_MODELS)) {
    sp <- stratified_split(y_train)
    dtr <- xgb.DMatrix(X_train[sp$train, , drop = FALSE], label = yb[sp$train])
    dvl <- xgb.DMatrix(X_train[sp$val,   , drop = FALSE], label = yb[sp$val])
    fit <- tryCatch(xgb.train(list(objective = "binary:logistic", eval_metric = "logloss",
                                   max_depth = 4, eta = 0.1, subsample = 0.8,
                                   colsample_bytree = 0.8, nthread = 1),
                              dtr, nrounds = 400, evals = list(val = dvl),
                              early_stopping_rounds = 20, verbose = 0),
                    error = function(e) NULL)
    if (is.null(fit)) next
    preds[, m] <- predict(fit, xgb.DMatrix(X_all))
  }
  score <- rowMeans(preds, na.rm = TRUE)

  df <- data.frame(sample = manifest$sample, study = manifest$study,
                   class = manifest$class, tier = manifest$tier,
                   in_training = manifest$sample %in% train_samp,
                   arm = arm, score = score[manifest$sample])
  score_rows[[arm]] <- df

  # Ordered testing on held-out tissues and comparison with full compendium
  run_jt <- function(sub, label) {
    sub$tier <- droplevels(factor(sub$tier, levels = LADDER_ORDER))
    keep <- names(which(table(sub$tier) >= MIN_GROUP))
    sub <- sub[sub$tier %in% keep, ]; sub$tier <- droplevels(sub$tier)
    if (nlevels(sub$tier) < 3) return(NULL)
    tt <- jt_test(sub$score, sub$tier)
    data.frame(arm = arm, scope = label, n_tiers = nlevels(sub$tier),
               tiers = paste(levels(sub$tier), collapse = "<"), n = nrow(sub),
               jt_statistic = tt$statistic, z = round(tt$z, 3),
               p_value = signif(tt$p_value, 4), n_perm = tt$n_perm)
  }

  r_held <- run_jt(df[!df$in_training, ], "held-out only")
  r_plus <- run_jt(df[!df$in_training | df$tier == "cancer", ], "plus cancer")
  test_rows[[paste0(arm, "_held")]] <- r_held
  test_rows[[paste0(arm, "_plus")]] <- r_plus

  # Within-series unconfounded pairwise contrasts
  for (st in unique(df$study)) {
    g <- df[df$study == st, ]
    tiers <- names(which(table(droplevels(factor(g$tier, levels = LADDER_ORDER))) >= MIN_GROUP))
    tiers <- LADDER_ORDER[LADDER_ORDER %in% tiers]
    if (length(tiers) < 2) next
    for (i in seq_len(length(tiers) - 1L)) for (j in seq(i + 1L, length(tiers))) {
      lo <- g$score[g$tier == tiers[i]]; hi <- g$score[g$tier == tiers[j]]
      held_lo <- !all(g$in_training[g$tier == tiers[i]])
      held_hi <- !all(g$in_training[g$tier == tiers[j]])
      if (!held_lo && !held_hi) next
      w <- suppressWarnings(wilcox.test(hi, lo, alternative = "greater"))
      ws_rows[[length(ws_rows) + 1]] <- data.frame(
        arm = arm, study = st, lower_tier = tiers[i], upper_tier = tiers[j],
        n_lower = length(lo), n_upper = length(hi),
        median_lower = round(median(lo), 4), median_upper = round(median(hi), 4),
        lower_held_out = held_lo, upper_held_out = held_hi,
        p_value = signif(w$p.value, 4))
    }
  }

  for (t in levels(df$tier)) {
    g <- df[df$tier == t, ]
    if (!nrow(g)) next
    for (st in unique(g$study))
      ts_rows[[length(ts_rows) + 1]] <- data.frame(
        arm = arm, tier = t, study = st, n = sum(g$study == st),
        median_score = round(median(g$score[g$study == st]), 4))
  }

  p <- ggplot(df, aes(tier, score, fill = tier)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
    labs(title = paste0("Tissue Progression Ladder - Arm ", arm),
         subtitle = paste0(ARM_LABEL[[arm]], " (trained on cancer vs. donor normal)"),
         x = NULL, y = "Predicted Carcinoma Probability Score") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "plain"),
          axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(fig_dir, sprintf("ladder_%s.pdf", arm)), p, width = 7.5, height = 5)
  ggsave(file.path(fig_dir, sprintf("ladder_%s.png", arm)), p, width = 7.5, height = 5, dpi = 300)
}

write.csv(do.call(rbind, score_rows), file.path(res_dir, "ladder_scores.csv"), row.names = FALSE)
if (length(test_rows)) write.csv(do.call(rbind, test_rows), file.path(res_dir, "ladder_test.csv"), row.names = FALSE)
if (length(ws_rows)) write.csv(do.call(rbind, ws_rows), file.path(res_dir, "ladder_within_study.csv"), row.names = FALSE)
if (length(ts_rows)) {
  ts_df <- do.call(rbind, ts_rows)
  ts_df <- ts_df[order(ts_df$arm, match(ts_df$tier, LADDER_ORDER), ts_df$study), ]
  write.csv(ts_df, file.path(res_dir, "ladder_tier_by_study.csv"), row.names = FALSE)
}