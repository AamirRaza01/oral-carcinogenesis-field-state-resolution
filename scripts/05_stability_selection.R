# ==============================================================================
# Script: 05_stability_selection.R
# Purpose: Bootstrap stability selection (500 iterations) via elastic net 
#          (alpha = 0.5) across 6 experimental configurations (3 preprocessing 
#          arms x 2 normal tissue definitions).
# Reference: Paper Sections III-C, III-D, and Table II
# ==============================================================================

source("00_oral_common.R")
set.seed(123)

library(glmnet)
library(parallel)

N_BOOT    <- 500
ALPHA     <- 0.5     # Elastic net mixing parameter
PI_THRESH <- 0.90    # Selection frequency threshold for stability
N_FOLDS   <- 5
N_WORKERS <- min(4L, max(1L, parallel::detectCores() - 1L))

manifest <- load_manifest()
grid <- run_grid()

# ------------------------------------------------------------------------------
# 1. Stability selection engine for an expression matrix and outcome
# ------------------------------------------------------------------------------
run_stability <- function(X, y, n_boot = N_BOOT) {
  genes <- colnames(X)

  boot_one <- function(b) {
    idx <- unlist(lapply(levels(y), function(lev) {
      i <- which(y == lev); sample(i, length(i), replace = TRUE)
    }))
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
    if (nlevels(droplevels(yb)) < 2) return(NULL)

    fold <- rep(NA_integer_, length(yb))
    for (lev in levels(yb)) {
      li <- which(yb == lev)
      fold[li] <- sample(rep(seq_len(N_FOLDS), length.out = length(li)))
    }
    fit <- tryCatch(cv.glmnet(Xb, yb, family = "binomial", alpha = ALPHA,
                              type.measure = "class", foldid = fold),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    co <- coef(fit, s = "lambda.1se")
    setdiff(rownames(co)[which(co[, 1] != 0)], "(Intercept)")
  }

  # Parallel execution: PSOCK cluster on Windows, forked processes on Unix
  res <- if (N_WORKERS <= 1L) {
    lapply(seq_len(n_boot), boot_one)
  } else if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(N_WORKERS)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
    parallel::clusterExport(cl, c("X", "y"), envir = environment())
    parallel::clusterExport(cl, c("ALPHA", "N_FOLDS"), envir = globalenv())
    parallel::clusterSetRNGStream(cl, 123)
    parallel::parLapply(cl, seq_len(n_boot), boot_one)
  } else {
    mclapply(seq_len(n_boot), boot_one, mc.cores = N_WORKERS, mc.set.seed = TRUE)
  }

  counts <- setNames(rep(0, length(genes)), genes)
  eff <- sum(!vapply(res, is.null, logical(1)))
  for (r in res) if (!is.null(r)) counts[r] <- counts[r] + 1
  list(freq = counts / eff, counts = counts, effective = eff)
}

# ------------------------------------------------------------------------------
# 2. Execute stability selection across configuration grid
# ------------------------------------------------------------------------------
all_freq <- list(); panels <- list(); summary_rows <- list()

for (i in seq_len(nrow(grid))) {
  arm <- grid$arm[i]; set <- grid$set[i]; run <- grid$run[i]
  message("\nProcessing configuration: ", run, " (", ARM_LABEL[[arm]], ", set: ", set, ")")

  expr <- load_arm(arm)
  samp <- intersect(set_samples(manifest, set), colnames(expr))
  X <- t(expr[, samp, drop = FALSE])
  y <- set_labels(manifest, samp)

  st <- run_stability(X, y)
  stable <- names(st$freq)[st$freq >= PI_THRESH]
  message(sprintf("  Stable genes (pi >= %.2f): %d", PI_THRESH, length(stable)))

  all_freq[[run]] <- st$freq
  panels[[run]]   <- stable
  summary_rows[[run]] <- data.frame(run = run, arm = arm, set = set,
                                    n_arrays = nrow(X), n_cancer = sum(y == "cancer"),
                                    n_normal = sum(y == "normal"),
                                    n_stable = length(stable), effective_boots = st$effective)

  write.csv(data.frame(gene = names(st$freq), pi = as.numeric(st$freq),
                       count = as.integer(st$counts)),
            file.path(res_dir, sprintf("stability_%s.csv", run)), row.names = FALSE)
}

saveRDS(panels,   file.path(res_dir, "stability_panels.rds"))
saveRDS(all_freq, file.path(res_dir, "stability_frequencies.rds"))
summ <- do.call(rbind, summary_rows); rownames(summ) <- NULL
write.csv(summ, file.path(res_dir, "stability_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Inter-run concordance analysis
# ------------------------------------------------------------------------------
# Evaluate gene ranking agreement via Spearman correlation of selection frequencies
# and discrete panel overlap via Jaccard similarity.
runs <- names(all_freq)
genes_common <- Reduce(intersect, lapply(all_freq, names))
agree <- list()

for (a in seq_along(runs)) for (b in seq_along(runs)) {
  if (b <= a) next
  ra <- runs[a]; rb <- runs[b]
  rho <- cor(all_freq[[ra]][genes_common], all_freq[[rb]][genes_common], method = "spearman")
  ia <- panels[[ra]]; ib <- panels[[rb]]
  jac <- if (length(union(ia, ib))) length(intersect(ia, ib)) / length(union(ia, ib)) else NA_real_
  agree[[length(agree) + 1]] <- data.frame(
    run_1 = ra, run_2 = rb, spearman_pi = round(rho, 4),
    jaccard_panel = round(jac, 4), n_shared = length(intersect(ia, ib)),
    n_only_1 = length(setdiff(ia, ib)), n_only_2 = length(setdiff(ib, ia)))
}

agree_df <- do.call(rbind, agree)
agree_df <- agree_df[order(-agree_df$spearman_pi), ]
write.csv(agree_df, file.path(res_dir, "stability_agreement.csv"), row.names = FALSE)

tally <- table(unlist(panels))
write.csv(data.frame(gene = names(tally), n_runs_stable = as.integer(tally)),
          file.path(res_dir, "stability_gene_tally.csv"), row.names = FALSE)