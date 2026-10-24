# ==============================================================================
# Script: 12_continuum_panel.R
# Purpose: Identify biomarker panels resolving the three-state oral continuum
#          (donor mucosa, field mucosa, carcinoma) via multinomial, one-vs-rest,
#          and pairwise stability selection across preprocessing arms.
# Reference: Paper Section III-C, Section IV-B, Table II
# ==============================================================================

source("00_oral_common.R")
set.seed(123)

suppressPackageStartupMessages({ library(glmnet); library(parallel) })

N_BOOT    <- 500
ALPHA     <- 0.5
PI_THRESH <- 0.90
N_FOLDS   <- 5
N_WORKERS <- min(4L, max(1L, parallel::detectCores() - 1L))

CONTINUUM <- c("donor", "field", "cancer")

CONFIGS <- c("multinomial",
             "ovr_field",
             "pair_donor_field", "pair_field_cancer")

manifest <- load_manifest()

# ------------------------------------------------------------------------------
# 1. Three-state continuum labeling (dysplasia held out)
# ------------------------------------------------------------------------------
continuum_label <- function(m) {
  lab <- rep(NA_character_, nrow(m))
  lab[m$class == "cancer"] <- "cancer"
  lab[m$class == "normal" & m$normal_provenance == "donor"] <- "donor"
  lab[m$class == "normal" & m$normal_provenance != "donor"] <- "field"
  factor(lab, levels = CONTINUUM)
}

manifest$continuum <- continuum_label(manifest)

# ------------------------------------------------------------------------------
# 2. Resampled stability selection engine
# ------------------------------------------------------------------------------
run_stability <- function(X, y, family, n_boot = N_BOOT) {
  genes <- colnames(X)

  boot_one <- function(b) {
    idx <- unlist(lapply(levels(y), function(lev) {
      i <- which(y == lev); if (!length(i)) return(integer(0))
      sample(i, length(i), replace = TRUE)
    }))
    Xb <- X[idx, , drop = FALSE]; yb <- droplevels(y[idx])
    if (nlevels(yb) < 2) return(NULL)

    fold <- rep(NA_integer_, length(yb))
    for (lev in levels(yb)) {
      li <- which(yb == lev)
      fold[li] <- sample(rep(seq_len(N_FOLDS), length.out = length(li)))
    }
    fit <- tryCatch(cv.glmnet(Xb, yb, family = family, alpha = ALPHA,
                              type.measure = "class", foldid = fold),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    co <- coef(fit, s = "lambda.1se")
    sel <- if (is.list(co))
      unique(unlist(lapply(co, function(m) rownames(m)[which(as.numeric(m) != 0)])))
    else rownames(co)[which(co[, 1] != 0)]
    setdiff(sel, "(Intercept)")
  }

  res <- if (N_WORKERS <= 1L) {
    lapply(seq_len(n_boot), boot_one)
  } else if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(N_WORKERS)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
    parallel::clusterExport(cl, c("X", "y", "family"), envir = environment())
    parallel::clusterExport(cl, c("ALPHA", "N_FOLDS"), envir = globalenv())
    parallel::clusterSetRNGStream(cl, 123)
    parallel::parLapply(cl, seq_len(n_boot), boot_one)
  } else {
    mclapply(seq_len(n_boot), boot_one, mc.cores = N_WORKERS, mc.set.seed = TRUE)
  }

  counts <- setNames(rep(0, length(genes)), genes)
  eff <- sum(!vapply(res, is.null, logical(1)))
  for (r in res) if (!is.null(r)) counts[r] <- counts[r] + 1
  list(freq = counts / eff, effective = eff)
}

config_data <- function(cfg, samples, lab) {
  if (cfg == "multinomial")
    return(list(keep = rep(TRUE, length(samples)), y = lab, family = "multinomial"))
  if (grepl("^ovr_", cfg)) {
    target <- sub("^ovr_", "", cfg)
    y <- factor(ifelse(lab == target, target, "rest"), levels = c("rest", target))
    return(list(keep = rep(TRUE, length(samples)), y = y, family = "binomial"))
  }
  parts <- strsplit(sub("^pair_", "", cfg), "_")[[1]]
  keep <- lab %in% parts
  list(keep = keep, y = factor(as.character(lab[keep]), levels = parts), family = "binomial")
}

# ------------------------------------------------------------------------------
# 3. Model fitting across configurations and arms
# ------------------------------------------------------------------------------
all_freq <- list(); panels <- list(); summary_rows <- list()

for (arm in names(ARMS)) {
  expr <- load_arm(arm)
  samp <- intersect(manifest$sample[!is.na(manifest$continuum)], colnames(expr))
  lab_all <- manifest$continuum[match(samp, manifest$sample)]

  for (cfg in CONFIGS) {
    run <- paste(cfg, arm, sep = "_")
    cd <- config_data(cfg, samp, lab_all)
    s <- samp[cd$keep]
    y <- droplevels(cd$y)
    if (nlevels(y) < 2) next

    X <- t(expr[, s, drop = FALSE])
    st <- run_stability(X, y, cd$family)
    stable <- names(st$freq)[st$freq >= PI_THRESH]

    all_freq[[run]] <- st$freq
    panels[[run]] <- stable
    summary_rows[[run]] <- data.frame(
      run = run, config = cfg, arm = arm, n_arrays = nrow(X),
      n_classes = nlevels(y), smallest_class = min(as.integer(table(y))),
      n_stable = length(stable), effective_boots = st$effective)

    write.csv(data.frame(gene = names(st$freq), pi = as.numeric(st$freq)),
              file.path(res_dir, sprintf("continuum_%s.csv", run)), row.names = FALSE)
  }
}

saveRDS(panels,   file.path(res_dir, "continuum_panels.rds"))
saveRDS(all_freq, file.path(res_dir, "continuum_frequencies.rds"))
write.csv(do.call(rbind, summary_rows), file.path(res_dir, "continuum_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Export configuration union and multinomial subsets
# ------------------------------------------------------------------------------
for (arm in names(ARMS)) {
  runs <- names(panels)[grepl(paste0("_", arm, "$"), names(panels))]
  if (!length(runs)) next
  uni <- sort(unique(unlist(panels[runs])))
  mn <- panels[[paste("multinomial", arm, sep = "_")]]
  writeLines(uni, file.path(res_dir, sprintf("continuum_union_%s.txt", arm)))
  if (length(mn)) writeLines(mn, file.path(res_dir, sprintf("continuum_multinomial_%s.txt", arm)))
}