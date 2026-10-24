# ==============================================================================
# Script: 09_external_validation.R
# Purpose: Independent cross-platform cohort validation on GSE25099 (Affymetrix
#          Human Exon 1.0 ST Array; GPL5175) using within-sample gene rank
#          transformations and Entrez Gene ID harmonization.
# Reference: Paper Section IV-I (Cross-Platform Transfer)
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

options(timeout = 3600)

suppressPackageStartupMessages({
  library(GEOquery); library(xgboost); library(pROC); library(ggplot2)
  library(AnnotationDbi); library(hgu133plus2.db)
})

EXT_GSE    <- "GSE25099"
EXT_GPL    <- "GPL5175"
TRAIN_SET  <- "DONOR"
BEST_RUN   <- "C_DONOR"
N_MODELS   <- 50

EXPECT_NORMAL <- 22
EXPECT_CANCER <- 57

manifest <- load_manifest()
ext_dir  <- file.path(base_dir, "external"); dir.create(ext_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load candidate panels
# ------------------------------------------------------------------------------
panel_file <- file.path(res_dir, "panel_for_ladder.txt")
stab_file  <- file.path(res_dir, "stability_panels.rds")
if (!file.exists(panel_file)) stop("Missing panel file: ", panel_file)
if (!file.exists(stab_file))  stop("Missing stability panels: ", stab_file)

panel_min <- readLines(panel_file); panel_min <- panel_min[nzchar(panel_min)]
stab_all  <- readRDS(stab_file)
if (!BEST_RUN %in% names(stab_all)) stop("BEST_RUN not found in stability_panels.rds")
panel_stable <- stab_all[[BEST_RUN]]

PANELS <- list(stable = panel_stable, minimum = panel_min)

# ------------------------------------------------------------------------------
# 2. Retrieve and parse validation cohort
# ------------------------------------------------------------------------------
message("Downloading validation dataset: ", EXT_GSE)
gse  <- getGEO(EXT_GSE, GSEMatrix = TRUE, getGPL = TRUE, destdir = ext_dir)
eset <- if (is.list(gse)) gse[[1]] else gse
ext_expr <- Biobase::exprs(eset)
pd <- Biobase::pData(eset)

if (max(ext_expr, na.rm = TRUE) > 100) {
  ext_expr <- log2(ext_expr + 1)
}

# Standardize sample phenotypic labels
lab_cols <- grep("title|source_name|characteristics", names(pd), value = TRUE)
txt <- tolower(apply(pd[, lab_cols, drop = FALSE], 1, paste, collapse = " | "))
ext_class <- ifelse(grepl("normal|control|non-?tumor|healthy", txt), "normal",
             ifelse(grepl("cancer|tumor|carcinoma|oscc", txt), "cancer", NA_character_))

n_norm <- sum(ext_class == "normal", na.rm = TRUE)
n_canc <- sum(ext_class == "cancer", na.rm = TRUE)
if (n_norm != EXPECT_NORMAL || n_canc != EXPECT_CANCER) {
  warning(sprintf("Sample label counts (%d normal, %d cancer) differ from expected (%d, %d).",
                  n_norm, n_canc, EXPECT_NORMAL, EXPECT_CANCER), call. = FALSE)
}

keep <- !is.na(ext_class)
ext_expr <- ext_expr[, keep, drop = FALSE]; ext_class <- ext_class[keep]

# ------------------------------------------------------------------------------
# 3. Map external platform probes to Entrez Gene IDs
# ------------------------------------------------------------------------------
fd <- Biobase::fData(eset)
ga_col <- grep("gene_assignment", names(fd), value = TRUE, ignore.case = TRUE)[1]

ext_entrez <- NULL; ext_symbol <- NULL
if (!is.na(ga_col)) {
  first_rec <- trimws(sapply(strsplit(as.character(fd[[ga_col]]), "///", fixed = TRUE), `[`, 1))
  parts <- strsplit(first_rec, "//", fixed = TRUE)
  ext_entrez <- trimws(sapply(parts, function(p) if (length(p) >= 5) p[5] else NA_character_))
  ext_symbol <- trimws(sapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_))
  ext_entrez[!grepl("^[0-9]+$", ext_entrez)] <- NA
  ext_symbol[!nzchar(ext_symbol)] <- NA
  names(ext_entrez) <- names(ext_symbol) <- rownames(fd)
}

# Fallback to Bioconductor annotation package if platform table is unannotated
if (is.null(ext_entrez) || sum(!is.na(ext_entrez)) < 1000) {
  if (!requireNamespace("huex10sttranscriptcluster.db", quietly = TRUE)) {
    stop("Platform table lacks annotation. Install huex10sttranscriptcluster.db.")
  }
  db <- huex10sttranscriptcluster.db::huex10sttranscriptcluster.db
  ext_entrez <- AnnotationDbi::mapIds(db, keys = rownames(ext_expr), column = "ENTREZID",
                                      keytype = "PROBEID", multiVals = "first")
  ext_symbol <- AnnotationDbi::mapIds(db, keys = rownames(ext_expr), column = "SYMBOL",
                                      keytype = "PROBEID", multiVals = "first")
}

ext_entrez <- ext_entrez[rownames(ext_expr)]
ext_symbol <- ext_symbol[rownames(ext_expr)]
ok <- !is.na(ext_entrez)
ext_expr <- ext_expr[ok, , drop = FALSE]
ext_entrez <- ext_entrez[ok]; ext_symbol <- ext_symbol[ok]

# Collapse redundant probe sets by maximal mean intensity
pm <- rowMeans(ext_expr, na.rm = TRUE)
pick <- vapply(split(seq_along(ext_entrez), ext_entrez),
               function(ix) if (length(ix) == 1L) ix else ix[which.max(pm[ix])], integer(1))
ext_by_entrez <- ext_expr[pick, , drop = FALSE]
rownames(ext_by_entrez) <- names(pick)
ext_sym_of_entrez <- setNames(ext_symbol[pick], names(pick))
all_ext_symbols <- unique(na.omit(ext_symbol))

# ------------------------------------------------------------------------------
# 4. Harmonize discovery and validation feature spaces
# ------------------------------------------------------------------------------
map_file <- file.path(out_dir, "shared_symbol_selection_map.csv")
if (!file.exists(map_file)) stop("Missing symbol selection map: ", map_file)
sel <- read.csv(map_file, stringsAsFactors = FALSE)
sel$entrez <- AnnotationDbi::mapIds(hgu133plus2.db, keys = sel$probe_id,
                                    column = "ENTREZID", keytype = "PROBEID",
                                    multiVals = "first")
sel <- sel[!is.na(sel$entrez) & !duplicated(sel$entrez), ]
sym_of_entrez <- setNames(sel$gene_symbol, sel$entrez)
entrez_of_sym <- setNames(sel$entrez, sel$gene_symbol)

shared_entrez <- intersect(rownames(ext_by_entrez), names(sym_of_entrez))
ext_gene <- ext_by_entrez[shared_entrez, , drop = FALSE]
rownames(ext_gene) <- sym_of_entrez[shared_entrez]

# Record cross-platform gene representation
cov_rows <- list()
for (pn in names(PANELS)) {
  for (g in PANELS[[pn]]) {
    e <- entrez_of_sym[g]
    hit <- !is.na(e) && e %in% rownames(ext_by_entrez)
    cov_rows[[length(cov_rows) + 1]] <- data.frame(
      panel = pn, gene = g,
      entrez = if (is.na(e)) NA_character_ else unname(e),
      on_external_platform = hit,
      external_symbol = if (hit) unname(ext_sym_of_entrez[[unname(e)]]) else NA_character_,
      symbol_match = g %in% all_ext_symbols)
  }
}
coverage <- do.call(rbind, cov_rows)
write.csv(coverage, file.path(res_dir, "external_gene_coverage.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Within-sample rank transformation and model transfer
# ------------------------------------------------------------------------------
rank_within_sample <- function(M) apply(M, 2, function(x) rank(x, na.last = "keep") / sum(!is.na(x)))

res_rows <- list(); cm_rows <- list()

for (arm in names(ARMS)) {
  expr   <- load_arm(arm)
  common <- intersect(rownames(expr), rownames(ext_gene))

  disc_rank <- rank_within_sample(expr[common, , drop = FALSE])
  ext_rank  <- rank_within_sample(ext_gene[common, , drop = FALSE])

  train_samp <- intersect(set_samples(manifest, TRAIN_SET), colnames(expr))
  y_train <- set_labels(manifest, train_samp)
  yb <- as.integer(y_train == "cancer")
  y_ext <- as.integer(ext_class == "cancer")

  for (pn in names(PANELS)) {
    genes <- intersect(PANELS[[pn]], common)
    if (!length(genes)) next

    X_train <- t(disc_rank[genes, train_samp, drop = FALSE])
    X_ext   <- t(ext_rank[genes, , drop = FALSE])

    preds <- matrix(NA_real_, nrow(X_ext), N_MODELS); ths <- rep(NA_real_, N_MODELS)
    for (m in seq_len(N_MODELS)) {
      sp  <- stratified_split(y_train)
      dtr <- xgb.DMatrix(X_train[sp$train, , drop = FALSE], label = yb[sp$train])
      dvl <- xgb.DMatrix(X_train[sp$val,   , drop = FALSE], label = yb[sp$val])
      fit <- tryCatch(xgb.train(list(objective = "binary:logistic", eval_metric = "logloss",
                                     max_depth = 4, eta = 0.1, subsample = 0.8,
                                     colsample_bytree = 1, nthread = 1),
                                dtr, nrounds = 400, evals = list(val = dvl),
                                early_stopping_rounds = 20, verbose = 0),
                      error = function(e) NULL)
      if (is.null(fit)) next
      preds[, m] <- predict(fit, xgb.DMatrix(X_ext))
      pv <- predict(fit, dvl)
      r  <- tryCatch(roc(yb[sp$val], pv, quiet = TRUE), error = function(e) NULL)
      ths[m] <- if (is.null(r)) 0.5 else
        r$thresholds[which.max(r$sensitivities + r$specificities - 1)]
    }
    if (all(is.na(preds))) next

    p_ext <- rowMeans(preds, na.rm = TRUE)
    th <- median(ths[is.finite(ths)], na.rm = TRUE); if (!is.finite(th)) th <- 0.5

    auc <- as.numeric(roc(y_ext, p_ext, quiet = TRUE)$auc)
    pred_cls <- ifelse(p_ext > th, "cancer", "normal")
    cm <- table(True = factor(ext_class, levels = c("normal", "cancer")),
                Predicted = factor(pred_cls, levels = c("normal", "cancer")))
    sens <- cm["cancer", "cancer"] / sum(cm["cancer", ])
    spec <- cm["normal", "normal"] / sum(cm["normal", ])

    key <- paste(pn, arm, sep = "_")
    res_rows[[key]] <- data.frame(
      panel = pn, arm = arm, arm_label = ARM_LABEL[[arm]],
      n_panel_genes = length(PANELS[[pn]]), n_transferred = length(genes),
      n_external = length(y_ext), auc = round(auc, 4),
      sensitivity = round(sens, 3), specificity = round(spec, 3), threshold = round(th, 3))
    cm_rows[[key]] <- data.frame(panel = pn, arm = arm, as.data.frame(cm))

    df <- data.frame(score = p_ext, class = ext_class)
    p <- ggplot(df, aes(class, score, fill = class)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.5) +
      geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
      geom_hline(yintercept = th, linetype = "dashed", color = "red") +
      labs(title = paste0("External validation on ", EXT_GSE, " (Arm ", arm, ")"),
           subtitle = sprintf("%s | %s panel (%d genes) | AUC %.3f",
                              ARM_LABEL[[arm]], pn, length(genes), auc),
           x = NULL, y = "Predicted cancer score") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none", plot.title = element_text(face = "plain"))
    ggsave(file.path(fig_dir, sprintf("external_%s_%s.pdf", pn, arm)), p, width = 6.5, height = 5)
    ggsave(file.path(fig_dir, sprintf("external_%s_%s.png", pn, arm)), p, width = 6.5, height = 5, dpi = 300)
  }
}

if (length(res_rows)) {
  out <- do.call(rbind, res_rows)
  out <- out[order(out$panel, -out$auc), ]
  write.csv(out, file.path(res_dir, "external_validation.csv"), row.names = FALSE)
  write.csv(do.call(rbind, cm_rows), file.path(res_dir, "external_confusion.csv"), row.names = FALSE)
}