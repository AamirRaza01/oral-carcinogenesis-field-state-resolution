# ==============================================================================
# Script: 03_build_preprocessing_arms.R
# Purpose: Generate preprocessing Arms A, B, and C:
#          - Arm A: Unadjusted fRMA baseline.
#          - Arm B: ComBat adjusted with mod = NULL (unguided batch removal).
#          - Arm C: ComBat adjusted with mod = ~class (protecting carcinoma vs.
#                   non-malignant contrast).
#          Applies a 75% probe variance filter and collapses to unique gene symbols.
# Reference: Paper Section III-B (Preprocessing Arms)
# ==============================================================================

source("00_oral_common.R")

library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

COMBAT_MOD_SCHEME <- "two_level"   # "two_level" or "three_level"
VAR_KEEP_FRACTION <- 0.75          # retain top 75% most variable probes

# ------------------------------------------------------------------------------
# 1. Load manifest and pooled fRMA matrix
# ------------------------------------------------------------------------------
manifest <- load_manifest()

master_file <- file.path(out_dir, "master_frma_pooled_463.rds")
master_expr <- readRDS(master_file)
master_expr <- master_expr[, manifest$sample, drop = FALSE]
message(sprintf("Master matrix: %d probes x %d arrays", nrow(master_expr), ncol(master_expr)))

batch_factor <- as.factor(manifest$study)
class_factor <- as.factor(manifest$class)

# Binary biological contrast vector (carcinoma vs. non-carcinoma)
class_binary <- factor(ifelse(manifest$class == "cancer", "cancer", "normal"),
                       levels = c("normal", "cancer"))

tab <- table(manifest$study, manifest$class)
single <- rownames(tab)[rowSums(tab > 0) < 2]
if (length(single))
  warning("Series contributing a single class: ", paste(single, collapse = ", "))

# ------------------------------------------------------------------------------
# 2. Build preprocessing arms
# ------------------------------------------------------------------------------
message("\n--> Arm A: fRMA baseline")
expr_arm_A_raw <- master_expr

message("--> Arm B: ComBat (mod = NULL)")
expr_arm_B_raw <- ComBat(dat = master_expr, batch = batch_factor, mod = NULL,
                         par.prior = TRUE, prior.plots = FALSE)

message("--> Arm C: ComBat (mod = ~class, scheme = ", COMBAT_MOD_SCHEME, ")")
design_class <- if (COMBAT_MOD_SCHEME == "two_level")
  model.matrix(~ class_binary) else model.matrix(~ class_factor)

expr_arm_C_raw <- ComBat(dat = master_expr, batch = batch_factor, mod = design_class,
                         par.prior = TRUE, prior.plots = FALSE)

saveRDS(expr_arm_A_raw, file.path(out_dir, "arm_A_probes_raw.rds"))
saveRDS(expr_arm_B_raw, file.path(out_dir, "arm_B_probes_raw.rds"))
saveRDS(expr_arm_C_raw, file.path(out_dir, "arm_C_probes_raw.rds"))

# ------------------------------------------------------------------------------
# 3. Define shared gene universe from uncorrected baseline (Arm A)
# ------------------------------------------------------------------------------
# Variance filtering and probe selection are evaluated once on Arm A and applied
# identically to all arms to maintain an invariant feature space across arms.
probe_vars <- apply(expr_arm_A_raw, 1, var)
var_threshold <- quantile(probe_vars, probs = 1 - VAR_KEEP_FRACTION)
keep_probes <- names(probe_vars)[probe_vars >= var_threshold]

symbol_map <- mapIds(hgu133plus2.db, keys = keep_probes, column = "SYMBOL",
                     keytype = "PROBEID", multiVals = "first")
ok <- !is.na(symbol_map) & symbol_map != ""
keep_probes <- keep_probes[ok]
keep_symbols <- symbol_map[ok]

probe_means <- rowMeans(expr_arm_A_raw[keep_probes, , drop = FALSE])
split_idx <- split(seq_along(keep_probes), keep_symbols)
selected <- vapply(split_idx, function(ix)
  if (length(ix) == 1L) ix else ix[which.max(probe_means[ix])], integer(1))

SHARED_PROBES  <- keep_probes[selected]
SHARED_SYMBOLS <- names(split_idx)
ord <- order(SHARED_SYMBOLS)
SHARED_PROBES  <- SHARED_PROBES[ord]
SHARED_SYMBOLS <- SHARED_SYMBOLS[ord]

write.csv(data.frame(gene_symbol = SHARED_SYMBOLS, probe_id = SHARED_PROBES,
                     mean_intensity_armA = probe_means[SHARED_PROBES]),
          file.path(out_dir, "shared_symbol_selection_map.csv"), row.names = FALSE)

collapse_arm <- function(expr_mat, arm_name) {
  out <- expr_mat[SHARED_PROBES, , drop = FALSE]
  rownames(out) <- SHARED_SYMBOLS
  message(sprintf("  %s: %d genes x %d arrays", arm_name, nrow(out), ncol(out)))
  out
}

# ------------------------------------------------------------------------------
# 4. Collapse onto shared gene universe and export
# ------------------------------------------------------------------------------
final_arm_A <- collapse_arm(expr_arm_A_raw, "Arm A (fRMA baseline)")
final_arm_B <- collapse_arm(expr_arm_B_raw, "Arm B (ComBat mod = NULL)")
final_arm_C <- collapse_arm(expr_arm_C_raw, "Arm C (ComBat mod = ~class)")

stopifnot(identical(rownames(final_arm_A), rownames(final_arm_B)),
          identical(rownames(final_arm_A), rownames(final_arm_C)),
          identical(colnames(final_arm_A), manifest$sample),
          identical(colnames(final_arm_B), manifest$sample),
          identical(colnames(final_arm_C), manifest$sample))

saveRDS(final_arm_A, file.path(out_dir, "arm_A_collapsed_genes.rds"))
saveRDS(final_arm_B, file.path(out_dir, "arm_B_collapsed_genes.rds"))
saveRDS(final_arm_C, file.path(out_dir, "arm_C_collapsed_genes.rds"))

sets <- lapply(c("FULL", "DONOR"), function(s) set_samples(manifest, s))
names(sets) <- c("FULL", "DONOR")
for (s in names(sets))
  writeLines(sets[[s]], file.path(out_dir, sprintf("samples_%s.txt", s)))