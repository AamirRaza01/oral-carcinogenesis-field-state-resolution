# ==============================================================================
# Script: 3.1_build_preprocessing_arms.R
# Purpose: Generate preprocessing Arms A, B, and C with ComBat, apply a 75%
#          variance filter, and collapse probe sets to unique gene symbols.
# ==============================================================================

base_dir <- Sys.getenv("ORAL_CANCER_BASE_DIR", unset = ".")
out_dir  <- file.path(base_dir, "processed_data")

library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

# ------------------------------------------------------------------------------
# 1. Load manifest and master pooled fRMA matrix
# ------------------------------------------------------------------------------
manifest_path <- file.path(base_dir, "Oral-Cancer-Samples-Metadata_corrected.csv")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

master_file <- file.path(out_dir, "master_frma_pooled_463.rds")
master_expr <- readRDS(master_file)

# Align sample columns with manifest order
master_expr <- master_expr[, manifest$sample, drop = FALSE]

batch_factor <- as.factor(manifest$study)
class_factor <- as.factor(manifest$class)

# ------------------------------------------------------------------------------
# 2. Construct raw expression matrices for Arms A, B, and C
# ------------------------------------------------------------------------------
message("\n--> Generating Arm A: fRMA baseline...")
expr_arm_A_raw <- master_expr

message("\n--> Generating Arm B: ComBat (mod = NULL)...")
expr_arm_B_raw <- ComBat(
  dat   = master_expr,
  batch = batch_factor,
  mod   = NULL,
  par.prior = TRUE,
  prior.plots = FALSE
)

message("\n--> Generating Arm C: ComBat (mod = ~class)...")
design_class <- model.matrix(~class_factor)
expr_arm_C_raw <- ComBat(
  dat   = master_expr,
  batch = batch_factor,
  mod   = design_class,
  par.prior = TRUE,
  prior.plots = FALSE
)

saveRDS(expr_arm_A_raw, file.path(out_dir, "arm_A_probes_raw.rds"))
saveRDS(expr_arm_B_raw, file.path(out_dir, "arm_B_probes_raw.rds"))
saveRDS(expr_arm_C_raw, file.path(out_dir, "arm_C_probes_raw.rds"))

# ------------------------------------------------------------------------------
# 3. Variance filtering and probe-to-symbol collapse
# ------------------------------------------------------------------------------
probe_ids <- rownames(master_expr)
symbol_map <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

process_arm <- function(expr_mat, arm_name) {
  message(sprintf("\nProcessing %s: Variance filtering and symbol collapsing...", arm_name))
  
  # Retain top 75% most variable probes
  probe_vars <- apply(expr_mat, 1, var)
  var_threshold <- quantile(probe_vars, probs = 0.25)
  keep_probes <- probe_vars >= var_threshold
  expr_filtered <- expr_mat[keep_probes, , drop = FALSE]
  
  # Remove unannotated probes
  arm_symbols <- symbol_map[rownames(expr_filtered)]
  valid_mask <- !is.na(arm_symbols) & (arm_symbols != "")
  expr_annotated <- expr_filtered[valid_mask, , drop = FALSE]
  matched_symbols <- arm_symbols[valid_mask]
  
  # Collapse probe sets by highest mean expression across all samples
  probe_means <- rowMeans(expr_annotated)
  split_indices <- split(seq_len(nrow(expr_annotated)), matched_symbols)
  selected_indices <- sapply(split_indices, function(idx) {
    if (length(idx) == 1) return(idx)
    idx[which.max(probe_means[idx])]
  })
  
  expr_collapsed <- expr_annotated[selected_indices, , drop = FALSE]
  rownames(expr_collapsed) <- names(split_indices)
  expr_collapsed <- expr_collapsed[order(rownames(expr_collapsed)), , drop = FALSE]
  
  return(expr_collapsed)
}

final_arm_A <- process_arm(expr_arm_A_raw, "Arm A (fRMA alone)")
final_arm_B <- process_arm(expr_arm_B_raw, "Arm B (fRMA + ComBat mod=NULL)")
final_arm_C <- process_arm(expr_arm_C_raw, "Arm C (fRMA + ComBat mod=class)")

saveRDS(final_arm_A, file.path(out_dir, "arm_A_collapsed_genes.rds"))
saveRDS(final_arm_B, file.path(out_dir, "arm_B_collapsed_genes.rds"))
saveRDS(final_arm_C, file.path(out_dir, "arm_C_collapsed_genes.rds"))

message("\nGenerated collapsed expression matrices in ", out_dir)