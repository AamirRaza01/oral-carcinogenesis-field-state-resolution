# ==============================================================================
# Script: 01_frma_normalization.R
# Purpose: Isolated per-series fRMA normalization and master compendium assembly.
# ==============================================================================

options(timeout = 3600)
options(download.file.method = "libcurl")

base_dir <- Sys.getenv("ORAL_CANCER_BASE_DIR", unset = ".")
temp_dir <- file.path(base_dir, "temp")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = temp_dir, TEMP = temp_dir, TMP = temp_dir)

library(GEOquery)
library(affy)
library(frma)
library(hgu133plus2frmavecs)

raw_dir  <- file.path(base_dir, "raw_data")
out_dir  <- file.path(base_dir, "processed_data")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(base_dir, "Oral-Cancer-Samples-Metadata_corrected.csv")
if (!file.exists(manifest_path)) {
  stop("Manifest not found: ", manifest_path)
}
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

studies_ordered <- c("GSE74530", "GSE138206", "GSE78060", "GSE6791", "GSE9844", "GSE42743", "GSE30784")

# Chunk size for memory-managed CEL file loading
CHUNK_SIZE <- 25

# ------------------------------------------------------------------------------
# 1. Process each series independently
# ------------------------------------------------------------------------------
for (geo_id in studies_ordered) {
  out_rds <- file.path(out_dir, paste0(geo_id, "_frma.rds"))
  
  if (file.exists(out_rds)) {
    message("--> [SKIP] Found existing normalized file for: ", geo_id)
    next
  }
  
  message("\nProcessing GEO series: ", geo_id)
  series_dir <- file.path(raw_dir, geo_id)
  dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)
  
  tar_files <- list.files(series_dir, pattern = "\\.tar$", full.names = TRUE)
  if (length(tar_files) == 0) {
    message("Downloading raw supplementary TAR for ", geo_id, "...")
    getGEOSuppFiles(geo_id, baseDir = raw_dir, makeDirectory = TRUE, filter_regex = "\\.tar$")
    tar_files <- list.files(series_dir, pattern = "\\.tar$", full.names = TRUE)
  }
  
  if (length(tar_files) > 0) {
    message("Extracting CEL archives...")
    untar(tar_files[1], exdir = series_dir)
    unlink(tar_files[1])
  }
  
  target_samples <- manifest$sample[manifest$study == geo_id]
  all_cels <- list.files(series_dir, pattern = "\\.CEL\\.gz$|\\.CEL$", ignore.case = TRUE, full.names = TRUE)
  
  matched_cels <- character(0)
  for (gsm in target_samples) {
    hit <- grep(gsm, all_cels, value = TRUE)
    if (length(hit) > 0) matched_cels <- c(matched_cels, hit[1])
  }
  
  message(sprintf("Matched %d CEL files for %d target samples.", length(matched_cels), length(target_samples)))
  if (length(matched_cels) != length(target_samples)) {
    stop("Mismatch in sample count for ", geo_id)
  }
  
  num_chunks <- ceiling(length(matched_cels) / CHUNK_SIZE)
  chunk_list <- list()
  
  for (chk in seq_len(num_chunks)) {
    start_idx <- ((chk - 1) * CHUNK_SIZE) + 1
    end_idx   <- min(chk * CHUNK_SIZE, length(matched_cels))
    current_files <- matched_cels[start_idx:end_idx]
    
    affy_obj <- ReadAffy(filenames = current_files)
    sampleNames(affy_obj) <- sub("(GSM[0-9]+).*", "\\1", sampleNames(affy_obj))
    
    frma_obj <- frma(affy_obj, verbose = FALSE)
    chunk_list[[chk]] <- exprs(frma_obj)
    
    rm(affy_obj, frma_obj)
    gc()
  }
  
  series_mat <- do.call(cbind, chunk_list)
  series_mat <- series_mat[order(rownames(series_mat)), target_samples, drop = FALSE]
  
  saveRDS(series_mat, file = out_rds)
  message("Saved: ", out_rds, " [", nrow(series_mat), " probes x ", ncol(series_mat), " samples]")
  
  # Remove temporary uncompressed raw files
  unlink(series_dir, recursive = TRUE)
  rm(series_mat, chunk_list)
  gc()
}

# ------------------------------------------------------------------------------
# 2. Master Assembly: Build 463-sample pooled expression matrix
# ------------------------------------------------------------------------------
message("\nAssembling master pooled matrix across all 7 studies...")

mat_list <- list()
for (geo_id in manifest$study[!duplicated(manifest$study)]) {
  rds_file <- file.path(out_dir, paste0(geo_id, "_frma.rds"))
  mat_list[[geo_id]] <- readRDS(rds_file)
}

ref_probes <- rownames(mat_list[[1]])
for (geo_id in names(mat_list)) {
  if (!identical(rownames(mat_list[[geo_id]]), ref_probes)) {
    stop("Probe misalignment detected in: ", geo_id)
  }
}

master_expr <- do.call(cbind, mat_list)
master_expr <- master_expr[, manifest$sample, drop = FALSE]

master_file <- file.path(out_dir, "master_frma_pooled_463.rds")
saveRDS(master_expr, file = master_file)

unlink(temp_dir, recursive = TRUE)
message("Saved master matrix: ", master_file)