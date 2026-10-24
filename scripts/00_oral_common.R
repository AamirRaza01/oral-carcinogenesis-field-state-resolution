# ==============================================================================
# Script: 00_oral_common.R
# Purpose: Shared project configuration, metadata loading, and sample cohort
#          definitions across all downstream analytical pipelines.
# ==============================================================================

base_dir <- Sys.getenv("ORAL_CANCER_BASE_DIR", unset = ".")
out_dir  <- file.path(base_dir, "processed_data")
fig_dir  <- file.path(base_dir, "figures")
res_dir  <- file.path(base_dir, "results")
for (d in c(out_dir, fig_dir, res_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Master sample metadata for the pooled compendium
manifest_path <- file.path(base_dir, "Oral-Cancer-Samples-Metadata_corrected.csv")

ARMS <- c(A = "arm_A_collapsed_genes.rds",
          B = "arm_B_collapsed_genes.rds",
          C = "arm_C_collapsed_genes.rds")

ARM_LABEL <- c(A = "fRMA baseline",
               B = "ComBat mod = NULL",
               C = "ComBat mod = ~class")

# ------------------------------------------------------------------------------
# Sample Cohort Definitions
#   FULL  : Carcinoma vs. all available non-malignant mucosa.
#   DONOR : Carcinoma vs. independent healthy donor mucosa only.
# Dysplasia samples are held out from training and evaluated downstream.
# ------------------------------------------------------------------------------
DONOR_PROVENANCE <- "donor"
LADDER_ORDER <- c("donor", "contralateral", "adjacent", "matched", "dysplasia", "cancer")

load_manifest <- function() {
  if (!file.exists(manifest_path)) stop("Manifest not found: ", manifest_path)
  m <- read.csv(manifest_path, stringsAsFactors = FALSE)
  names(m) <- tolower(names(m))
  stopifnot(all(c("sample", "study", "class") %in% names(m)))
  if (!"normal_provenance" %in% names(m)) m$normal_provenance <- NA_character_
  m$normal_provenance[m$normal_provenance %in% c("NA", "")] <- NA_character_
  m$class <- tolower(trimws(m$class))
  if (!"batch" %in% names(m)) m$batch <- match(m$study, unique(m$study))

  # Assign ordinal progression tiers for the ladder analysis;
  # 'matched-unresolved' is grouped into 'matched' controls.
  m$tier <- ifelse(m$class == "cancer", "cancer",
            ifelse(m$class == "dysplasia", "dysplasia",
            ifelse(m$normal_provenance %in% "donor", "donor",
            ifelse(m$normal_provenance %in% "contralateral", "contralateral",
            ifelse(m$normal_provenance %in% "adjacent", "adjacent", "matched")))))
  m$tier <- factor(m$tier, levels = LADDER_ORDER)
  m
}

set_samples <- function(manifest, set = c("FULL", "DONOR")) {
  set <- match.arg(set)
  is_cancer <- manifest$class == "cancer"
  is_normal <- manifest$class == "normal"
  keep <- if (set == "FULL") is_cancer | is_normal
          else is_cancer | (is_normal & manifest$normal_provenance %in% DONOR_PROVENANCE)
  manifest$sample[keep]
}

set_labels <- function(manifest, samples) {
  cl <- manifest$class[match(samples, manifest$sample)]
  factor(ifelse(cl == "cancer", "cancer", "normal"), levels = c("normal", "cancer"))
}

load_arm <- function(arm) {
  f <- file.path(out_dir, ARMS[[arm]])
  if (!file.exists(f)) stop("Arm matrix not found: ", f)
  readRDS(f)
}

run_grid <- function() {
  g <- expand.grid(arm = names(ARMS), set = c("FULL", "DONOR"),
                   stringsAsFactors = FALSE)
  g$run <- paste0(g$arm, "_", g$set)
  g[order(g$set, g$arm), ]
}

# Stratified train / validation / test partition (64 / 16 / 20 split within class)
stratified_split <- function(y, train_frac = 0.64, val_frac = 0.16) {
  tr <- va <- te <- integer(0)
  for (lev in levels(y)) {
    idx <- which(y == lev); shuf <- sample(idx); n <- length(idx)
    n_tr <- max(1, round(n * train_frac)); n_va <- max(1, round(n * val_frac))
    if (n_tr + n_va >= n) n_va <- max(1, n - n_tr - 1)
    tr <- c(tr, shuf[1:n_tr])
    va <- c(va, shuf[(n_tr + 1):(n_tr + n_va)])
    te <- c(te, shuf[(n_tr + n_va + 1):n])
  }
  list(train = tr, val = va, test = te)
}

message("Loaded 00_oral_common.R | base_dir = ", base_dir)