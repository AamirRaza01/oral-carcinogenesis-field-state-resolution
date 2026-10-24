# ==============================================================================
# Script: 14_panel_grounding.R
# Purpose: Functional grounding of panel genes against somatic driver mutation
#          status (TCGA-HNSC oral cavity subset) and inferred transcription
#          factor activity via DoRothEA/VIPER regulons.
# Reference: Paper Section IV-H, Table IV
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

suppressPackageStartupMessages({ library(data.table) })

TCGA_DIR   <- file.path(base_dir, "tcga_hnsc")
PANEL_FILE <- file.path(res_dir, "continuum_multinomial_A.txt")
MIN_GROUP  <- 5
N_NULL     <- 1000
FDR        <- 0.05
HYPERMUT   <- 500

DRIVERS_PUBLISHED <- c("TP53", "CDKN2A", "PIK3CA", "NOTCH1", "CASP8", "HRAS",
                       "FAT1", "AJUBA", "NSD1", "NFE2L2", "TRAF3")

DRIVER_SOURCE    <- "union"
ONCODRIVE_FDR    <- 0.1
ONCODRIVE_MINMUT <- 5

read_lines_clean <- function(f) { x <- readLines(f); x[nzchar(x)] }

if (!file.exists(PANEL_FILE)) stop("Missing panel file: ", PANEL_FILE)
panel <- read_lines_clean(PANEL_FILE)

# ------------------------------------------------------------------------------
# 1. Load oral cavity cohort expression and phenotypic metadata
# ------------------------------------------------------------------------------
expr_file  <- file.path(TCGA_DIR, "hnsc_tmm_log2cpm_by_symbol.csv")
pheno_file <- file.path(TCGA_DIR, "hnsc_pheno.csv")
oral_file  <- file.path(TCGA_DIR, "hnsc_oral_cavity_samples.txt")
for (f in c(expr_file, pheno_file)) if (!file.exists(f)) stop("Missing ", f)

E <- as.matrix(fread(expr_file), rownames = 1)
ph <- fread(pheno_file)
oral <- if (file.exists(oral_file)) read_lines_clean(oral_file) else
  ph$sample_barcode[ph$hnsc_site == "Oral cavity"]

tum <- intersect(colnames(E), ph$sample_barcode[grepl("Tumor", ph$sample_type)])
tum <- intersect(tum, oral)
E <- E[, tum, drop = FALSE]

# ------------------------------------------------------------------------------
# 2. Somatic mutations (MC3 consensus calls)
# ------------------------------------------------------------------------------
maf_dt <- NULL
if (requireNamespace("TCGAmutations", quietly = TRUE)) {
  hnsc_maf <- TCGAmutations::tcga_load(study = "HNSC")
  maf_dt <- data.table::as.data.table(hnsc_maf@data)
} else {
  f <- Sys.glob(file.path(TCGA_DIR, "*maf*"))[1]
  if (is.na(f)) stop("No MAF file identified in ", TCGA_DIR)
  maf_dt <- fread(f)
}

maf <- maf_dt[Variant_Classification != "Silent"]
BARCODE_N <- 15
maf[, sample := substr(Tumor_Sample_Barcode, 1, BARCODE_N)]
colnames(E) <- substr(colnames(E), 1, BARCODE_N)
E <- E[, !duplicated(colnames(E)), drop = FALSE]
maf <- maf[sample %in% colnames(E)]
E <- E[, colnames(E) %in% unique(maf$sample), drop = FALSE]

tmb <- maf[, .N, by = sample]
hyper <- tmb$sample[tmb$N > HYPERMUT]

mut_status <- function(gene) colnames(E) %in% unique(maf$sample[maf$Hugo_Symbol == gene])

# ------------------------------------------------------------------------------
# 2b. Driver mutation identification via OncodriveCLUST
# ------------------------------------------------------------------------------
DRIVERS_ONCODRIVE <- character(0)
if (exists("hnsc_maf") && requireNamespace("maftools", quietly = TRUE)) {
  oral_tsb <- unique(maf_dt$Tumor_Sample_Barcode[
    substr(maf_dt$Tumor_Sample_Barcode, 1, BARCODE_N) %in% colnames(E)])
  oral_maf <- maftools::subsetMaf(hnsc_maf, tsb = oral_tsb, mafObj = TRUE)
  
  aacol <- intersect(c("HGVSp_Short", "Protein_Change", "amino_acid_change", "AAChange"),
                     names(maf_dt))[1]
  od <- tryCatch(maftools::oncodrive(maf = oral_maf, AACol = aacol,
                                     minMut = ONCODRIVE_MINMUT, pvalMethod = "zscore"),
                 error = function(e) NULL)
  
  if (!is.null(od)) {
    od <- data.table::as.data.table(od)
    write.csv(od, file.path(res_dir, "oncodrive_oral.csv"), row.names = FALSE)
    DRIVERS_ONCODRIVE <- od$Hugo_Symbol[od$fdr < ONCODRIVE_FDR]
    
    pdf(file.path(fig_dir, "oncodrive_oral.pdf"), width = 6, height = 5)
    try(maftools::plotOncodrive(res = od, fdrCutOff = ONCODRIVE_FDR, useFraction = TRUE), silent = TRUE)
    dev.off()
    
    png(file.path(fig_dir, "oncodrive_oral.png"), width = 6, height = 5, units = "in", res = 300)
    try(maftools::plotOncodrive(res = od, fdrCutOff = ONCODRIVE_FDR, useFraction = TRUE), silent = TRUE)
    dev.off()
  }
}

DRIVERS <- switch(DRIVER_SOURCE,
                  published = DRIVERS_PUBLISHED,
                  oncodrive = DRIVERS_ONCODRIVE,
                  union     = union(DRIVERS_PUBLISHED, DRIVERS_ONCODRIVE))
DRIVERS <- DRIVERS[DRIVERS %in% unique(maf$Hugo_Symbol)]

# ------------------------------------------------------------------------------
# 3. Mutation frequencies in panel and driver genes
# ------------------------------------------------------------------------------
freq <- data.frame(gene = panel,
                   n_mutated = vapply(panel, function(g) sum(mut_status(g)), integer(1)),
                   n_tumors = ncol(E))
freq$percent <- round(100 * freq$n_mutated / freq$n_tumors, 2)
write.csv(freq, file.path(res_dir, "panel_mutation_frequency.csv"), row.names = FALSE)

drv_freq <- data.frame(driver = DRIVERS,
                       n_mutated = vapply(DRIVERS, function(g) sum(mut_status(g)), integer(1)))
drv_freq$percent <- round(100 * drv_freq$n_mutated / ncol(E), 1)
drv_freq$testable <- drv_freq$n_mutated >= MIN_GROUP &
  (ncol(E) - drv_freq$n_mutated) >= MIN_GROUP
write.csv(drv_freq, file.path(res_dir, "driver_mutation_frequency.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Association testing: panel gene expression vs. driver mutation status
# ------------------------------------------------------------------------------
test_set <- function(genes, samples, label) {
  Es <- E[intersect(genes, rownames(E)), samples, drop = FALSE]
  rows <- list()
  for (d in DRIVERS) {
    mstat <- colnames(E) %in% unique(maf$sample[maf$Hugo_Symbol == d])
    names(mstat) <- colnames(E)
    m <- mstat[samples]
    if (sum(m) < MIN_GROUP || sum(!m) < MIN_GROUP) next
    for (g in rownames(Es)) {
      a <- Es[g, m]; b <- Es[g, !m]
      w <- suppressWarnings(wilcox.test(a, b))
      rows[[length(rows) + 1]] <- data.frame(
        stratum = label, gene = g, driver = d,
        n_mutant = sum(m), n_wildtype = sum(!m),
        log2fc = round(mean(a) - mean(b), 4), p_value = w$p.value)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out
}

strata <- list(all = colnames(E), non_hypermutated = setdiff(colnames(E), hyper))
pair_rows <- list(); driver_summary_rows <- list(); unified_summary_rows <- list()
panel_in <- intersect(panel, rownames(E))

for (nm in names(strata)) {
  s <- strata[[nm]]
  if (length(s) < 2 * MIN_GROUP) next
  obs <- test_set(panel_in, s, nm)
  if (is.null(obs)) next
  pair_rows[[nm]] <- obs
  sig <- obs[obs$p_adj < FDR, ]
  n_genes_resp <- length(unique(sig$gene))

  driver_summary_rows[[nm]] <- data.frame(
    stratum = nm, n_tumors = length(s), n_pairs = nrow(obs),
    n_significant = nrow(sig), n_genes_responding = n_genes_resp,
    n_genes_tested = length(unique(obs$gene)),
    median_abs_log2fc = round(median(abs(sig$log2fc)), 3),
    max_abs_log2fc = round(max(abs(sig$log2fc)), 3))

  unified_summary_rows[[paste0("driver_", nm)]] <- data.frame(
    layer = "driver", stratum = nm, n_tumors = length(s), n_pairs = nrow(obs),
    n_significant = nrow(sig), n_genes_responding = n_genes_resp,
    n_genes_tested = length(unique(obs$gene)),
    metric_type = "log2FC",
    median_effect = round(median(abs(sig$log2fc)), 3),
    max_effect = round(max(abs(sig$log2fc)), 3))
}

if (length(pair_rows)) {
  write.csv(do.call(rbind, pair_rows), file.path(res_dir, "panel_driver_pairs.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 5. Transcription factor regulon analysis (DoRothEA/VIPER)
# ------------------------------------------------------------------------------
tf_ok <- requireNamespace("viper", quietly = TRUE) && requireNamespace("dorothea", quietly = TRUE)
if (tf_ok) {
  reg <- dorothea::dorothea_hs
  reg <- reg[reg$confidence %in% c("A", "B", "C"), ]
  regulons <- split(reg, reg$tf)
  regulons <- regulons[vapply(regulons, nrow, integer(1)) >= 20]
  vip <- viper::viper(E, lapply(regulons, function(r)
    list(tfmode = setNames(r$mor, r$target), likelihood = rep(1, nrow(r)))),
    method = "scale", minsize = 20, verbose = FALSE)

  tf_test <- function(genes) {
    rows <- list()
    for (g in intersect(genes, rownames(E))) for (tf in rownames(vip)) {
      ct <- suppressWarnings(cor.test(E[g, ], vip[tf, ], method = "spearman"))
      rows[[length(rows) + 1]] <- data.frame(gene = g, tf = tf,
                                             rho = round(unname(ct$estimate), 4),
                                             p_value = ct$p.value)
    }
    out <- do.call(rbind, rows)
    out$p_adj <- p.adjust(out$p_value, method = "BH")
    out
  }

  obs <- tf_test(panel_in)
  sig <- obs[obs$p_adj < FDR, ]
  strong <- sig[abs(sig$rho) >= 0.4, ]
  write.csv(obs, file.path(res_dir, "panel_tf_pairs.csv"), row.names = FALSE)

  reg_summary <- data.frame(
    stratum = "oral tumors", n_tumors = ncol(E), n_pairs = nrow(obs),
    n_significant = nrow(sig), n_strong = nrow(strong),
    n_genes_responding = length(unique(sig$gene)),
    n_genes_tested = length(unique(obs$gene)),
    median_abs_rho = round(median(abs(sig$rho)), 3),
    max_abs_rho = round(max(abs(sig$rho)), 3))
  write.csv(reg_summary, file.path(res_dir, "panel_tf_summary.csv"), row.names = FALSE)

  unified_summary_rows[["regulator"]] <- data.frame(
    layer = "regulator", stratum = "oral tumors", n_tumors = ncol(E), n_pairs = nrow(obs),
    n_significant = nrow(sig), n_genes_responding = length(unique(sig$gene)),
    n_genes_tested = length(unique(obs$gene)),
    metric_type = "Spearman_rho",
    median_effect = round(median(abs(sig$rho)), 3),
    max_effect = round(max(abs(sig$rho)), 3))

  if (nrow(strong)) {
    top <- as.data.frame(table(strong$tf), stringsAsFactors = FALSE)
    names(top) <- c("tf", "n")
    top <- top[order(-top$n), ]
    write.csv(top, file.path(res_dir, "panel_tf_top_regulators.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 6. Export summary reports
# ------------------------------------------------------------------------------
if (length(driver_summary_rows)) {
  write.csv(do.call(rbind, driver_summary_rows), file.path(res_dir, "panel_driver_summary.csv"), row.names = FALSE)
}
if (length(unified_summary_rows)) {
  write.csv(do.call(rbind, unified_summary_rows), file.path(res_dir, "panel_grounding_summary.csv"), row.names = FALSE)
}