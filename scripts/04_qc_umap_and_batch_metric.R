# ==============================================================================
# Script: 04_qc_umap_and_batch_metric.R
# Purpose: Quantify inter-series mixing and phenotype preservation via average
#          silhouette width in principal component space, and generate comparative
#          UMAP projections across preprocessing arms.
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

suppressPackageStartupMessages({
  library(uwot); library(ggplot2); library(patchwork); library(cluster)
})

N_NEIGHBORS <- 20
MIN_DIST    <- 0.1
N_PC        <- 20
USE_PROBES  <- TRUE

PROBE_FILES <- c(A = "arm_A_probes_raw.rds",
                 B = "arm_B_probes_raw.rds",
                 C = "arm_C_probes_raw.rds")

manifest <- load_manifest()

load_matrix <- function(arm) {
  f <- file.path(out_dir, if (USE_PROBES) PROBE_FILES[[arm]] else ARMS[[arm]])
  if (!file.exists(f)) stop("Missing ", f)
  m <- readRDS(f)
  m[, manifest$sample, drop = FALSE]
}

# ------------------------------------------------------------------------------
# 1. Silhouette metrics and UMAP embeddings per arm
# ------------------------------------------------------------------------------
umap_rows <- list(); sil_rows <- list()

for (arm in names(ARMS)) {
  message("\nProcessing Arm ", arm, " (", ARM_LABEL[[arm]], ")")
  m <- load_matrix(arm)
  X <- t(m)
  rm(m); invisible(gc())

  pc <- prcomp(X, center = TRUE, scale. = FALSE, rank. = N_PC)
  var_pc <- round(100 * sum(pc$sdev[seq_len(N_PC)]^2) / sum(pc$sdev^2), 1)

  d <- dist(pc$x[, seq_len(N_PC)])
  for (lab in c("study", "class")) {
    g <- as.integer(as.factor(manifest[[lab]]))
    s <- summary(silhouette(g, d))$avg.width
    sil_rows[[paste(arm, lab)]] <- data.frame(
      arm = arm, arm_label = ARM_LABEL[[arm]], label = lab,
      n_groups = length(unique(g)), avg_silhouette = round(s, 4),
      pc_variance_percent = var_pc)
  }
  rm(d); invisible(gc())

  u <- uwot::umap(X, n_neighbors = N_NEIGHBORS, min_dist = MIN_DIST,
                  metric = "euclidean", n_threads = 2)
  rm(X); invisible(gc())

  umap_rows[[arm]] <- data.frame(UMAP1 = u[, 1], UMAP2 = u[, 2],
                                 study = manifest$study, class = manifest$class,
                                 arm = ARM_LABEL[[arm]])
}

sil <- do.call(rbind, sil_rows)
write.csv(sil, file.path(res_dir, "batch_silhouette.csv"), row.names = FALSE)

emb <- do.call(rbind, umap_rows)
emb$arm <- factor(emb$arm, levels = unname(ARM_LABEL[names(ARMS)]))
write.csv(emb, file.path(res_dir, "umap_coordinates.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Publication Figure 1 Generation
# ------------------------------------------------------------------------------
theme_paper <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "#0b0b0b", linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(size = base, colour = "#0b0b0b"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(size = base + 1, colour = "#0b0b0b"))
}

# Color-blind safe palette (Okabe and Ito)
study_colors <- c("#0072b2", "#e69f00", "#009e73", "#cc79a7",
                  "#56b4e9", "#d55e00", "#f0e442")
names(study_colors) <- sort(unique(emb$study))
class_colors <- c(cancer = "#eb6834", normal = "#1baf7a", dysplasia = "#2a78d6")

lab_study <- sprintf("(a) series label, average silhouette width %s",
                     paste(sprintf("%.2f", sil$avg_silhouette[sil$label == "study"]),
                           collapse = ", "))
lab_class <- sprintf("(b) tissue class, average silhouette width %s",
                     paste(sprintf("%.2f", sil$avg_silhouette[sil$label == "class"]),
                           collapse = ", "))

p_study <- ggplot(emb, aes(UMAP1, UMAP2, colour = study)) +
  geom_point(size = 1.1, alpha = 0.8) +
  facet_wrap(~arm, nrow = 1, scales = "free") +
  scale_colour_manual(values = study_colors) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.6, alpha = 1))) +
  labs(title = lab_study) + theme_paper()

p_class <- ggplot(emb, aes(UMAP1, UMAP2, colour = class)) +
  geom_point(size = 1.1, alpha = 0.8) +
  facet_wrap(~arm, nrow = 1, scales = "free") +
  scale_colour_manual(values = class_colors) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.6, alpha = 1))) +
  labs(title = lab_class) + theme_paper()

fig <- p_study / p_class
ggsave(file.path(fig_dir, "Fig1_UMAP_Arms.pdf"), fig, width = 10, height = 6.4)
ggsave(file.path(fig_dir, "Fig1_UMAP_Arms.png"), fig, width = 10, height = 6.4, dpi = 300)