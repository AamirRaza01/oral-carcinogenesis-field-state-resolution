# ==============================================================================
# Script: 4.1_quality_control_figures.R
# Purpose: High-resolution standalone and combined UMAP projections alongside
#          per-array variance boxplots and paired Wilcoxon statistics.
# ==============================================================================

base_dir <- Sys.getenv("ORAL_CANCER_BASE_DIR", unset = ".")
out_dir  <- file.path(base_dir, "processed_data")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(uwot)
library(ggplot2)
library(patchwork)
library(scales)

# ------------------------------------------------------------------------------
# 1. Load data and verify metadata synchronization
# ------------------------------------------------------------------------------
manifest_path <- file.path(base_dir, "Oral-Cancer-Samples-Metadata_corrected.csv")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

arm_A <- readRDS(file.path(out_dir, "arm_A_collapsed_genes.rds"))
arm_B <- readRDS(file.path(out_dir, "arm_B_collapsed_genes.rds"))
arm_C <- readRDS(file.path(out_dir, "arm_C_collapsed_genes.rds"))

stopifnot(identical(colnames(arm_A), manifest$sample))
stopifnot(identical(colnames(arm_B), manifest$sample))
stopifnot(identical(colnames(arm_C), manifest$sample))

manifest$class <- factor(manifest$class, levels = c("cancer", "normal", "dysplasia"))
manifest$study <- factor(manifest$study)

# ------------------------------------------------------------------------------
# 2. Graphics styling definitions
# ------------------------------------------------------------------------------
theme_journal <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.title         = element_text(face = "plain", size = rel(1.15), color = "#222222", hjust = 0),
      plot.subtitle      = element_text(face = "plain", size = rel(0.95), color = "#555555", hjust = 0, margin = margin(b = 8)),
      axis.title         = element_text(face = "plain", size = rel(1.0), color = "#333333"),
      axis.text          = element_text(face = "plain", size = rel(0.9), color = "#444444"),
      strip.background   = element_rect(fill = "#EFEFEF", color = "#444444", linewidth = 0.6),
      strip.text         = element_text(face = "plain", size = rel(0.95), color = "#222222", margin = margin(4, 4, 4, 4)),
      panel.grid.major   = element_line(color = "#ECECEC", linewidth = 0.45),
      panel.grid.minor   = element_line(color = "#F5F5F5", linewidth = 0.25),
      panel.border       = element_rect(color = "#444444", fill = NA, linewidth = 0.6),
      legend.position    = "bottom",
      legend.title       = element_text(face = "plain", size = rel(0.9), color = "#333333"),
      legend.text        = element_text(face = "plain", size = rel(0.85), color = "#444444"),
      legend.key         = element_blank(),
      legend.background  = element_blank()
    )
}

study_palette <- c(
  "GSE30784"  = "#E64B35",
  "GSE42743"  = "#4DBBD5",
  "GSE9844"   = "#00A087",
  "GSE6791"   = "#3C5488",
  "GSE78060"  = "#F39B7F",
  "GSE138206" = "#8491B4",
  "GSE74530"  = "#91D1C2"
)

class_palette <- c(
  "cancer"    = "#E76F51",
  "normal"    = "#2A9D8F",
  "dysplasia" = "#7209B7"
)

# ------------------------------------------------------------------------------
# 3. Compute UMAP embeddings
# ------------------------------------------------------------------------------
set.seed(42)
calc_umap <- function(mat) {
  as.data.frame(uwot::umap(t(mat), n_neighbors = 15, min_dist = 0.3, metric = "cosine", n_threads = 2))
}

u_A <- calc_umap(arm_A)
u_B <- calc_umap(arm_B)
u_C <- calc_umap(arm_C)
colnames(u_A) <- colnames(u_B) <- colnames(u_C) <- c("UMAP1", "UMAP2")

df_A <- cbind(u_A, manifest, arm = "Arm A: fRMA (Unadjusted)")
df_B <- cbind(u_B, manifest, arm = "Arm B: ComBat (mod = NULL)")
df_C <- cbind(u_C, manifest, arm = "Arm C: ComBat (mod = ~class)")

combined_umap_df <- rbind(df_A, df_B, df_C)
combined_umap_df$arm <- factor(combined_umap_df$arm, levels = c(
  "Arm A: fRMA (Unadjusted)",
  "Arm B: ComBat (mod = NULL)",
  "Arm C: ComBat (mod = ~class)"
))

# ------------------------------------------------------------------------------
# 4. Generate standalone UMAP panels
# ------------------------------------------------------------------------------
arms_list <- list(
  list(df = df_A, prefix = "ArmA_fRMA", title_batch = "Arm A (fRMA) - Study Batch Structure", title_bio = "Arm A (fRMA) - Phenotype Separation"),
  list(df = df_B, prefix = "ArmB_ComBat_NULL", title_batch = "Arm B (ComBat mod=NULL) - Study Batch Structure", title_bio = "Arm B (ComBat mod=NULL) - Phenotype Separation"),
  list(df = df_C, prefix = "ArmC_ComBat_class", title_batch = "Arm C (ComBat mod=~class) - Study Batch Structure", title_bio = "Arm C (ComBat mod=~class) - Phenotype Separation")
)

for (item in arms_list) {
  p_s <- ggplot(item$df, aes(x = UMAP1, y = UMAP2, color = study)) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_manual(values = study_palette, name = "GEO Series") +
    guides(color = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
    labs(title = item$title_batch, x = "UMAP 1", y = "UMAP 2") +
    theme_journal()
  
  p_c <- ggplot(item$df, aes(x = UMAP1, y = UMAP2, color = class)) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_manual(values = class_palette, name = "Diagnostic Class", labels = c("Cancer", "Normal", "Dysplasia")) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1))) +
    labs(title = item$title_bio, x = "UMAP 1", y = "UMAP 2") +
    theme_journal()
  
  p_arm_both <- (p_s | p_c) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  ggsave(file.path(fig_dir, paste0("Fig_UMAP_", item$prefix, "_DualPanel.pdf")), p_arm_both, width = 10.5, height = 5.2)
  ggsave(file.path(fig_dir, paste0("Fig_UMAP_", item$prefix, "_DualPanel.png")), p_arm_both, width = 10.5, height = 5.2, dpi = 300)
}

# ------------------------------------------------------------------------------
# 5. Multi-panel comparative overview
# ------------------------------------------------------------------------------
p_all_study <- ggplot(combined_umap_df, aes(x = UMAP1, y = UMAP2, color = study)) +
  geom_point(alpha = 0.75, size = 1.4) +
  facet_wrap(~arm, nrow = 1) +
  scale_color_manual(values = study_palette, name = "GEO Series") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1))) +
  labs(title = "A  Cross-Study Batch Integration Across Preprocessing Arms", x = "UMAP 1", y = "UMAP 2") +
  theme_journal()

p_all_class <- ggplot(combined_umap_df, aes(x = UMAP1, y = UMAP2, color = class)) +
  geom_point(alpha = 0.75, size = 1.4) +
  facet_wrap(~arm, nrow = 1) +
  scale_color_manual(values = class_palette, name = "Diagnostic Class", labels = c("Cancer", "Normal", "Dysplasia")) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1))) +
  labs(title = "B  Biological Phenotype Preservation Across Preprocessing Arms", x = "UMAP 1", y = "UMAP 2") +
  theme_journal()

fig_umap_master <- (p_all_study / p_all_class) + plot_layout(heights = c(1, 1))

ggsave(file.path(fig_dir, "Fig1_UMAP_Arms_Comparative_Master.pdf"), fig_umap_master, width = 12.0, height = 8.5)
ggsave(file.path(fig_dir, "Fig1_UMAP_Arms_Comparative_Master.png"), fig_umap_master, width = 12.0, height = 8.5, dpi = 300)

# ------------------------------------------------------------------------------
# 6. Array variance analysis and paired Wilcoxon signed-rank tests
# ------------------------------------------------------------------------------
var_A <- apply(arm_A, 2, var)
var_B <- apply(arm_B, 2, var)
var_C <- apply(arm_C, 2, var)

test_AB <- wilcox.test(var_A, var_B, paired = TRUE)
test_AC <- wilcox.test(var_A, var_C, paired = TRUE)
test_BC <- wilcox.test(var_B, var_C, paired = TRUE)

format_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2 x 10^-16")
  return(sprintf("p = %.2e", p))
}

stats_caption <- paste0(
  "Paired Wilcoxon signed-rank test (N = 463 paired arrays):\n",
  "Arm A vs. Arm B: ", format_p(test_AB$p.value),
  "  |  Arm A vs. Arm C: ", format_p(test_AC$p.value),
  "  |  Arm B vs. Arm C: ", format_p(test_BC$p.value)
)

var_df <- data.frame(
  sample   = rep(manifest$sample, 3),
  study    = rep(manifest$study, 3),
  class    = rep(manifest$class, 3),
  variance = c(var_A, var_B, var_C),
  arm      = factor(
    rep(c("Arm A: fRMA", "Arm B: ComBat NULL", "Arm C: ComBat ~class"), each = nrow(manifest)),
    levels = c("Arm A: fRMA", "Arm B: ComBat NULL", "Arm C: ComBat ~class")
  )
)

arm_colors <- c(
  "Arm A: fRMA"             = "#E76F51",
  "Arm B: ComBat NULL"      = "#2A9D8F",
  "Arm C: ComBat ~class"    = "#457B9D"
)

# Global array variance boxplot
p_var_global <- ggplot(var_df, aes(x = arm, y = variance, fill = arm)) +
  geom_boxplot(width = 0.42, outlier.shape = NA, color = "#111111", linewidth = 0.55, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.18, seed = 42), size = 0.75, alpha = 0.25, color = "#111111") +
  scale_fill_manual(values = arm_colors) +
  labs(title = "Per-Array Expression Variance Across Preprocessing Arms",
       subtitle = stats_caption, x = "Preprocessing Arm", y = "Array Expression Variance") +
  theme_journal(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Fig2_Array_Variance_Global_Wilcoxon.pdf"), p_var_global, width = 7.0, height = 5.5)
ggsave(file.path(fig_dir, "Fig2_Array_Variance_Global_Wilcoxon.png"), p_var_global, width = 7.0, height = 5.5, dpi = 300)

# Study-stratified array variance distributions
p_var_by_study <- ggplot(var_df, aes(x = arm, y = variance, fill = arm)) +
  geom_boxplot(width = 0.50, outlier.shape = NA, color = "#111111", linewidth = 0.5, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.18, seed = 42), size = 0.65, alpha = 0.28, color = "#111111") +
  facet_wrap(~study, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = arm_colors) +
  labs(title = "Per-Array Variance Distribution Stratified by GEO Series",
       x = "Preprocessing Arm", y = "Array Expression Variance") +
  theme_journal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(fig_dir, "Fig2_Array_Variance_By_Study_Faceted.pdf"), p_var_by_study, width = 11.5, height = 7.0)
ggsave(file.path(fig_dir, "Fig2_Array_Variance_By_Study_Faceted.png"), p_var_by_study, width = 11.5, height = 7.0, dpi = 300)