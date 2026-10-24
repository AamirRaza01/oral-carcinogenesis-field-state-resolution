# ==============================================================================
# Script: 11_field_axis.R
# Purpose: Quantify transcriptome-wide standardized effect sizes along the 
#          carcinoma axis (tumor vs. donor normal) and field axis (patient-derived
#          normal vs. donor normal). Evaluate selection bias under normal pooling.
# Reference: Paper Section IV-G, Figure 3
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

suppressPackageStartupMessages({ library(ggplot2) })

REUSE_CACHED  <- TRUE
TOP_CARCINOMA <- 500     # Top carcinoma-effect genes for concordant shift testing
D_MIN         <- 0.5     # Effect size threshold for stratum selection
N_PERM        <- 10000

manifest <- load_manifest()

donor <- manifest$sample[manifest$class == "normal" &
                         manifest$normal_provenance == "donor"]
patient <- manifest$sample[manifest$class == "normal" &
                           manifest$normal_provenance != "donor"]
tumor <- manifest$sample[manifest$class == "cancer"]

# Standardized mean difference with pooled standard deviation
smd <- function(M, a, b) {
  ma <- rowMeans(M[, a, drop = FALSE]); mb <- rowMeans(M[, b, drop = FALSE])
  va <- rowSums((M[, a, drop = FALSE] - ma)^2) / (length(a) - 1)
  vb <- rowSums((M[, b, drop = FALSE] - mb)^2) / (length(b) - 1)
  sp <- sqrt(((length(a) - 1) * va + (length(b) - 1) * vb) /
             (length(a) + length(b) - 2))
  sp[sp <= 0] <- NA
  (ma - mb) / sp
}

sel_rows <- list(); panel_rows <- list()
axes_all <- list(); stat_rows <- list(); point_rows <- list()

cache_ok <- REUSE_CACHED && all(vapply(names(ARMS), function(a)
  file.exists(file.path(res_dir, sprintf("gene_axes_%s.csv", a))), logical(1)))

for (arm in names(ARMS)) {
  message("\nProcessing Arm ", arm, " (", ARM_LABEL[[arm]], ")")
  axf <- file.path(res_dir, sprintf("gene_axes_%s.csv", arm))

  if (cache_ok) {
    ax <- read.csv(axf, stringsAsFactors = FALSE)
    ax <- ax[stats::complete.cases(ax[, c("d_carcinoma", "d_field")]), ]
  } else {
    expr <- load_arm(arm)
    d_carc  <- smd(expr, intersect(tumor, colnames(expr)),   intersect(donor, colnames(expr)))
    d_field <- smd(expr, intersect(patient, colnames(expr)), intersect(donor, colnames(expr)))
    ax <- data.frame(gene = rownames(expr), d_carcinoma = d_carc, d_field = d_field)
    ax <- ax[stats::complete.cases(ax), ]
    write.csv(ax, axf, row.names = FALSE)
  }

  rho <- cor(ax$d_carcinoma, ax$d_field, method = "spearman")
  slope <- unname(coef(lm(d_field ~ d_carcinoma, data = ax))[2])

  strong <- ax[order(-abs(ax$d_carcinoma)), ][seq_len(min(TOP_CARCINOMA, nrow(ax))), ]
  shared <- mean(abs(strong$d_field) >= 0.5 & sign(strong$d_field) == sign(strong$d_carcinoma))

  ax$arm <- ARM_LABEL[[arm]]
  axes_all[[arm]] <- ax
  stat_rows[[arm]] <- data.frame(arm = arm, arm_label = ARM_LABEL[[arm]],
                                 rho = rho, slope = slope, shared = shared)

  # Evaluate stability selection frequency against field effect size
  if (!cache_ok) {
    stratum <- ax[abs(ax$d_carcinoma) >= D_MIN, ]
    for (set in c("DONOR", "FULL")) {
      run <- paste(arm, set, sep = "_")
      f <- file.path(res_dir, sprintf("stability_%s.csv", run))
      if (!file.exists(f)) next
      pi_tab <- read.csv(f, stringsAsFactors = FALSE)
      m <- merge(stratum, pi_tab[, c("gene", "pi")], by = "gene")
      if (nrow(m) < 50) next

      r_field <- cor(m$pi, abs(m$d_field), method = "spearman")
      r_carc  <- cor(m$pi, abs(m$d_carcinoma), method = "spearman")
      sel_rows[[run]] <- data.frame(arm = arm, set = set, run = run,
                                    n_stratum = nrow(m),
                                    rho_pi_field = round(r_field, 4),
                                    rho_pi_carcinoma = round(r_carc, 4),
                                    axis_correlation = round(rho, 4),
                                    shared_fraction = round(shared, 4))
    }
  }

  # Test panel features against carcinoma-matched permutation null
  panels <- readRDS(file.path(res_dir, "stability_panels.rds"))
  for (set in c("DONOR", "FULL")) {
    run <- paste(arm, set, sep = "_")
    if (!run %in% names(panels)) next
    g <- intersect(panels[[run]], ax$gene)
    if (length(g) < 2) next
    point_rows[[run]] <- data.frame(arm_label = ARM_LABEL[[arm]], set = set, gene = g,
                                    d_carcinoma = ax$d_carcinoma[match(g, ax$gene)],
                                    d_field = ax$d_field[match(g, ax$gene)])
    if (cache_ok) next
    obs <- mean(abs(ax$d_field[match(g, ax$gene)]))

    pool <- lapply(g, function(gn) {
      target <- ax$d_carcinoma[ax$gene == gn]
      ax$gene[order(abs(ax$d_carcinoma - target))][seq_len(100)]
    })
    null <- replicate(N_PERM, mean(abs(ax$d_field[match(
      vapply(pool, sample, character(1), size = 1), ax$gene)])))
    z <- (obs - mean(null)) / sd(null)
    p_more <- (1 + sum(null >= obs)) / (N_PERM + 1)
    p_less <- (1 + sum(null <= obs)) / (N_PERM + 1)

    panel_rows[[run]] <- data.frame(arm = arm, set = set, run = run, n_genes = length(g),
                                    mean_abs_field = round(obs, 4),
                                    null_mean = round(mean(null), 4),
                                    z = round(z, 3),
                                    p_more_field = signif(p_more, 4),
                                    p_less_field = signif(p_less, 4))
  }
}

# ------------------------------------------------------------------------------
# Publication Figure 3: Carcinoma vs. Field Axis Scatter and Densities
# ------------------------------------------------------------------------------
INK <- "#0b0b0b"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"

if (length(axes_all)) {
  A <- do.call(rbind, axes_all)
  S <- do.call(rbind, stat_rows)
  P <- if (length(point_rows)) do.call(rbind, point_rows) else NULL
  arm_levels <- unname(ARM_LABEL[names(axes_all)])
  A$arm <- factor(A$arm, levels = arm_levels)
  S$arm_label <- factor(S$arm_label, levels = arm_levels)
  if (!is.null(P)) {
    P$arm_label <- factor(P$arm_label, levels = arm_levels)
    P$set <- factor(P$set, levels = c("DONOR", "FULL"))
  }

  xr <- max(abs(A$d_carcinoma)); yr <- max(abs(A$d_field))
  S$lab <- sprintf("rho %.2f   slope %.2f\n%.0f%% shifted",
                   S$rho, S$slope, 100 * S$shared)

  p <- ggplot(A, aes(d_carcinoma, d_field)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = AQUA, alpha = 0.07) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.5, fill = AQUA, alpha = 0.07) +
    geom_bin2d(bins = 80) +
    scale_fill_gradient(low = "#e8e8e8", high = INK, trans = "log10", guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = INK) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = INK) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.45, colour = INK) +
    geom_abline(data = S, aes(slope = slope, intercept = 0), linewidth = 0.7, colour = BLUE) +
    facet_wrap(~arm_label, nrow = 1) +
    coord_cartesian(xlim = c(-xr, xr), ylim = c(-yr, yr)) +
    labs(x = "carcinoma axis, tumor against donor normal",
         y = "field axis\npatient-derived normal against donor normal") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = INK, linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(colour = INK, size = 10.5),
          legend.position = "bottom", legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(r = 6), lineheight = 1.05),
          plot.margin = margin(6, 8, 6, 10))

  if (!is.null(P)) {
    p <- p + geom_point(data = P, aes(shape = set, colour = set),
                        size = 2.1, stroke = 0.6, fill = "white") +
      scale_shape_manual(values = c(DONOR = 21, FULL = 22),
                         labels = c("donor-set panel gene", "pooled-set panel gene")) +
      scale_colour_manual(values = c(DONOR = BLUE, FULL = ORANGE),
                          labels = c("donor-set panel gene", "pooled-set panel gene"))
    lab <- P[P$set == "DONOR" & abs(P$d_carcinoma) > 1.4, ]
    if (nrow(lab)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = lab, aes(label = gene), size = 2.5, colour = INK,
          segment.colour = "#999999", segment.size = 0.25,
          min.segment.length = 0, max.overlaps = Inf, box.padding = 0.45)
      } else {
        p <- p + geom_text(data = lab, aes(label = gene), size = 2.5, colour = INK,
                           hjust = -0.15, vjust = -0.4)
      }
    }
  }

  p <- p + geom_text(data = S, aes(x = -xr * 0.97, y = yr * 0.97, label = lab),
                     hjust = 0, vjust = 1, size = 2.9, colour = INK, inherit.aes = FALSE)

  ggsave(file.path(fig_dir, "field_axis.pdf"), p, width = 11.5, height = 4.3)
  ggsave(file.path(fig_dir, "field_axis.png"), p, width = 11.5, height = 4.3, dpi = 300)
  write.csv(S[, c("arm", "arm_label", "rho", "slope", "shared")],
            file.path(res_dir, "field_axis_summary.csv"), row.names = FALSE)
}

if (length(sel_rows)) {
  write.csv(do.call(rbind, sel_rows), file.path(res_dir, "field_axis_selection.csv"), row.names = FALSE)
}
if (length(panel_rows)) {
  write.csv(do.call(rbind, panel_rows), file.path(res_dir, "field_axis_panels.csv"), row.names = FALSE)
}