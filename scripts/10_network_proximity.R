# ==============================================================================
# Script: 10_network_proximity.R
# Purpose: Quantify shortest-path network proximity between selected panel genes
#          and known HNSCC driver genes within the human physical interactome
#          against a degree-matched permutation null model.
# ==============================================================================

source("00_oral_common.R")
set.seed(42)

suppressPackageStartupMessages({
  library(igraph); library(data.table)
})

options(timeout = 3600)

STRING_SCORE  <- 700     # High-confidence physical interaction cutoff
N_NULL        <- 1000    # Degree-matched permutation iterations
N_DEGREE_BINS <- 20

net_dir <- file.path(base_dir, "network"); dir.create(net_dir, showWarnings = FALSE)

# Curated significantly mutated genes in HNSCC (TCGA Network, Nature 2015)
DRIVERS <- c("TP53", "FAT1", "CDKN2A", "PIK3CA", "NOTCH1", "KMT2D", "NSD1",
             "CASP8", "AJUBA", "HLA-A", "TGFBR2", "HRAS", "FBXW7", "RB1",
             "PIK3R1", "TRAF3", "NFE2L2", "CUL3", "PTEN", "EP300", "ZNF750",
             "RAC1")

# ------------------------------------------------------------------------------
# 1. Construct physical protein-protein interaction network
# ------------------------------------------------------------------------------
links_url <- "https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz"
info_url  <- "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
links_gz <- file.path(net_dir, basename(links_url))
info_gz  <- file.path(net_dir, basename(info_url))

if (!file.exists(links_gz)) download.file(links_url, links_gz, mode = "wb")
if (!file.exists(info_gz))  download.file(info_url,  info_gz,  mode = "wb")

links <- fread(links_gz)
info  <- fread(info_gz)
setnames(info, old = names(info)[1:2], new = c("protein_id", "symbol"))

links <- links[combined_score >= STRING_SCORE]
id2sym <- setNames(info$symbol, info$protein_id)
links[, `:=`(a = id2sym[protein1], b = id2sym[protein2])]
links <- links[!is.na(a) & !is.na(b) & a != b]

g <- graph_from_data_frame(unique(links[, .(a, b)]), directed = FALSE)
g <- simplify(g)
comp <- components(g)
g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))

deg <- degree(g)
all_genes <- names(deg)
drivers_in <- intersect(DRIVERS, all_genes)
if (length(drivers_in) < 5) stop("Insufficient driver genes mapped to network.")

# ------------------------------------------------------------------------------
# 2. Closest-distance network proximity with degree-matched null
# ------------------------------------------------------------------------------
D <- distances(g, v = drivers_in)
min_to_driver <- apply(D, 2, function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else min(x) })

# Stratify genes into degree bins for matched permutation
brk <- unique(quantile(deg, probs = seq(0, 1, length.out = N_DEGREE_BINS + 1)))
deg_bin <- cut(deg, breaks = brk, include.lowest = TRUE, labels = FALSE)
bin_members <- split(all_genes, deg_bin)

proximity <- function(genes) mean(min_to_driver[genes], na.rm = TRUE)

null_draw <- function(genes) {
  bins <- deg_bin[match(genes, all_genes)]
  vapply(bins, function(b) sample(bin_members[[as.character(b)]], 1L), character(1))
}

panels <- readRDS(file.path(res_dir, "stability_panels.rds"))
rows <- list(); gene_rows <- list()

for (run in names(panels)) {
  p <- panels[[run]]
  mapped <- intersect(p, all_genes)
  if (!length(mapped)) next

  obs <- proximity(mapped)
  null <- replicate(N_NULL, proximity(null_draw(mapped)))
  z <- (obs - mean(null)) / sd(null)
  pv <- (1 + sum(null <= obs)) / (N_NULL + 1)

  rows[[run]] <- data.frame(
    run = run, arm = substr(run, 1, 1), set = sub("^._", "", run),
    n_panel = length(p), n_mapped = length(mapped),
    proximity = round(obs, 4), null_mean = round(mean(null), 4),
    null_sd = round(sd(null), 4), z = round(z, 3), p_value = signif(pv, 4))

  for (gn in mapped) {
    d <- D[, gn]
    d <- d[is.finite(d)]
    nearest <- if (length(d)) names(which.min(d)) else NA_character_
    gene_rows[[paste(run, gn)]] <- data.frame(
      run = run, gene = gn, degree = unname(deg[gn]),
      distance_to_nearest_driver = unname(min_to_driver[gn]),
      nearest_driver = nearest)
  }
}

prox <- do.call(rbind, rows)
write.csv(prox, file.path(res_dir, "network_proximity.csv"), row.names = FALSE)
write.csv(do.call(rbind, gene_rows), file.path(res_dir, "panel_gene_to_driver.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Concordance comparison between cohort definitions
# ------------------------------------------------------------------------------
d_prox <- prox$proximity[prox$set == "DONOR"]
f_prox <- prox$proximity[prox$set == "FULL"]
if (length(d_prox) >= 2 && length(f_prox) >= 2) {
  suppressWarnings(wilcox.test(d_prox, f_prox))
}