# ============================================================
# pipeline/04d_stability_figures.R
#
# Draws the stability figure for the illustrative rod partition (scDHA + Louvain) from the tables written
# by steps 04b and 04c; it computes nothing new. Panel A shows the cluster sizes at CHOSEN_K for each
# scDHA training seed, panel B the pairwise ARI between those seed partitions, and panel C the
# within-embedding bootstrap stability (mean bootstrap ARI and minimum per-cluster Jaccard) at each k.
#
# Paper:    Supplementary Figure S6
# Inputs:   output/tables/04c_seed_sizes.csv
#           output/tables/04c_seed_ari_matrix.csv
#           output/tables/04b_chosen_partition_stability.csv
# Outputs:  output/paper/Supplementary_Figure_S6.pdf and .png
# Run:      Rscript run_pipeline.R 04d
# ============================================================

source("config.R")
source("R/utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })

tdir <- TABLE_DIR
fig_out <- PAPER_DIR
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

req <- file.path(tdir, c("04c_seed_sizes.csv", "04c_seed_ari_matrix.csv",
                         "04b_chosen_partition_stability.csv"))
if (!all(file.exists(req))) {
  stop("Missing stability tables; run stages 04b and 04c first:\n  ",
       paste(req[!file.exists(req)], collapse = "\n  "))
}

# ---- Panel A: cluster sizes at CHOSEN_K per scDHA training seed, stacked by size rank ----
sizes <- read.csv(file.path(tdir, "04c_seed_sizes.csv"))
sz_long <- sizes %>%
  tidyr::separate_rows(sizes, sep = " / ") %>%
  group_by(seed) %>% mutate(rank = factor(row_number())) %>% ungroup() %>%
  mutate(n = as.numeric(sizes), seed = factor(seed))
pa <- ggplot(sz_long, aes(x = seed, y = n, fill = rank)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = "cluster\n(size rank)") +
  labs(x = "scDHA training seed", y = "cells",
       title = sprintf("A  Seed sensitivity (k = %d)", CHOSEN_K),
       subtitle = "cluster sizes across scDHA training seeds") +
  pub_theme()

# ---- Panel B: pairwise ARI between seed partitions ----
ariM <- read.csv(file.path(tdir, "04c_seed_ari_matrix.csv"), check.names = FALSE)
seeds <- ariM$seed
mat <- as.matrix(ariM[, -1]); rownames(mat) <- seeds; colnames(mat) <- seeds   # CSV header names are R-mangled (X1, ...)
ari_long <- as.data.frame(as.table(mat)); colnames(ari_long) <- c("seed_a", "seed_b", "ARI")
pb <- ggplot(ari_long, aes(seed_a, seed_b, fill = ARI)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI)), size = 3) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                       midpoint = 0.5, limits = c(0, 1), name = "ARI") +
  labs(x = "seed", y = "seed", title = "B  Cross-seed partition agreement",
       subtitle = "pairwise adjusted Rand index") +
  pub_theme() + theme(legend.position = "right")

# ---- Panel C: within-embedding bootstrap stability per k (step 04b) ----
boot <- read.csv(file.path(tdir, "04b_chosen_partition_stability.csv"))
boot_long <- boot %>%
  select(target_k, mean_ari, jacc_min) %>%
  tidyr::pivot_longer(c(mean_ari, jacc_min), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, mean_ari = "mean bootstrap ARI",
                         jacc_min = "min per-cluster Jaccard"))
pc <- ggplot(boot_long, aes(factor(target_k), value, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "grey50") +   # Stage A Jaccard threshold
  scale_fill_manual(values = c("mean bootstrap ARI" = "#7570B3",
                               "min per-cluster Jaccard" = "#1B9E77"), name = NULL) +
  ylim(0, 1) +
  labs(x = "target k", y = "stability",
       title = "C  Within-embedding bootstrap stability",
       subtitle = sprintf("%d bootstrap resamples of 80%% of cells; dashed line, Stage-A threshold 0.6", JACCARD_BOOTSTRAPS)) +
  pub_theme() + theme(legend.position = "top")

# A and B side by side, C below
fig <- (pa | pb) / pc + plot_layout(heights = c(1, 0.9))
save_fig(fig, "Supplementary_Figure_S6", fig_out, width = 11, height = 9)
message("Saved stability supplementary figure to ", fig_out)
