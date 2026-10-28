# ============================================================
# pipeline/08b_bc_positive_control_figure.R
#
# Draws Figure 4, the bipolar cell positive control. It clusters the reference (seed 1)
# scDHA embedding from step 07 at k = K_SHOW and shows: (A) a UMAP of that embedding coloured
# by subcluster; (B) R(k), the mean pairwise ARI across the five scDHA embeddings, against k,
# with the rule's k* (K_RULE) marked; (C) seed-1 silhouette, ROGUE and size balance against k;
# (D) subcluster sizes at k = K_SHOW.
#
# Paper:    Figure 4
# Inputs:   output/tables/12_verdicts.csv (bipolar k*)
#           output/rds/07_bc_positive_control.rds          (step 07 with CELLTYPE=BC)
#           output/tables/07_bc_stability_by_k_extended.csv (step 07b; if absent,
#           output/tables/07_bc_stability_by_k.csv from step 07)
# Outputs:  output/paper/Figure_4.{pdf,png}
# Run:      Rscript run_pipeline.R 08b   (or: Rscript pipeline/08b_bc_positive_control_figure.R)
# ============================================================

source("config.R"); source("R/utils.R"); source("R/decision_rule.R")
suppressPackageStartupMessages({
  library(Seurat); library(uwot); library(ggplot2); library(dplyr)
  library(tidyr); library(patchwork)
})

# The rule's bipolar k* in the scDHA family, from step 12; the figure shows the partition at that k
ver <- read.csv(file.path(TABLE_DIR, "12_verdicts.csv"))
K_RULE <- ver$k_star[ver$case == "retina" & ver$population == "BC" & ver$embedding == "scdha"]  # marked in panels B, C
stopifnot(length(K_RULE) == 1, !is.na(K_RULE))
K_SHOW <- K_RULE  # k of the displayed partition (panels A and D)
REF_SEED <- "1" # reference scDHA embedding for the displayed UMAP

frds <- file.path(RDS_DIR, "07_bc_positive_control.rds")
# Prefer the extended sweep from step 07b; fall back to the step 07 table
fstb <- file.path(TABLE_DIR, "07_bc_stability_by_k_extended.csv")
if (!file.exists(fstb)) fstb <- file.path(TABLE_DIR, "07_bc_stability_by_k.csv")
if (!file.exists(frds) || !file.exists(fstb))
  stop("Need stage-07 BC outputs: run `07` (and optionally `07b`) first.")

pc  <- readRDS(frds)
lat <- pc$seed_lat[[REF_SEED]]         # Harmony-corrected scDHA latent (reference seed)
stb <- read.csv(fstb)

# ---- Louvain partition at the displayed k (same search as step 07) ----
# louvain_to_k(latent, target_k): Louvain labels (named by cell) whose cluster count is closest to
# target_k over a fixed 60-point resolution grid
louvain_to_k <- function(latent, target_k) {
  o <- CreateSeuratObject(matrix(0, 2, nrow(latent), dimnames = list(c("d1", "d2"), rownames(latent))))
  o[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA")
  o <- FindNeighbors(o, reduction = "lat", dims = 1:ncol(latent), verbose = FALSE)
  bl <- NULL; bd <- Inf
  for (r in exp(seq(log(1e-4), log(10), length.out = 60))) {
    set.seed(SEED); o <- FindClusters(o, resolution = r, algorithm = 1, verbose = FALSE)
    lab <- as.integer(Idents(o)); nk <- length(unique(lab))
    if (abs(nk - target_k) < bd) { bd <- abs(nk - target_k); bl <- lab }; if (bd == 0) break
  }
  names(bl) <- rownames(latent); bl
}

lab <- louvain_to_k(lat, K_SHOW)
# Relabel clusters by descending size (cluster 1 = largest) so colours do not depend on Louvain IDs
ord   <- names(sort(table(lab), decreasing = TRUE))
remap <- setNames(seq_along(ord), ord)
sub   <- factor(remap[as.character(lab)], levels = seq_along(ord))

# ---- Panel A: UMAP of the reference embedding, coloured by subcluster ----
set.seed(SEED)
um <- uwot::umap(lat, n_neighbors = 30, min_dist = 0.3, verbose = FALSE)
umd <- data.frame(UMAP1 = um[, 1], UMAP2 = um[, 2], Subcluster = sub)
pal <- setNames(scales::hue_pal()(nlevels(sub)), levels(sub))

pa <- ggplot(umd, aes(UMAP1, UMAP2, color = Subcluster)) +
  geom_point(size = 0.35, alpha = 0.75) +
  scale_color_manual(values = pal, name = "Subcluster",
                     guide = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 2)) +
  labs(title = sprintf("A  Bipolar subclusters (k = %d)", K_SHOW)) +
  pub_theme() +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(size = 8), legend.key.height = unit(3.4, "mm"))

# ---- Panel B: R(k) with a band from the minimum pairwise ARI to its reflection above the mean ----
# Filled points: k that pass the Stage B validity gates (as in sccb_verdict); open points: excluded k
stb$status <- factor(ifelse(sccb_valid_k(stb), "Valid k", "Excluded by Stage-B gates"),
                     levels = c("Valid k", "Excluded by Stage-B gates"))
pb <- ggplot(stb, aes(k, mean_cross_seed_ari)) +
  geom_vline(xintercept = K_RULE, linetype = "dotted", color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_ribbon(aes(ymin = min_cross_seed_ari, ymax = mean_cross_seed_ari),   # minimum pairwise ARI up to R(k)
              fill = "#1B9E77", alpha = 0.15, color = NA) +
  geom_line(color = "#1B9E77", linewidth = 1) +
  geom_point(aes(shape = status), color = "#1B9E77", fill = "white", size = 2, stroke = 0.8) +
  scale_shape_manual(values = c("Valid k" = 16, "Excluded by Stage-B gates" = 1), name = NULL) +
  annotate("text", x = K_RULE - 0.3, y = 0.05, label = sprintf("k* = %d", K_RULE),
           size = 2.6, color = "grey40", hjust = 1) +
  scale_x_continuous(breaks = seq(2, 18, 2)) + ylim(0, 1) +
  labs(x = "Cluster number k", y = "R(k): mean pairwise ARI",
       title = "B  Across-embedding reproducibility") +
  pub_theme() + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# ---- Panel C: seed-1 quality metrics against k ----
ql <- stb %>%
  select(k, silhouette, rogue, balance) %>%
  pivot_longer(-k, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("silhouette", "rogue", "balance"),
                         labels = c("silhouette", "ROGUE purity", "size balance (1 - Gini)")))
qcol <- c("silhouette" = "#3288BD", "ROGUE purity" = "#66C2A5",
          "size balance (1 - Gini)" = "#8073AC")
pc_ <- ggplot(ql, aes(k, value, color = metric)) +
  geom_vline(xintercept = K_RULE, linetype = "dotted", color = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 1) + geom_point(size = 1.8) +
  scale_color_manual(values = qcol, name = NULL) +
  scale_x_continuous(breaks = seq(2, 18, 2)) + ylim(0, 1) +
  labs(x = "cluster number (k)", y = "metric value",
       title = "C  Cluster quality") +
  pub_theme() + theme(plot.title = element_text(face = "bold"),
                      legend.position = "top", legend.text = element_text(size = 8))

# ---- Panel D: subcluster sizes at the displayed partition ----
sz <- as.data.frame(table(sub)); colnames(sz) <- c("Subcluster", "n")
sz$pct <- round(100 * sz$n / sum(sz$n), 1)
pd <- ggplot(sz, aes(Subcluster, n, fill = Subcluster)) +
  geom_col(width = 0.78) +
  geom_text(aes(label = n), vjust = -0.35, size = 2.3) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Subcluster", y = "Number of cells",
       title = sprintf("D  Subcluster sizes (k = %d)", K_SHOW)) +
  pub_theme() + theme(plot.title = element_text(face = "bold"))

# ---- Assemble and save ----
fig <- (pa | (pb / pc_)) / pd + plot_layout(heights = c(1.35, 0.75))
save_fig(fig, "Figure_4",
         PAPER_DIR, width = 11, height = 10)
message(sprintf("Saved BC positive-control figure (k=%d; sizes %s).",
                K_SHOW, paste(sz$n, collapse = "/")))
