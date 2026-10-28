# ============================================================
# pipeline/08_discovery_comparison_figure.R
#
# Draws Figure 6, which compares bipolar cells (positive control) with rods on the scDHA
# embedding family. Panel A plots R(k), the mean pairwise ARI across the five scDHA embeddings,
# against k, with a band from the minimum pairwise ARI to its reflection above the mean.
# Panel B plots the seed-1 silhouette against k for the same partitions.
#
# Paper:    Figure 6
# Inputs:   output/tables/07_bc_stability_by_k_extended.csv  (step 07b; if absent,
#           output/tables/07_bc_stability_by_k.csv from step 07)
#           output/tables/07_rods_stability_by_k.csv         (step 07 with CELLTYPE=Rods)
# Outputs:  output/paper/Figure_6.{pdf,png}
# Run:      Rscript run_pipeline.R 08   (or: Rscript pipeline/08_discovery_comparison_figure.R)
# ============================================================

source("config.R"); source("R/utils.R"); source("R/decision_rule.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })

tdir <- TABLE_DIR
# Prefer the extended bipolar sweep from step 07b; fall back to the step 07 table
fbc <- file.path(tdir, "07_bc_stability_by_k_extended.csv")
if (!file.exists(fbc)) fbc <- file.path(tdir, "07_bc_stability_by_k.csv")
frd <- file.path(tdir, "07_rods_stability_by_k.csv")
if (!file.exists(fbc) || !file.exists(frd))
  stop("Need stage-07 outputs for both BC and Rods: run steps 07_bc and 07_rods first.")

# ---- Combine the two curves ----
bc <- read.csv(fbc); bc$cell_type <- "Bipolar cells"; bc$valid <- sccb_valid_k(bc)
rd <- read.csv(frd); rd$cell_type <- "Rods (no established subtypes)"; rd$valid <- sccb_valid_k(rd)
d  <- bind_rows(bc, rd)
d$status <- factor(ifelse(d$valid, "Valid k", "Excluded by Stage-B gates"), levels = c("Valid k", "Excluded by Stage-B gates"))
cols <- c("Bipolar cells"="#1B9E77", "Rods (no established subtypes)"="#D95F02")

# ---- Panel A: R(k) with a band from the minimum pairwise ARI to its reflection above the mean ----
# Filled points: k that pass the Stage B validity gates (as in sccb_verdict); open points: excluded k
pa <- ggplot(d, aes(k, mean_cross_seed_ari, color=cell_type)) +
  geom_hline(yintercept=0.8, linetype="dashed", color="grey40", linewidth=0.3) +
  geom_ribbon(data=subset(d, !is.na(min_cross_seed_ari)),   # minimum pairwise ARI up to R(k)
              aes(ymin=min_cross_seed_ari, ymax=mean_cross_seed_ari, fill=cell_type),
              alpha=0.12, color=NA, inherit.aes=TRUE) +
  geom_line(linewidth=1) + geom_point(aes(shape=status), fill="white", size=2, stroke=0.8) +
  scale_shape_manual(values=c("Valid k"=16, "Excluded by Stage-B gates"=1), name=NULL) +
  scale_color_manual(values=cols, name=NULL) + scale_fill_manual(values=cols, name=NULL, guide="none") +
  scale_x_continuous(breaks=seq(2,18,2)) + ylim(0,1) +
  labs(x="Cluster number k", y="R(k): mean pairwise ARI",
       title="A  Across-embedding reproducibility vs number of clusters") +
  pub_theme() + theme(legend.position="top")

# ---- Panel B: seed-1 silhouette ----
pb <- ggplot(d, aes(k, silhouette, color=cell_type)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  scale_color_manual(values=cols, name=NULL) + scale_x_continuous(breaks=seq(2,18,2)) +
  labs(x="cluster number (k)", y="silhouette width",
       title="B  Cluster quality vs number of clusters") +
  pub_theme() + theme(legend.position="none")

# ---- Assemble and save ----
fig <- pa / pb + plot_layout(heights=c(1,0.9))
save_fig(fig, "Figure_6", PAPER_DIR, width=8, height=8)
message("Saved Figure 6 (discovery comparison).")
