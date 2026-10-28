# ============================================================
# pipeline/10_synthetic_figure.R
#
# Draws Figure 7, the validation on synthetic (Splatter) data. Panel A shows UMAPs of
# representative simulated datasets coloured by their true group (null: one population;
# structured: three groups; rare: one small group). Panel B plots the ARI of the seed-1
# partition against the true groups for each k, and panel C plots R(k), the mean pairwise ARI
# across the five scDHA embeddings. Panel A is left out if the step 11 table is missing.
#
# Paper:    Figure 7
# Inputs:   output/tables/11_synthetic_umap_coords.csv (step 11; optional, for panel A)
#           output/tables/09_synthetic_validation.csv  (step 09)
# Outputs:  output/paper/Figure_7.{pdf,png}
# Run:      Rscript run_pipeline.R 10   (or: Rscript pipeline/10_synthetic_figure.R)
# ============================================================

source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

fu <- file.path(TABLE_DIR, "11_synthetic_umap_coords.csv")
fv <- file.path(TABLE_DIR, "09_synthetic_validation.csv")
if (!file.exists(fv)) stop("Run stage 09 first.")

# ---- Panel A: UMAPs of the simulated data, one facet per scenario in step 11 order ----
gcols <- c("group 1"="#B0B7BE","group 2"="#3288BD","group 3"="#E7298A")
if (file.exists(fu)) {
  um <- read.csv(fu)
  um$scenario <- factor(um$scenario, levels=unique(um$scenario[order(um$order)]))
  pa <- ggplot(um, aes(UMAP1, UMAP2, color=group)) +
    geom_point(size=0.25, alpha=0.6) +
    facet_wrap(~scenario, nrow=1, scales="free") +
    scale_color_manual(values=gcols, name="true group", guide=guide_legend(override.aes=list(size=2.5,alpha=1))) +
    labs(title="A  Synthetic datasets (colored by ground-truth group)") +
    pub_theme() + theme(axis.text=element_blank(), axis.ticks=element_blank(),
      axis.title=element_blank(), legend.position="right",
      panel.border=element_rect(color="grey85", fill=NA), strip.text=element_text(size=9))
} else { pa <- NULL }

# ---- Panels B and C: recovery of the truth and R(k) against k ----
# Display labels per scenario; the two null replicates share one label and colour
d <- read.csv(fv)
lab <- c(null_A="Nulls (2 datasets)", null_B="Nulls (2 datasets)",
         weak_3="Weak (de = 0.03)", medium_3="Medium (de = 0.10)",
         strong_3="Strong (de = 0.30)", rare_3="Rare 5% group (de = 0.15)")
d$label <- factor(lab[d$scenario],
  levels=c("Nulls (2 datasets)","Weak (de = 0.03)","Rare 5% group (de = 0.15)","Medium (de = 0.10)","Strong (de = 0.30)"))
cols <- c("Nulls (2 datasets)"="grey60","Weak (de = 0.03)"="#66C2A5","Rare 5% group (de = 0.15)"="#E7298A",
          "Medium (de = 0.10)"="#3288BD","Strong (de = 0.30)"="#1B9E77")
# Only k values that the reference (first-seed) partition actually reached are plotted
d <- d[d$actual_k == d$k, ]

pb <- ggplot(d, aes(k, ari_vs_truth, color=label, group=scenario)) +
  geom_vline(xintercept=3, linetype="dashed", color="grey70") +
  geom_line(linewidth=0.9, alpha=0.9) + geom_point(size=1.8) +
  scale_color_manual(values=cols, name=NULL, guide=guide_legend(nrow=1)) + scale_x_continuous(breaks=2:6) + ylim(0,1) +
  labs(x="Cluster number k", y="ARI vs. ground truth",
       title="B  Recovery of true structure",
       subtitle="dashed line: true k = 3") +
  pub_theme() + theme(legend.text=element_text(size=8))

pc <- ggplot(d, aes(k, mean_cross_seed_ari, color=label, group=scenario)) +
  geom_vline(xintercept=3, linetype="dashed", color="grey70") +
  geom_hline(yintercept=0.8, linetype="dotted", color="grey30") +
  geom_line(linewidth=0.9, alpha=0.9) + geom_point(size=1.8) +
  scale_color_manual(values=cols, name=NULL, guide=guide_legend(nrow=1)) + scale_x_continuous(breaks=2:6) + ylim(0,1) +
  labs(x="Cluster number k", y="R(k): mean pairwise ARI across embeddings",
       title="C  Across-embedding reproducibility", subtitle="dotted line: \u03c4 = 0.8") +
  pub_theme()

# ---- Assemble and save: panels B and C share one legend across the full row width ----
bc <- (pb | pc) + plot_layout(guides="collect") &
      theme(legend.position="bottom", legend.text=element_text(size=8))
fig <- if (!is.null(pa)) pa / bc + plot_layout(heights=c(0.8,1)) else bc
save_fig(fig, "Figure_7", PAPER_DIR, width=11, height=9)
message("Saved Figure 7 (synthetic validation: data UMAPs + recovery).")
