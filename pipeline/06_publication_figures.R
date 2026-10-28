# ============================================================
# pipeline/06_publication_figures.R
#
# Draws paper Figures 2, 3 and 5 and Supplementary Figures S1, S2 and S5 from cached outputs of earlier steps,
# and writes the supplementary score and composition tables. Figure 2 shows QC, the integrated clusters,
# marker genes and the reference cell types. Figure 3 is the Phase 1 (whole retina) Stage A screen: metric
# heatmaps with Stage A exclusions hatched, and the reported partition (Stage B k* from step 07_all, then
# the top-ranked admissible configuration at k*) compared with the reference annotation. Figure 5 is the
# Phase 2 rod Stage A screen and the illustrative rod partition (scDHA + Louvain at CHOSEN_K). Phase 1
# ARI/NMI are taken from the benchmark table as written, not recomputed here.
#
# Paper:    Figure 2, Figure 3, Figure 5, Supplementary Figures S1, S2 and S5 (in the code, objects fig1,
#           fig2, fig3 are saved as Figure_2, Figure_3, Figure_5)
# Inputs:   output/rds/02_annotated.rds, output/rds/02_rods.rds
#           output_full_dataset/tables/poc_scores_matrix_4metric.csv, poc_scores_matrix.csv,
#             poc_all_cluster_assignments.csv, poc_harmony_pca_umap.csv (step 03; jacc_min from 12c_p1)
#           output/tables/07_all_stability_by_k.csv (step 07_all; Phase 1 k*)
#           output/tables/03_scores_matrix_4metric.csv, 03_scores_matrix.csv, 03_all_cluster_assignments.csv,
#             03_harmony_scdha_latent.csv (step 04; jacc_min from 12c_rods)
#           output/tables/04_deg_counts.csv, 04_all_markers.csv (step 05; computed here if absent)
#           output/tables/04_harmony_scdha_umap.csv (rod UMAP from step 05; computed here if absent)
# Outputs:  output/paper/Figure_2, Figure_3, Figure_5 (.png and .pdf)
#           output/paper/Supplementary_Figure_S1, S2, S5 (.png and .pdf)
#           output/paper/source_*.csv (the data behind the figures: Phase 1 and rod scores, rod markers,
#             cell types per sample, rod cluster composition)
# Run:      Rscript run_pipeline.R 06
# ============================================================
source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({
  library(Seurat); library(uwot)
  library(ggplot2); library(patchwork)
  library(dplyr); library(tidyr)
  library(RColorBrewer)
  library(aricode)
  library(cowplot)
})
set.seed(SEED)

# Main and supplementary figures both go to output/paper/
OUT <- PAPER_DIR
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SUPP <- PAPER_DIR
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

T0 <- Sys.time()
# Takes a message; prints it with the elapsed time
log_step <- function(m) message(sprintf("[%s] %s",
                                          format(Sys.time() - T0, digits = 3), m))

# Takes a plot, directory, file stem and size in inches; writes <dir>/<name>.png (600 dpi) and .pdf
save_both <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".png")), p,
          width = w, height = h, units = "in", dpi = 600, bg = "white",
          limitsize = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
          width = w, height = h, units = "in", bg = "white",
          limitsize = FALSE)
  log_step(sprintf("  Saved %s/%s.{png,pdf}", dir, name))
}

# Takes a CSV path (first column cell id); returns a numeric matrix with cell ids as row names
read_lat <- function(p) {
  d <- read.csv(p, check.names = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- d[[1]]; m
}

# ---- Shared aesthetics: method display names and order, colour palettes ----
clean_names <- c("PCA_KMeans"   = "PCA + K-Means",
                 "PCA_GMM"      = "PCA + GMM",
                 "scDHA_KMeans" = "scDHA + K-Means",
                 "scDHA_LL"     = "scDHA + Louvain",
                 "PCA_LL"       = "PCA + Louvain",
                 "scDHA_Leiden" = "scDHA + Leiden",
                 "PCA_Leiden"   = "PCA + Leiden")
method_order <- c("PCA + K-Means", "PCA + GMM",
                   "PCA + Louvain", "PCA + Leiden",
                   "scDHA + K-Means",
                   "scDHA + Louvain", "scDHA + Leiden")
unified_palette <- c("#440154", "#3B528B", "#21908C",
                      "#5DC863", "#FDE725")
# Metric heatmaps (Figures 3A-G and 5A-E): one shared 0-1 scale, ColorBrewer YlGnBu (colorblind- and
# print-safe), so a given colour means the same value in every panel; one colour bar per figure
heat_palette <- RColorBrewer::brewer.pal(9, "YlGnBu")
heat_scale <- function() scale_fill_gradientn(colors = heat_palette, limits = c(0, 1), oob = scales::squish,
                                              na.value = "grey85", breaks = seq(0, 1, 0.25), name = NULL)
# shared_bar(width_cm): the figure's single colour bar, as a grob
shared_bar <- function(width_cm = 6) {
  p <- ggplot(data.frame(x = 1:2, y = 1, v = c(0, 1)), aes(x, y, fill = v)) + geom_tile() + heat_scale() +
    labs(fill = "Metric value (same scale in every panel)") +
    theme_void() + theme(legend.position = "bottom", legend.title.position = "top",
                         legend.title = element_text(size = 11, hjust = 0.5),
                         legend.key.width = unit(width_cm / 5, "cm"), legend.key.height = unit(0.4, "cm"),
                         legend.text = element_text(size = 10))
  bar <- cowplot::get_plot_component(p, "guide-box-bottom", return_all = FALSE)
  cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label("Metric value (one scale for all panels)", size = 11),
                     bar, ncol = 1, rel_heights = c(0.35, 1))
}
ct_pal <- c("Rods" = "#E41A1C", "Cones" = "#377EB8",
             "BC" = "#4DAF4A", "Muller Glia" = "#984EA3",
             "RGC/Neuron" = "#FF7F00", "AC" = "#A65628",
             "RPE" = "#F781BF", "Microglia" = "#999999",
             "Unassigned" = "#CCCCCC")
# display label with the umlaut, matching the relabelled cell types below
names(ct_pal)[names(ct_pal) == "Muller Glia"] <- "M\u00fcller glia"
cl_pal <- setNames(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                       "#A65628","#F781BF","#999999","#66C2A5"),
                    as.character(1:9))

# Takes cluster labels; returns a factor renumbered by cluster size (1 = largest)
size_order_remap <- function(lab) {
  tbl <- sort(table(lab), decreasing = TRUE)
  remap <- setNames(seq_along(tbl), names(tbl))
  factor(remap[as.character(lab)], levels = seq_along(tbl))
}

# The Phase 1 reported partition is chosen by the decision rule in the Figure 3 block below (Stage B k*
# from 07_all_stability_by_k.csv, then the top-ranked Stage A admissible configuration at k*).
source("R/decision_rule.R")

# ============================================================
# Shared inputs: annotated whole-retina object and its UMAP metadata
# ============================================================
log_step("\n========== Loading shared inputs ==========")
obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
# display label: the panels and tables print "Müller glia" (stored label "Muller Glia")
MG <- "M\u00fcller glia"
obj@meta.data$cell_type <- sub("^Muller Glia$", MG, as.character(obj@meta.data$cell_type))
DefaultAssay(obj) <- "SCT"

ump <- Embeddings(obj, "umap")
d_meta <- data.frame(UMAP1 = ump[, 1], UMAP2 = ump[, 2],
                      sample    = obj$sample,
                      time_label = obj$time_label,
                      cluster   = factor(obj$SCT_snn_res.0.6),
                      cell_type = as.character(obj$cell_type))
d_meta <- d_meta[sample(nrow(d_meta)), ]   # shuffle rows so no sample or cluster is always drawn on top

sample_pal <- setNames(
  RColorBrewer::brewer.pal(min(8, length(unique(d_meta$sample))), "Set2")[
    seq_along(unique(d_meta$sample))],
  sort(unique(d_meta$sample)))

n_cl <- length(levels(d_meta$cluster))
cluster_pal <- colorRampPalette(
  RColorBrewer::brewer.pal(11, "Spectral"))(n_cl)
names(cluster_pal) <- levels(d_meta$cluster)

# ============================================================
# Paper Figure 2 (object fig1, saved as Figure_2): QC violins, cluster UMAP, marker dot plot,
# cell-type UMAP and cell-type composition
# ============================================================
log_step("\n========== Figure 2: QC, clusters, markers, cell types ==========")

# Figure 2A: QC violins per sample (UMI, genes, % ribosomal).
# Step 01 filters cells on UMI, gene, % mitochondrial and % hemoglobin thresholds and then removes the
# mitochondrial and hemoglobin genes, so those percentages are about 0 on the saved object and are not shown.
qc_long <- obj@meta.data %>%
  select(sample, nCount_RNA, nFeature_RNA, percent.rb) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("nCount_RNA", "nFeature_RNA", "percent.rb"),
                          labels = c("UMI Counts", "Genes Detected", "% Ribosomal")))

f1a <- ggplot(qc_long, aes(sample, value, fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.size = 0.15, alpha = 0.55, fill = "white") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = sample_pal, guide = "none") +
  pub_theme(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text  = element_text(face = "bold", size = 11),
        plot.title  = element_text(face = "bold", size = 13),
        plot.margin = margin(2, 4, 2, 4, "mm")) +
  labs(title = "A", x = NULL, y = NULL)

# Figure 2B: UMAP coloured by Harmony cluster (SCT_snn_res.0.6)
f1b <- ggplot(d_meta, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.18, alpha = 0.7) +
  scale_color_manual(values = cluster_pal, name = "Cluster", drop = FALSE) +
  pub_theme(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 9),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.30, "cm"),
        plot.margin = margin(2, 4, 2, 4, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1),
                               ncol = 2)) +
  labs(title = "B", x = "UMAP1", y = "UMAP2")

# Figure 2C: canonical marker dot plot by reference cell type
ct_levels_dp <- c("Rods", "Cones", "BC", MG,
                   "RGC/Neuron", "AC", "RPE", "Microglia")
markers_use <- c("Rho", "Nrl", "Pde6h", "Opn1sw",
                  "Vsx2", "Isl1", "Snhg11", "Slc32a1",
                  "Pax6", "Vim", "Apoe", "Rpe65", "Tmem119")
markers_use <- intersect(markers_use, rownames(obj))

obj_dp <- obj
obj_dp@meta.data$cell_type_dp <-
  factor(as.character(obj_dp$cell_type), levels = rev(ct_levels_dp))
Idents(obj_dp) <- "cell_type_dp"
obj_dp <- subset(obj_dp, idents = ct_levels_dp)
DefaultAssay(obj_dp) <- "SCT"

f1c <- DotPlot(obj_dp, features = markers_use, assay = "SCT",
                cols = c("lightgrey", "#3B287C"), dot.scale = 5) +
  RotatedAxis() +
  pub_theme(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
        axis.text.y   = element_text(size = 9),
        axis.title.x  = element_text(face = "bold", size = 10),
        axis.title.y  = element_text(face = "bold", size = 10),
        legend.position = "right",
        legend.title  = element_text(face = "bold", size = 9),
        legend.text   = element_text(size = 8),
        legend.box    = "vertical",
        plot.title    = element_text(face = "bold", size = 13),
        plot.margin   = margin(2, 4, 2, 4, "mm")) +
  labs(title = "C", x = "Marker Genes", y = "Cell Type")

# Figure 2D: UMAP coloured by reference cell type
ct_levels <- c("Rods", "Cones", "BC", MG,
                "RGC/Neuron", "AC", "RPE", "Microglia", "Unassigned")
ct_levels <- ct_levels[ct_levels %in% as.character(obj$cell_type)]   # no empty "Unassigned" legend entry
d_meta$cell_type <- factor(d_meta$cell_type, levels = ct_levels)

f1d <- ggplot(d_meta, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.20, alpha = 0.75) +
  scale_color_manual(values = ct_pal, drop = FALSE, name = "Cell Type") +
  pub_theme(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 10),
        legend.text     = element_text(size = 9),
        plot.margin     = margin(2, 4, 2, 4, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "D", x = "UMAP1", y = "UMAP2")

# Figure 2E: cell-type composition per sample, plus an Overall column after a dashed divider
prop_df <- obj@meta.data %>%
  count(sample, cell_type) %>%
  group_by(sample) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
overall <- obj@meta.data %>%
  count(cell_type) %>%
  mutate(sample = "Overall", pct = 100 * n / sum(n))
prop_full <- bind_rows(prop_df, overall) %>%
  mutate(cell_type = factor(cell_type, levels = ct_levels),
         sample = factor(sample,
                          levels = c(sort(unique(obj$sample)), "Overall")))
n_samples <- length(unique(obj$sample))
dash_x <- n_samples + 0.5

f1e <- ggplot(prop_full, aes(sample, pct, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.78) +
  geom_vline(xintercept = dash_x, linetype = "dashed",
              color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = ct_pal, drop = FALSE, name = "Cell Type") +
  scale_y_continuous(expand = c(0, 0)) +
  pub_theme(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        plot.title   = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 10),
        legend.text  = element_text(size = 9),
        plot.margin  = margin(2, 4, 2, 4, "mm")) +
  labs(title = "E", x = NULL, y = "Proportion (%)")

# Assemble Figure 2: A on top, B | C, then D | E
fig1 <- f1a / (f1b | f1c) / (f1d | f1e) +
  plot_layout(heights = c(0.85, 1.35, 1.10))
save_both(fig1, OUT, "Figure_2", 14, 16)

# ============================================================
# Paper Figure 3 (object fig2, saved as Figure_3): Phase 1 whole-retina Stage A heatmaps and evidence
# ============================================================
log_step("\n========== Figure 3: Phase 1 (heatmap + evidence) ==========")

poc4 <- read.csv("output_full_dataset/tables/poc_scores_matrix_4metric.csv",
                  stringsAsFactors = FALSE)
pocR <- read.csv("output_full_dataset/tables/poc_scores_matrix.csv",
                  stringsAsFactors = FALSE)
poc <- merge(poc4, pocR[, c("method","target_k","silhouette")],
              by = c("method","target_k"))

# Blank (NA) the tiles where clustering did not reach its target k (target_k != actual_k)
poc_miss <- !is.na(poc$actual_k) & poc$actual_k != poc$target_k
metric_cols_poc <- c("rogue", "balance", "jacc_min", "marker_spec",
                      "silhouette", "ari", "nmi")
for (mc in metric_cols_poc) {
  if (mc %in% colnames(poc)) poc[poc_miss, mc] <- NA
}
log_step(sprintf("  Masked %d / %d Phase-1 cells (target_k != actual_k)",
                  sum(poc_miss), nrow(poc)))

# ---- Decision rule: mark Stage A exclusions (hatched) and pick the reported partition ----
# Stage B k* comes from the whole-retina R(k) curve of step 07_all; if that table is missing, k = 8 is used.
# The reported configuration is the top-ranked admissible one at k*, or PCA + Louvain if none is admissible.
sel_poc <- sccb_select(poc)
poc$excluded <- !sel_poc$admissible[match(paste(poc$method, poc$target_k),
                                          paste(sel_poc$method, sel_poc$target_k))]
curve_all_f <- file.path(TABLE_DIR, "07_all_stability_by_k.csv")   # step 07_all (CELLTYPE=ALL)
POC_CHOSEN_K <- if (file.exists(curve_all_f)) sccb_verdict(read.csv(curve_all_f))$k_star else NA
if (is.na(POC_CHOSEN_K)) { POC_CHOSEN_K <- 8
  log_step("  WARNING: no Phase-1 Stage-B curve (07_all); falling back to k = 8") }
at_k <- sccb_select(poc[poc$target_k == POC_CHOSEN_K, ])
POC_CHOSEN_METHOD <- if (at_k$admissible[1]) at_k$method[1] else "PCA_LL"
POC_CHOSEN_LABEL <- sprintf("%s, k* = %d", clean_names[[POC_CHOSEN_METHOD]], POC_CHOSEN_K)
log_step(sprintf("  Stage A excludes %d / %d; reported partition: %s",
                 sum(poc$excluded), nrow(poc), POC_CHOSEN_LABEL))

# ---- Phase 1 partitions and reference labels, aligned by cell id (for panels H to K) ----
parts_poc <- read.csv("output_full_dataset/tables/poc_all_cluster_assignments.csv",
                       check.names = FALSE)
gt <- as.character(obj$cell_type); names(gt) <- colnames(obj)
ids <- parts_poc$cell_id
gt_aligned <- gt[ids]

# ---- Phase 1 ARI/NMI are read from the benchmark table (poc_scores_matrix_4metric.csv) ----
# They are not recomputed here, so the figure always matches the score table reported in the paper.
log_step(sprintf("  Using benchmark-table ARI/NMI (%s k=%d ARI = %.4f)",
                  POC_CHOSEN_METHOD, POC_CHOSEN_K,
                  suppressWarnings(poc$ari[poc$method == POC_CHOSEN_METHOD &
                                           poc$target_k == POC_CHOSEN_K][1])))

# ---- Figure 3A to 3G: metric heatmaps (method x k) ----
# Takes the score table, a metric column and a title; returns a heatmap with Stage A exclusions hatched
# and the reported partition outlined
mk_panel_poc <- function(df, metric_col, title_str) {
  d <- df %>%
    mutate(method_label = factor(clean_names[method],
                                  levels = rev(method_order)),
           k_label = factor(target_k, levels = sort(unique(target_k))),
           is_chosen = method == POC_CHOSEN_METHOD &
                        target_k == POC_CHOSEN_K,
           excluded = if ("excluded" %in% names(df)) excluded else FALSE,
           value = .data[[metric_col]],
           label = ifelse(is.na(value), "—",
                           sprintf("%.2f", value)))
  vmin <- min(d$value, na.rm = TRUE)
  vmax <- max(d$value, na.rm = TRUE)
  d$value_norm <- (d$value - vmin) / max((vmax - vmin), 1e-9)
  d$txt_col <- ifelse(!is.na(d$value_norm) & d$value_norm > 0.55,
                       "grey10", "white")
  # three diagonal hatch strokes per Stage A excluded tile (tile centre +/- 0.5)
  hx <- subset(d, excluded)
  hx <- data.frame(x = as.integer(hx$k_label), y = as.integer(hx$method_label))
  hatch <- do.call(rbind, lapply(list(c(-.5, -.5, .5, .5), c(0, -.5, .5, 0), c(-.5, 0, 0, .5)),
    function(o) if (nrow(hx)) data.frame(x = hx$x + o[1], y = hx$y + o[2],
                                         xend = hx$x + o[3], yend = hx$y + o[4])))
  ggplot(d, aes(k_label, method_label)) +
    geom_tile(aes(fill = value), color = "grey20", linewidth = 0.3) +
    # hatch strokes: a white halo under a black line, so they show on light and dark tiles
    { if (!is.null(hatch)) geom_segment(data = hatch, aes(x = x, y = y, xend = xend, yend = yend),
                                        inherit.aes = FALSE, color = "white", linewidth = 1.0) } +
    { if (!is.null(hatch)) geom_segment(data = hatch, aes(x = x, y = y, xend = xend, yend = yend),
                                        inherit.aes = FALSE, color = "black", linewidth = 0.35) } +
    geom_tile(data = subset(d, is_chosen), fill = NA,
              color = "black", linewidth = 1.4) +
    heat_scale() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    pub_theme(base_size = 11) +
    theme(axis.title.y = element_blank(), axis.line = element_blank(),
          axis.title.x = element_text(size = 10, margin = margin(t = 1)),
          axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none",
          legend.key.height = unit(0.30, "cm"),
          legend.key.width  = unit(1.20, "cm"),
          legend.text = element_text(size = 9),
          panel.grid = element_blank(),
          plot.title = element_text(size = 13, face = "bold",
                                       hjust = 0.5, margin = margin(b = 3)),
          plot.margin = margin(1, 1, 1, 1, "mm")) +
    labs(title = title_str, x = "k")
}

# Initial titles; the final panel letters are set further down
# ROGUE of the unsplit population (step 03), shown in the panel title when available
rb_f <- "output_full_dataset/tables/poc_rogue_baseline.csv"
rogue_title <- if (file.exists(rb_f)) sprintf("ROGUE purity (baseline %.2f)", read.csv(rb_f)$rogue_baseline) else "ROGUE purity"
p1_a <- mk_panel_poc(poc, "rogue",       paste("e —", rogue_title))
p1_b <- mk_panel_poc(poc, "balance",     "f — Cluster-size balance (1 - Gini)")
p1_c <- mk_panel_poc(poc, "jacc_min",    "g — Per-cluster Jaccard min stability")
p1_d <- mk_panel_poc(poc, "marker_spec", "h — Marker specificity")
p1_e <- mk_panel_poc(poc, "silhouette",  "i — Silhouette width")
p1_f <- mk_panel_poc(poc, "ari",         "j — ARI vs reference annotation")
p1_g <- mk_panel_poc(poc, "nmi",         "k — NMI vs reference annotation")

# ---- Figure 3H to 3K: evidence for the reported partition ----
# lab_pll8 holds the reported partition (POC_CHOSEN_METHOD at POC_CHOSEN_K), renumbered by size
ump_pca_f <- "output_full_dataset/tables/poc_harmony_pca_umap.csv"
if (!file.exists(ump_pca_f)) {   # UMAP of the whole-retina Harmony-PCA latent from step 03, computed once and cached
  lat_p1 <- read_lat("output_full_dataset/tables/poc_harmony_pca_latent.csv")
  set.seed(SEED)
  u <- uwot::umap(lat_p1, n_neighbors = 30, min_dist = 0.3, n_threads = N_CORES, verbose = FALSE)
  write.csv(data.frame(cell_id = rownames(lat_p1), UMAP1 = u[, 1], UMAP2 = u[, 2]), ump_pca_f, row.names = FALSE)
}
ump_pca <- read_lat(ump_pca_f)
ump_pca <- ump_pca[ids, , drop = FALSE]
raw_pll8 <- parts_poc[[paste0(POC_CHOSEN_METHOD, "_k", POC_CHOSEN_K)]]; names(raw_pll8) <- parts_poc$cell_id
lab_pll8 <- size_order_remap(raw_pll8[ids])
gt_v <- factor(gt_aligned, levels = ct_levels)

# Takes a UMAP matrix, group labels, palette and titles; returns a UMAP scatter coloured by group
mk_umap <- function(ump, grp, palette, title, subtitle) {
  d <- data.frame(UMAP1 = ump[, 1], UMAP2 = ump[, 2], group = grp)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = group)) +
    geom_point(size = 0.30, alpha = 0.75) +
    scale_color_manual(values = palette, drop = FALSE) +
    pub_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey25"),
          legend.title = element_blank(),
          legend.text  = element_text(size = 8),
          legend.position = "bottom",
          legend.box = "horizontal",
          plot.margin = margin(2, 2, 2, 2, "mm")) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                 nrow = 1)) +
    labs(title = title, subtitle = subtitle)
}

f2a <- mk_umap(ump_pca, lab_pll8, cl_pal, "a",
                paste("Pipeline partition:", POC_CHOSEN_LABEL))
f2b <- mk_umap(ump_pca, gt_v, ct_pal, "b",
                "Reference annotation")

# Figure 3J: confusion matrix, % of each pipeline cluster in each reference cell type
ct_v <- gt_aligned
tab <- table(GroundTruth = ct_v, Cluster = lab_pll8)
pct <- as.data.frame(round(100 * t(t(tab) / colSums(tab)), 1))
pct$Cluster <- factor(pct$Cluster, levels = sort(unique(pct$Cluster)))
pct$GroundTruth <- factor(pct$GroundTruth, levels = ct_levels)

f2c <- ggplot(pct, aes(Cluster, GroundTruth, fill = Freq)) +
  geom_tile(color = "grey50", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", Freq),
                  color = ifelse(Freq > 50, "white", "grey20")),
            size = 2.8, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradient(low = "white", high = "#3B528B",
                       limits = c(0, 100), name = "% of cluster") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_equal() +
  pub_theme(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 9, margin = margin(t = 3)),
        axis.title.y = element_text(face = "bold", size = 9, margin = margin(r = 3)),
        legend.position = "right",
        legend.key.height = unit(1.0, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  labs(title = "C  Confusion matrix",
       x = "Pipeline cluster (1 = largest)", y = "Cell type (reference annotation)")

# Figure 3K: ARI vs the reference annotation for every method at the reported k*
ari_at_k <- poc %>%
  filter(target_k == POC_CHOSEN_K, !is.na(ari)) %>%
  mutate(method_label = factor(clean_names[method], levels = method_order))

f2d <- ggplot(ari_at_k,
               aes(reorder(method_label, ari), ari, fill = method_label)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", ari)),            # value only (axis says ARI), so labels fit
            hjust = -0.1, size = 3.0, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 1.25), expand = c(0, 0),
                      breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  labs(title = sprintf("D  Cell-type recovery at chosen k = %d", POC_CHOSEN_K),
       y = sprintf("ARI vs cell-type annotation (k = %d)", POC_CHOSEN_K))

# Final panel letters: label-free metrics A to E, post hoc ARI/NMI F to G, evidence H to K
p1_a <- p1_a + labs(title = paste("A ", rogue_title))
p1_b <- p1_b + labs(title = "B  Cluster-size balance (1 - Gini)")
p1_c <- p1_c + labs(title = "C  Minimum bootstrap Jaccard")
p1_d <- p1_d + labs(title = "D  Marker specificity")
p1_e <- p1_e + labs(title = "E  Silhouette width")
p1_f <- p1_f + labs(title = "F  ARI vs annotation (post hoc)")
p1_g <- p1_g + labs(title = "G  NMI vs annotation (post hoc)")
f2a  <- f2a  + labs(title = "H", subtitle = paste("Pipeline partition:", POC_CHOSEN_LABEL))
f2b  <- f2b  + labs(title = "I",
                      subtitle = "Reference annotation")
f2c  <- f2c  + labs(title = "J  Confusion matrix",
                      x = "Pipeline cluster (1 = largest)",
                      y = "Cell type (reference annotation)")
f2d  <- f2d  + labs(title = sprintf("K  Cell-type recovery at chosen k = %d",
                                        POC_CHOSEN_K),
                      y = sprintf("ARI vs cell-type annotation (k = %d)",
                                    POC_CHOSEN_K))

# Assemble Figure 3 with cowplot::plot_grid, which sizes each grid cell by rel_widths regardless of axis
# labels and legends, so all heatmaps get the same width (patchwork lets that chrome change column widths).
# Layout: 7 metric heatmaps in two rows of four (A to D, E to G), the four evidence panels (H to K) below.
# tile_key(hatched, horizontal): legend for the heatmap tile marks (hatched = excluded by Stage A,
# gray = achieved k differs from the target or the value could not be computed, boxed = the highlighted partition)
tile_key <- function(hatched = TRUE, horizontal = FALSE, box_label = "Partition reported by the rule (boxed)") {
  lab <- c(if (hatched) "Excluded by Stage A (hatched)",
           "Target k not reached or value not computable (gray)",
           box_label)
  kind <- c(if (hatched) "hatch", "gray", "box")
  n <- length(lab)
  d <- if (horizontal) data.frame(x = c(1, 14, 27)[seq_len(n)], y = 1, lab, kind) else data.frame(x = 1, y = rev(seq_len(n)), lab, kind)
  hx <- subset(d, kind == "hatch")
  g <- ggplot(d, aes(x, y)) +
    geom_tile(data = subset(d, kind != "gray"), fill = "white", colour = "grey20", linewidth = 0.3, width = 0.8, height = 0.6) +
    geom_tile(data = subset(d, kind == "gray"), fill = "grey85", colour = "grey20", linewidth = 0.3, width = 0.8, height = 0.6) +
    geom_tile(data = subset(d, kind == "box"), fill = NA, colour = "black", linewidth = 1.2, width = 0.8, height = 0.6) +
    { if (nrow(hx)) geom_segment(data = data.frame(x = hx$x + c(-.4, -.4, 0), y = hx$y + c(-.3, 0, -.3),
                                                   xend = hx$x + c(.4, 0, .4), yend = hx$y + c(.3, .3, 0)),
                                 aes(x = x, y = y, xend = xend, yend = yend), colour = "black", linewidth = 0.35) } +
    geom_text(aes(x = x + 0.6, label = lab), hjust = 0, size = 3.4) +
    coord_cartesian(xlim = if (horizontal) c(0.5, 27) else c(0.5, 9), ylim = c(0.4, max(d$y) + 0.6), clip = "off") +
    theme_void() + theme(plot.margin = margin(4, 4, 4, 4, "mm"))
  g
}

hm_row1 <- cowplot::plot_grid(p1_a, p1_b, p1_c, p1_d, ncol = 4, align = "hv", axis = "tblr")
hm_row2 <- cowplot::plot_grid(p1_e, p1_f, p1_g,
                              cowplot::plot_grid(tile_key(hatched = TRUE), shared_bar(7), ncol = 1, rel_heights = c(1, 0.45)),
                              ncol = 4, align = "h", axis = "tb")
heatmap_block <- cowplot::plot_grid(hm_row1, hm_row2, ncol = 1, align = "v", axis = "lr")

evidence_block <- cowplot::plot_grid(f2a, f2b, f2c, f2d, ncol = 4, align = "hv", axis = "tblr")

fig2 <- cowplot::plot_grid(
  heatmap_block, evidence_block,
  ncol = 1, rel_heights = c(2, 1.15)
)
fig2 <- fig2 + theme(plot.margin = margin(4, 4, 4, 10, "mm"))   # room for the row labels of the first column
save_both(fig2, OUT, "Figure_3", 22, 16)

# ============================================================
# Paper Figure 5 (object fig3, saved as Figure_5): Phase 2 rod Stage A heatmaps and the illustrative
# rod partition (CHOSEN_METHOD at CHOSEN_K)
# ============================================================
log_step("\n========== Figure 5: rods (heatmap + subclusters) ==========")

rod4 <- read.csv("output/tables/03_scores_matrix_4metric.csv",
                  stringsAsFactors = FALSE)
rodR <- read.csv("output/tables/03_scores_matrix.csv",
                  stringsAsFactors = FALSE)
rod <- merge(rod4, rodR[, c("method","target_k","silhouette")],
              by = c("method","target_k"))
# Blank (NA) the tiles where clustering did not reach its target k
rod_miss <- !is.na(rod$actual_k) & rod$actual_k != rod$target_k
metric_cols_rod <- c("rogue", "balance", "jacc_min", "marker_spec", "silhouette")
for (mc in metric_cols_rod) {
  if (mc %in% colnames(rod)) rod[rod_miss, mc] <- NA
}
log_step(sprintf("  Masked %d / %d Phase-2 cells (target_k != actual_k)",
                  sum(rod_miss), nrow(rod)))

# ---- Figure 5A to 5E: rod metric heatmaps (method x k) ----
# Takes the rod score table, a metric column and a title; returns a heatmap with the illustrative
# partition outlined
mk_panel_rod <- function(df, metric_col, title_str) {
  d <- df %>%
    mutate(method_label = factor(clean_names[method],
                                  levels = rev(method_order)),
           k_label = factor(target_k, levels = sort(unique(target_k))),
           is_chosen = method == CHOSEN_METHOD & target_k == CHOSEN_K,
           value = .data[[metric_col]],
           label = ifelse(is.na(value), "—",
                           sprintf("%.2f", value)))
  vmin <- min(d$value, na.rm = TRUE)
  vmax <- max(d$value, na.rm = TRUE)
  d$value_norm <- (d$value - vmin) / max((vmax - vmin), 1e-9)
  d$txt_col <- ifelse(!is.na(d$value_norm) & d$value_norm > 0.55,
                       "grey10", "white")
  ggplot(d, aes(k_label, method_label)) +
    geom_tile(aes(fill = value), color = "grey20", linewidth = 0.3) +
    geom_tile(data = subset(d, is_chosen), fill = NA,
              color = "black", linewidth = 1.6) +
    heat_scale() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    pub_theme(base_size = 12) +
    theme(axis.title.y = element_blank(), axis.line = element_blank(),
          axis.title.x = element_text(size = 11, margin = margin(t = 1)),
          axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 11),
          legend.position = "none",
          legend.key.height = unit(0.30, "cm"),
          legend.key.width  = unit(1.20, "cm"),
          legend.text = element_text(size = 10),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14, face = "bold",
                                       hjust = 0.5, margin = margin(b = 3)),
          plot.margin = margin(1, 1, 1, 1, "mm")) +
    labs(title = title_str, x = "k")
}

# Initial titles; the final panel letters are set further down
p2_a <- mk_panel_rod(rod, "rogue",       "h — ROGUE purity")
p2_b <- mk_panel_rod(rod, "balance",     "i — Cluster-size balance (1 - Gini)")
p2_c <- mk_panel_rod(rod, "jacc_min",    "j — Per-cluster Jaccard min stability")
p2_d <- mk_panel_rod(rod, "marker_spec", "k — Marker specificity")
p2_e <- mk_panel_rod(rod, "silhouette",  "l — Silhouette width")

# ---- Figure 5F to 5L: illustrative rod partition ----
# Cone-cleaned rods in the step 04 cell order
har_rod <- read_lat("output/tables/03_harmony_scdha_latent.csv")
parts_rod <- read.csv("output/tables/03_all_cluster_assignments.csv",
                       check.names = FALSE)
rods <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
DefaultAssay(rods) <- "RNA"
rods[["RNA"]] <- JoinLayers(rods[["RNA"]])
rod_ids <- parts_rod$cell_id
rods_clean <- subset(rods, cells = rod_ids); rods_clean <- rods_clean[, rod_ids]
rm(rods); gc()

# Rod UMAP: reuse the cached file if present, otherwise compute it (same settings as step 05) and cache it
rod_ump_path <- "output/tables/04_harmony_scdha_umap.csv"
if (file.exists(rod_ump_path)) {
  ump_rod <- read_lat(rod_ump_path)[rod_ids, , drop = FALSE]
} else {
  set.seed(SEED)
  ump_rod <- uwot::umap(har_rod[rod_ids, ],
                          n_neighbors = 30, min_dist = 0.3,
                          n_threads = N_CORES, verbose = FALSE)
  rownames(ump_rod) <- rod_ids; colnames(ump_rod) <- c("UMAP1","UMAP2")
  write.csv(data.frame(cell_id = rownames(ump_rod), ump_rod),
            rod_ump_path, row.names = FALSE)
}

partition_col <- paste0(CHOSEN_METHOD, "_k", CHOSEN_K)
raw_lab <- parts_rod[[partition_col]]; names(raw_lab) <- parts_rod$cell_id
lab_v <- size_order_remap(raw_lab[rod_ids])
time_v <- rods_clean$time_label[rod_ids]
# shared axis limits so the per-timepoint panels are directly comparable with the full UMAP
xlim_all <- range(ump_rod[, 1]) + c(-0.5, 0.5)
ylim_all <- range(ump_rod[, 2]) + c(-0.5, 0.5)

# Takes cell indices, labels and a title; returns the large rod UMAP with a cluster legend
mk_rod_umap_big <- function(idx, lab, title) {
  d <- data.frame(UMAP1 = ump_rod[idx, 1], UMAP2 = ump_rod[idx, 2],
                   cluster = lab)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(size = 0.45, alpha = 0.8) +
    scale_color_manual(values = cl_pal, drop = FALSE, name = "Subcluster") +
    coord_fixed(xlim = xlim_all, ylim = ylim_all) +
    pub_theme(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          legend.position = "right",
          legend.title  = element_text(face = "bold", size = 10),
          legend.text   = element_text(size = 9),
          plot.margin   = margin(3, 3, 3, 3, "mm")) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    labs(title = title, x = "UMAP1", y = "UMAP2")
}

# Takes cell indices, labels and a title; returns a small rod UMAP without legend (one per timepoint)
mk_rod_umap_small <- function(idx, lab, title) {
  d <- data.frame(UMAP1 = ump_rod[idx, 1], UMAP2 = ump_rod[idx, 2],
                   cluster = lab)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(size = 0.30, alpha = 0.75) +
    scale_color_manual(values = cl_pal, drop = FALSE, guide = "none") +
    coord_fixed(xlim = xlim_all, ylim = ylim_all) +
    pub_theme(base_size = 8) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          axis.title    = element_text(size = 7),
          axis.text     = element_text(size = 6),
          plot.margin   = margin(2, 2, 2, 2, "mm")) +
    labs(title = title, x = "UMAP1", y = "UMAP2")
}

cluster_sizes <- as.numeric(table(lab_v))
n_clusters <- length(cluster_sizes)

f3a <- mk_rod_umap_big(seq_along(rod_ids), lab_v,
   sprintf("a — Rod Photoreceptor Subclusters (k = %d)", CHOSEN_K))

# One small UMAP per developmental timepoint (Figure 5G to 5J)
tps <- intersect(c("P15", "P49", "P75", "P91"), unique(time_v))
tp_panels <- lapply(seq_along(tps), function(i) {
  tp <- tps[i]
  idx <- which(time_v == tp)
  mk_rod_umap_small(idx, lab_v[idx], sprintf("%s — %s", letters[i + 1], tp))
})

# Per-cluster marker counts (Figure 5L): read a cached count table, or run FindAllMarkers here
# with the same settings as step 05
deg_path <- file.path(TABLE_DIR, "04_deg_counts.csv")
if (file.exists(deg_path)) {
  deg_counts <- read.csv(deg_path, stringsAsFactors = FALSE)
  deg_counts$cluster <- as.character(deg_counts$cluster)
} else {
  log_step("  Computing FindAllMarkers inline (no cached DEG counts)")
  rcd <- rods_clean
  DefaultAssay(rcd) <- "RNA"
  rcd[["RNA"]] <- JoinLayers(rcd[["RNA"]])
  rcd <- NormalizeData(rcd, verbose = FALSE)
  cl_ids <- as.character(lab_v)
  names(cl_ids) <- rod_ids
  rcd@meta.data$cluster_v <-
    factor(cl_ids[colnames(rcd)],
            levels = as.character(seq_along(cluster_sizes)))
  Idents(rcd) <- "cluster_v"
  markers_inline <- FindAllMarkers(rcd, only.pos = TRUE,
                                    min.pct = 0.10, logfc.threshold = 0.25,
                                    verbose = FALSE)
  markers_sig <- markers_inline %>% filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))
  deg_counts <- markers_sig %>% count(cluster, name = "n_sig_markers") %>%
    mutate(cluster = as.character(cluster)) %>% arrange(cluster)
  write.csv(deg_counts, deg_path,
            row.names = FALSE)
  write.csv(markers_sig, file.path(TABLE_DIR, "04_all_markers.csv"),
            row.names = FALSE)
  rm(rcd); gc()
}
deg_full <- data.frame(cluster = as.character(seq_len(n_clusters)))
deg_full <- merge(deg_full, deg_counts, by = "cluster", all.x = TRUE)
deg_full$n_sig_markers[is.na(deg_full$n_sig_markers)] <- 0
deg_full$cluster <- factor(deg_full$cluster,
                            levels = as.character(seq_len(n_clusters)))

# Figure 5K: cluster sizes
f3_size <- ggplot(
  data.frame(cluster = factor(seq_along(cluster_sizes)),
              n = cluster_sizes,
              pct = round(100 * cluster_sizes / length(rod_ids), 1)),
  aes(cluster, n, fill = cluster)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%s\n(%s%%)",
                                 format(n, big.mark = ","), pct)),
            vjust = -0.25, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cl_pal, guide = "none") +
  pub_theme(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.margin = margin(2, 3, 2, 3, "mm")) +
  expand_limits(y = max(cluster_sizes) * 1.22) +
  labs(title = "F  Subcluster sizes",
       x = "Subcluster", y = "Number of cells")

# Figure 5L: number of significant cluster-defining markers per cluster
f3_deg <- ggplot(deg_full, aes(cluster, n_sig_markers, fill = cluster)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = format(n_sig_markers, big.mark = ",")),
            vjust = -0.25, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cl_pal, guide = "none") +
  pub_theme(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.margin = margin(2, 3, 2, 3, "mm")) +
  expand_limits(y = max(deg_full$n_sig_markers) * 1.22) +
  labs(title = "G  Per-subcluster marker counts",
       x = "Subcluster", y = "Significant cluster-defining DEGs")

# Final panel letters: heatmaps A to E, illustrative partition F to L
p2_a <- p2_a + labs(title = "A  ROGUE purity")
p2_b <- p2_b + labs(title = "B  Balance (1 − Gini)")   # short title so it fits the panel
p2_c <- p2_c + labs(title = "C  Bootstrap Jaccard")   # short title so it fits the panel
p2_d <- p2_d + labs(title = "D  Marker specificity")
p2_e <- p2_e + labs(title = "E  Silhouette width")
f3a  <- f3a  + labs(title = sprintf("F  Rod subclusters (k = %d)",
                                        CHOSEN_K))
tp_panels[[1]] <- tp_panels[[1]] + labs(title = sprintf("G  %s", tps[1]))
tp_panels[[2]] <- tp_panels[[2]] + labs(title = sprintf("H  %s", tps[2]))
tp_panels[[3]] <- tp_panels[[3]] + labs(title = sprintf("I  %s", tps[3]))
tp_panels[[4]] <- tp_panels[[4]] + labs(title = sprintf("J  %s", tps[4]))
f3_size <- f3_size + labs(title = "K  Subcluster sizes")
f3_deg  <- f3_deg  + labs(title = "L  Per-subcluster marker counts")

# Assemble Figure 5 with cowplot::plot_grid (equal heatmap widths, as for Figure 3).
# Layout: the 5 metric heatmaps (A to E) in one row; below, F | G to J (2 x 2) | K over L.
heatmap_block_p2 <- cowplot::plot_grid(p2_a, p2_b, p2_c, p2_d, p2_e,
                                         ncol = 5, align = "hv", axis = "tblr")
heatmap_block_p2 <- cowplot::plot_grid(heatmap_block_p2,
                                         cowplot::plot_grid(tile_key(hatched = FALSE, horizontal = TRUE, box_label = "Illustrative partition (boxed)"),
                                                            shared_bar(8), ncol = 2, rel_widths = c(1.6, 1)),
                                         ncol = 1, rel_heights = c(1, 0.16))

mid_grid <- cowplot::plot_grid(
  tp_panels[[1]], tp_panels[[2]],
  tp_panels[[3]], tp_panels[[4]],
  ncol = 2)
right_col <- cowplot::plot_grid(f3_size, f3_deg, ncol = 1)
subcluster_block <- cowplot::plot_grid(
  f3a, mid_grid, right_col,
  ncol = 3, rel_widths = c(2.4, 2.0, 1.1)
)

fig3 <- cowplot::plot_grid(
  heatmap_block_p2, subcluster_block,
  ncol = 1, rel_heights = c(1, 1.15)
)
fig3 <- fig3 + theme(plot.margin = margin(4, 4, 4, 10, "mm"))   # room for the row labels of the first column
save_both(fig3, OUT, "Figure_5", 20, 13)

# ============================================================
# Supplementary Figures S1, S2 and S5 and supplementary tables
# ============================================================
log_step("\n========== Supplementary Figures S1, S2 and S5 and source tables ==========")

# ------------------------------------------------------------
# Supplementary Figure S1: extended QC and Harmony batch mixing
# ------------------------------------------------------------
log_step("Supplementary_Figure_S1: extended QC + by-sample UMAP")

# S1 A: cells per sample after QC
sa_cells <- obj@meta.data %>% count(sample, name = "n_cells")
sa_cells$sample <- factor(sa_cells$sample, levels = sort(sa_cells$sample))
s1a <- ggplot(sa_cells, aes(sample, n_cells, fill = sample)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n_cells, big.mark = ",")),
             vjust = -0.3, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = sample_pal, guide = "none") +
  expand_limits(y = max(sa_cells$n_cells) * 1.12) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "A  Cells per sample (post-filter)",
       x = NULL, y = "Cells retained after QC")

# S1 B: scDblFinder singlet/doublet counts per sample (placeholder panel if the column is absent)
if ("scDblFinder.class" %in% colnames(obj@meta.data)) {
  dbl_long <- obj@meta.data %>%
    count(sample, scDblFinder.class) %>%
    mutate(scDblFinder.class = factor(scDblFinder.class,
                                        levels = c("singlet", "doublet")))
  s1b <- ggplot(dbl_long, aes(sample, n, fill = scDblFinder.class)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(singlet = "#3B528B", doublet = "#FDE725"),
                       name = "scDblFinder") +
    pub_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "B  Doublet flags per sample",
         x = NULL, y = "Cell count")
} else {
  s1b <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, size = 4,
              label = "scDblFinder.class not in metadata") +
    labs(title = "B  Doublet flags per sample")
}

# S1 C: per-sample QC violins (UMI, genes, % ribosomal). Mitochondrial and hemoglobin genes are removed in
# step 01 after filtering, so percent.mt and percent.hb are about 0 on the saved object and are not shown.
qc_metrics_long <- obj@meta.data %>%
  select(sample, nCount_RNA, nFeature_RNA, percent.rb) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("nCount_RNA", "nFeature_RNA", "percent.rb"),
                          labels = c("UMI counts", "Genes detected",
                                      "% ribosomal")))

s1c <- ggplot(qc_metrics_long, aes(sample, value, fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.size = 0.15, alpha = 0.55, fill = "white") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = sample_pal, guide = "none") +
  pub_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "C  Per-sample QC distributions (UMI, genes, % ribosomal)",
       x = NULL, y = NULL)

# S1 D: Harmony UMAP coloured by sample (batch-mixing check)
s1d <- ggplot(d_meta, aes(UMAP1, UMAP2, color = sample)) +
  geom_point(size = 0.18, alpha = 0.6) +
  scale_color_manual(values = sample_pal, name = "Sample") +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "D  Harmony-corrected UMAP, colored by sample",
       x = "UMAP1", y = "UMAP2")

s_fig_s1 <- (s1a | s1b) / s1c / s1d +
  plot_layout(heights = c(0.9, 0.9, 1.4))
save_both(s_fig_s1, SUPP, "Supplementary_Figure_S1", 14, 14)

# ------------------------------------------------------------
# Supplementary Figure S5: composition of the illustrative rod partition by sample and timepoint
# ------------------------------------------------------------
log_step("Supplementary_Figure_S5: rod subcluster composition")

# S2 A and B: proportion of each cluster of the illustrative rod partition within every sample and every
# timepoint (cluster levels are fixed to 1 to 3, the clusters at CHOSEN_K = 3).
# comp_tab holds the underlying counts and is written to output/paper/source_rod_cluster_composition.csv.
comp_s <- data.frame(sample  = as.character(rods_clean$sample[rod_ids]),
                     cluster = factor(as.character(lab_v), levels = as.character(seq_len(CHOSEN_K))))
s2a <- ggplot(comp_s, aes(sample, fill = cluster)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_fill_manual(values = cl_pal, drop = FALSE, name = "Subcluster") +
  pub_theme(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A - Per-sample subcluster composition", x = NULL, y = "Proportion")

comp_t <- data.frame(time_label = as.character(time_v),
                     cluster = factor(as.character(lab_v), levels = as.character(seq_len(CHOSEN_K))))
s2b <- ggplot(comp_t, aes(time_label, fill = cluster)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_fill_manual(values = cl_pal, drop = FALSE, name = "Subcluster") +
  pub_theme(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13)) +
  labs(title = "B - Per-timepoint subcluster composition", x = NULL, y = "Proportion")

comp_tab <- rbind(
  cbind(level = "sample",    as.data.frame(table(group = comp_s$sample,     subcluster = comp_s$cluster))),
  cbind(level = "timepoint", as.data.frame(table(group = comp_t$time_label, subcluster = comp_t$cluster))))
comp_tab <- comp_tab[comp_tab$Freq > 0, ]
colnames(comp_tab)[colnames(comp_tab) == "Freq"] <- "n_cells"

s_fig_s2 <- (s2a | s2b) + plot_layout(guides = "collect")
save_both(s_fig_s2, SUPP, "Supplementary_Figure_S5", 12, 5)

# ------------------------------------------------------------
# Supplementary Figure S2: Phase 1 ARI across k, resolution plateaus, and sizes of the reported partition
# ------------------------------------------------------------
log_step("Supplementary_Figure_S2: Phase 1 supporting")

# S3 A: Phase 1 ARI vs the reference annotation across k, one line per method
meth_pal <- setNames(RColorBrewer::brewer.pal(7, "Set2"), method_order)   # same colour per method in S3 A and B
ari_long <- poc %>%
  filter(!is.na(ari)) %>%
  mutate(method_label = factor(clean_names[method], levels = method_order))

s3a <- ggplot(ari_long, aes(target_k, ari, color = method_label,
                              group = method_label)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_color_manual(values = meth_pal, name = "Method") +
  scale_x_continuous(breaks = sort(unique(ari_long$target_k))) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right") +
  labs(title = "A  Phase 1 ARI vs reference annotation, by (method, k)",
       x = "Target k", y = "ARI")

# S3 B: target k vs the k actually reached by the graph methods (points off the diagonal missed it)
plateau <- pocR %>%
  filter(method %in% c("PCA_LL", "PCA_Leiden", "scDHA_LL", "scDHA_Leiden")) %>%
  mutate(method_label = clean_names[method],
         hit = target_k == actual_k)

s3b <- ggplot(plateau, aes(target_k, actual_k, color = method_label,
                             shape = hit, group = method_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = meth_pal, name = "Method") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), breaks = c("TRUE", "FALSE"),
                      labels = c("yes", "no"), name = "Target k reached") +
  scale_x_continuous(breaks = sort(unique(plateau$target_k))) +
  scale_y_continuous(breaks = sort(unique(plateau$actual_k))) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right") +
  labs(title = "B  Resolution-plateau diagnostic",
       x = "Target k", y = "Actual k achieved")

# S3 C: cluster sizes of the Phase 1 reported partition
sizes_chosen <- as.numeric(table(lab_pll8))
size_df <- data.frame(cluster = factor(seq_along(sizes_chosen)),
                       n = sizes_chosen,
                       pct = round(100 * sizes_chosen / sum(sizes_chosen), 1))

s3c <- ggplot(size_df, aes(cluster, n, fill = cluster)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%s%%)",
                                  format(n, big.mark = ","), pct)),
             vjust = -0.25, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cl_pal, guide = "none") +
  expand_limits(y = max(sizes_chosen) * 1.18) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "C  Cluster sizes of the reported partition", subtitle = POC_CHOSEN_LABEL,
       x = "Pipeline cluster (1 = largest)", y = "Cells")

s_fig_s3 <- s3a / (s3b | s3c) +
  plot_layout(heights = c(0.9, 1.0))
save_both(s_fig_s3, SUPP, "Supplementary_Figure_S2", 14, 10)

# ============================================================
# Source tables (written to output/paper/)
# ============================================================
log_step("\nSource tables")

# Source table: Phase 1 scores for every (method, k)
poc_export <- poc %>%
  mutate(method_label = clean_names[method]) %>%
  select(method, method_label, target_k, actual_k,
          rogue, balance, jacc_min, marker_spec, silhouette, ari, nmi,
          any_of("sizes_str"))
write.csv(poc_export, file.path(SUPP, "source_phase1_scores.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (Phase 1 scores, %d rows)",
                  SUPP, "source_phase1_scores.csv", nrow(poc_export)))

# Source table: rod scores for every (method, k)
rod_export <- rod %>%
  mutate(method_label = clean_names[method]) %>%
  select(method, method_label, target_k, actual_k,
          rogue, balance, jacc_min, marker_spec, silhouette,
          any_of("sizes_str"))
write.csv(rod_export, file.path(SUPP, "source_rod_scores.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (Phase 2 scores, %d rows)",
                  SUPP, "source_rod_scores.csv", nrow(rod_export)))

# Source table: cluster-defining markers of the illustrative rod partition (copied only if the file exists)
all_markers_path <- file.path(TABLE_DIR, "04_all_markers.csv")
if (file.exists(all_markers_path)) {
  file.copy(all_markers_path, file.path(SUPP, "source_rod_markers.csv"),
             overwrite = TRUE)
  ms <- read.csv(all_markers_path)
  log_step(sprintf("  Wrote %s/%s (%d cluster-defining DEGs)",
                    SUPP, "source_rod_markers.csv", nrow(ms)))
}

# Source table: cell-type composition per sample
ct_per_sample <- obj@meta.data %>%
  count(sample, cell_type) %>%
  group_by(sample) %>%
  mutate(pct = round(100 * n / sum(n), 2)) %>%
  ungroup() %>%
  arrange(sample, desc(n))
write.csv(ct_per_sample, file.path(SUPP, "source_celltype_by_sample.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (cell-type composition per sample)",
                  SUPP, "source_celltype_by_sample.csv"))

# Source table: illustrative rod partition composition per sample and per timepoint (data for Supplementary Figure S5)
write.csv(comp_tab, file.path(SUPP, "source_rod_cluster_composition.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (rod subcluster composition per sample/timepoint)",
                  SUPP, "source_rod_cluster_composition.csv"))

log_step(sprintf("\n=== Step 06 done: %s ===",
                  format(Sys.time() - T0, digits = 3)))
