# ============================================================
# pipeline/02_cell_type_annotation.R
#
# Assigns a cell class to every cell of the integrated object. The reference label comes from
# the hand annotation of the 34 Harmony clusters at resolution 0.6
# (data/annotation/harmony_cluster_annotation.csv). An automated assignment is also computed
# for comparison: each cluster gets the CELLTYPE_PANELS marker panel (config.R) with the
# highest cluster-mean module score, or "Unassigned" if that score is below CELLTYPE_MIN_SCORE;
# it is stored as cell_type_auto. The step then saves the rod subset used in Phase 2.
#
# Paper:    Methods (reference annotation); labels used as ground truth in Phase 1
# Inputs:   output/rds/01_integrated.rds
#           data/annotation/harmony_cluster_annotation.csv
# Outputs:  output/rds/02_annotated.rds            (cell_type = reference, cell_type_auto = module score)
#           output/rds/02_rods.rds
#           output/tables/02_celltype_assignments.csv  (automated per-cluster scores)
#           output/tables/02_reference_annotation.csv  (per-cell reference and automated labels)
#           output/tables/02_celltype_summary.csv
#           output/figures/fig2_celltypes/
# Run:      Rscript run_pipeline.R 02
# ============================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(ggplot2); library(tibble); library(stringr); library(purrr)
  library(patchwork)
})
set.seed(SEED)

# ============================================================
# Part 1: Load and module-score every cell-type panel
# ============================================================
message("=== Part 1: Loading object + scoring marker panels ===")

obj <- readRDS(file.path(RDS_DIR, "01_integrated.rds"))
DefaultAssay(obj) <- "SCT"
message("  Cells: ", ncol(obj), " | clusters: ",
        length(levels(Idents(obj))))

panel_score_cols <- character(0)
for (ct in names(CELLTYPE_PANELS)) {
  genes_present <- intersect(CELLTYPE_PANELS[[ct]], rownames(obj))
  # A module score from fewer than two genes is not meaningful; skip the panel
  if (length(genes_present) < 2) {
    message("  Skipping ", ct, ": only ", length(genes_present),
            " markers present")
    next
  }
  score_col <- paste0("score_", gsub("[^A-Za-z0-9]", "_", ct))
  obj <- AddModuleScore(obj, features = list(genes_present),
                         name = paste0(score_col, "_"),
                         seed = SEED)
  # AddModuleScore appends "1" to the name; rename to score_<panel>
  obj@meta.data[[score_col]] <- obj@meta.data[[paste0(score_col, "_1")]]
  obj@meta.data[[paste0(score_col, "_1")]] <- NULL
  panel_score_cols <- c(panel_score_cols, score_col)
  message("  Scored ", ct, " (", length(genes_present), " markers)")
}

# ============================================================
# Part 2: Automated cluster-mean assignment
# ============================================================
message("\n=== Part 2: Cluster-mean assignment ===")

cluster_means <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(panel_score_cols), mean, .names = "{.col}"),
             n = n(), .groups = "drop")

# Map score column names back to panel names; each cluster takes the panel with the highest mean score
panel_to_ct <- setNames(names(CELLTYPE_PANELS),
                         paste0("score_",
                                 gsub("[^A-Za-z0-9]", "_",
                                       names(CELLTYPE_PANELS))))

cluster_assignment <- cluster_means %>%
  rowwise() %>%
  mutate(top_score = max(c_across(all_of(panel_score_cols)), na.rm = TRUE),
         top_panel = panel_score_cols[which.max(c_across(
                       all_of(panel_score_cols)))],
         cell_type = ifelse(top_score >= CELLTYPE_MIN_SCORE,
                             panel_to_ct[top_panel], "Unassigned")) %>%
  ungroup() %>%
  select(seurat_clusters, n, top_panel, top_score, cell_type,
          everything())

write.csv(cluster_assignment,
          file.path(TABLE_DIR, "02_celltype_assignments.csv"),
          row.names = FALSE)
message("  Cluster assignments:")
print(cluster_assignment %>% select(seurat_clusters, n, cell_type, top_score),
       n = Inf)

# Apply the automated cluster labels to cells. AddMetaData with a named vector is used
# because the `obj$col <- vec` setter can fail with "No cell overlap" on unnamed vectors.
ct_map <- setNames(cluster_assignment$cell_type,
                    as.character(cluster_assignment$seurat_clusters))
ct_vec <- ct_map[as.character(obj$seurat_clusters)]
names(ct_vec) <- colnames(obj)
obj <- AddMetaData(obj, metadata = ct_vec, col.name = "cell_type")
obj$cell_type <- factor(obj$cell_type,
                         levels = c(names(CELLTYPE_PANELS), "Unassigned"))

# ============================================================
# Part 3: Reference annotation, save annotated object and rod subset
# ============================================================
message("\n=== Part 3: Saving annotated object + rod subset ===")

# ---- Reference annotation ----
# Each Harmony cluster (resolution 0.6, the 34 clusters above) is assigned by hand to one of the
# eight classes from the canonical marker panels (per-cluster evidence:
# pipeline/02c_harmony_annotation_evidence.R; assignments and deciding markers:
# data/annotation/harmony_cluster_annotation.csv). The automated assignment is kept as cell_type_auto.
ann_path <- "data/annotation/harmony_cluster_annotation.csv"
# The CSV is keyed by res 0.6 cluster ID, so seurat_clusters must be exactly that clustering
stopifnot(file.exists(ann_path), identical(as.character(obj$seurat_clusters), as.character(obj[[paste0("SCT_snn_res.", DEFAULT_RESOLUTION)]][, 1])))
ann <- read.csv(ann_path, stringsAsFactors = FALSE)
stopifnot(setequal(as.character(ann$cluster), levels(obj$seurat_clusters)), all(ann$cell_type %in% names(CELLTYPE_PANELS)))
obj@meta.data$cell_type_auto <- as.character(obj$cell_type)
obj@meta.data$cell_type <- setNames(ann$cell_type, ann$cluster)[as.character(obj$seurat_clusters)]
write.csv(data.frame(cell_id = colnames(obj), harmony_cluster = as.character(obj$seurat_clusters), cell_type = obj$cell_type,
                     cell_type_auto = obj$cell_type_auto), file.path(TABLE_DIR, "02_reference_annotation.csv"), row.names = FALSE)
message(sprintf("  Reference annotation from %s: %d cells; automated rule agrees for %.1f%%", ann_path, ncol(obj),
                100 * mean(obj$cell_type == obj$cell_type_auto)))

# ---- Cell-type composition table (reference labels) ----
ct_counts <- table(obj$cell_type)
ct_summary <- data.frame(cell_type = names(ct_counts),
                          n         = as.numeric(ct_counts),
                          pct       = round(100 * as.numeric(ct_counts) /
                                              sum(ct_counts), 2))
write.csv(ct_summary, file.path(TABLE_DIR, "02_celltype_summary.csv"),
          row.names = FALSE)
message("\n  Cell-type composition (reference annotation):")
print(ct_summary)

saveRDS(obj, file.path(RDS_DIR, "02_annotated.rds"))
message("\n  Saved: output/rds/02_annotated.rds")

# Rod subset defined by the reference annotation (input to Phase 2)
rods <- subset(obj, cells = colnames(obj)[obj$cell_type == "Rods"])
saveRDS(rods, file.path(RDS_DIR, "02_rods.rds"))
message("  Saved: output/rds/02_rods.rds  (", ncol(rods), " rod cells)")

# ============================================================
# Part 4: Cell-type figures
# ============================================================
message("\n=== Part 4: Figures ===")

# UMAP colored by cell type
p_ct <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
                 cols = CELLTYPE_COLORS, label = TRUE, label.size = 3,
                 repel = TRUE) +
  pub_theme() +
  labs(title = paste0("Harmony UMAP: annotated cell types (n = ",
                       format(ncol(obj), big.mark = ","), ")"))
save_fig(p_ct, "umap_celltypes", FIG_CT, width = 8.5, height = 7)

# Stacked bar of cell-type proportions per sample
prop_df <- obj@meta.data %>%
  count(sample, cell_type) %>%
  group_by(sample) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
all_row <- obj@meta.data %>%
  count(cell_type) %>%
  mutate(sample = "Overall", pct = 100 * n / sum(n))
prop_full <- bind_rows(prop_df, all_row)

p_bar <- ggplot(prop_full, aes(x = sample, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = CELLTYPE_COLORS, drop = FALSE) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Cell-type proportions by sample",
       x = "Sample", y = "Percent of cells")
save_fig(p_bar, "celltype_proportions", FIG_CT, width = 9, height = 6)

# Marker dot plot organized by cell type
DefaultAssay(obj) <- "SCT"
ordered_markers <- unique(unlist(CELLTYPE_PANELS))
ordered_markers <- intersect(ordered_markers, rownames(obj))
p_dot <- DotPlot(obj, features = ordered_markers, group.by = "cell_type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  pub_theme() +
  labs(title = "Canonical markers by annotated cell type",
       x = "Gene", y = "Cell type")
save_fig(p_dot, "dotplot_markers_by_celltype", FIG_CT, width = 16, height = 6)

message("\n=== 02_cell_type_annotation.R complete ===")
