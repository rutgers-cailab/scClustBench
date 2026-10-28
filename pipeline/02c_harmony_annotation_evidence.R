# ============================================================
# pipeline/02c_harmony_annotation_evidence.R
#
# Writes the per-cluster evidence used to annotate the 34 Harmony clusters (resolution 0.6) by hand.
# For every cluster it reports the canonical marker panels (mean expression and percent of cells
# expressing each panel gene), the panel module scores, QC summaries, the top data-driven markers,
# and the labels already in the object (for comparison only). The hand assignments made from these
# tables are stored in data/annotation/harmony_cluster_annotation.csv. The step only reads the object.
#
# Paper:    Methods (reference annotation)
# Inputs:   output/rds/02_annotated.rds
# Outputs:  output/tables/02c_cluster_panel_genes.csv
#           output/tables/02c_cluster_scores_qc.csv
#           output/tables/02c_cluster_by_current_label.csv
#           output/tables/02c_cluster_by_auto_label.csv
#           output/tables/02c_cluster_by_sample.csv
#           output/tables/02c_cluster_top_markers.csv
# Run:      Rscript run_pipeline.R 02c
# ============================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("config.R")
obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"
cl <- obj$seurat_clusters
message("clusters: ", nlevels(cl), "; cells: ", ncol(obj))

# ---- Panel genes: mean normalised expression and percent expressing, per cluster ----
genes <- intersect(unique(unlist(CELLTYPE_PANELS)), rownames(obj))
panel_of <- setNames(rep(names(CELLTYPE_PANELS), lengths(CELLTYPE_PANELS)), unlist(CELLTYPE_PANELS))
X <- GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
dot <- do.call(rbind, lapply(levels(cl), function(k) {
  idx <- which(cl == k); x <- X[, idx, drop = FALSE]
  data.frame(cluster = k, gene = genes, panel = panel_of[genes],
             mean_expr = Matrix::rowMeans(x), pct_expr = 100 * Matrix::rowMeans(x > 0))
}))
write.csv(dot, file.path(TABLE_DIR, "02c_cluster_panel_genes.csv"), row.names = FALSE)

# ---- Panel module scores per cluster (same scoring as step 02), with QC summaries ----
for (p in names(CELLTYPE_PANELS)) {
  g <- intersect(CELLTYPE_PANELS[[p]], rownames(obj))
  obj <- AddModuleScore(obj, features = list(g), name = paste0("ms_", make.names(p)), seed = 1)
}
ms <- grep("^ms_", colnames(obj@meta.data), value = TRUE)
md <- obj@meta.data; md$cluster <- cl
summ <- md %>% group_by(cluster) %>%
  summarise(n = n(), median_umi = median(nCount_RNA), median_genes = median(nFeature_RNA),
            across(all_of(ms), mean), .groups = "drop")
write.csv(summ, file.path(TABLE_DIR, "02c_cluster_scores_qc.csv"), row.names = FALSE)

# ---- Labels and samples per cluster (comparison only) ----
lab <- if ("cell_type" %in% colnames(md)) md$cell_type else NA
write.csv(as.data.frame.matrix(table(cluster = cl, label_current = lab)), file.path(TABLE_DIR, "02c_cluster_by_current_label.csv"))
if ("cell_type_auto" %in% colnames(md))
  write.csv(as.data.frame.matrix(table(cluster = cl, label_auto = md$cell_type_auto)), file.path(TABLE_DIR, "02c_cluster_by_auto_label.csv"))
write.csv(as.data.frame.matrix(table(cluster = cl, sample = md$sample)), file.path(TABLE_DIR, "02c_cluster_by_sample.csv"))

# ---- Top data-driven markers per cluster (at most 400 cells per cluster, for speed) ----
Idents(obj) <- cl
# PrepSCTFindMarkers puts the per-sample SCT models on a common scale before testing
if ("SCT" %in% Assays(obj)) obj <- PrepSCTFindMarkers(obj, verbose = FALSE)
mk <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, max.cells.per.ident = 400,
                     random.seed = 1, verbose = FALSE)
top <- mk %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15, with_ties = FALSE)
write.csv(top, file.path(TABLE_DIR, "02c_cluster_top_markers.csv"), row.names = FALSE)
message("02c done")
