# ============================================================
# pipeline/01_qc_harmony_integration.R
#
# Builds the integrated mouse retina object from the six public GEO samples. Each sample is
# loaded, cleaned of doublets with scDblFinder, filtered on cell-level QC thresholds, then
# all samples are normalised with SCTransform (one model per sample), reduced by PCA and
# batch-corrected with Harmony on the sample label. Harmony is the only batch correction
# applied to the full dataset; neighbours, clusters and the UMAP are all computed on it.
#
# Paper:    Methods (data processing); the res 0.6 clusters are the units of the reference annotation
# Inputs:   raw 10X directories listed in config.R SAMPLES (under SCCB_RAW_DATA)
# Outputs:  output/rds/01_integrated.rds          (Seurat object with "harmony" reduction)
#           output/tables/01_doublet_summary.csv
#           output/tables/01_qc_filter_summary.csv
#           output/figures/fig1_qc/                (QC violins, diagnostic UMAPs, marker dot plot)
# Run:      Rscript run_pipeline.R 01
# Options:  SCCB_RAW_DATA (folder holding the GEO sample folders; default: parent of the repository)
# ============================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(sctransform)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
})

# SCTransform passes large objects to future workers; raise the export limit to 32 GB
options(future.globals.maxSize = 32000 * 1024^2)
set.seed(SEED)

# ============================================================
# Part 1: Load 10X data and tag metadata
# ============================================================
message("=== Part 1: Loading 10X data ===")

obj_list <- list()
for (i in seq_len(nrow(SAMPLES))) {
  s <- SAMPLES[i, ]
  message("  Loading ", s$sample, " from ", s$dir)

  counts <- read10x_robust(s$dir)
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
  obj$sample     <- s$sample
  obj$animal     <- s$animal
  obj$condition  <- s$condition
  obj$source     <- s$source
  obj$time_label <- s$time_label
  # Prefix barcodes with the sample ID so they stay unique after merging
  obj <- RenameCells(obj, add.cell.id = s$sample)
  Idents(obj) <- "condition"
  obj_list[[s$sample]] <- obj
  message("    Cells: ", ncol(obj), " | Genes: ", nrow(obj))
}

# ============================================================
# Part 2: Doublet removal (per-sample scDblFinder)
# ============================================================
message("\n=== Part 2: Doublet removal ===")

singlet_list   <- list()
doublet_summary <- data.frame()

for (nm in names(obj_list)) {
  message("  Processing ", nm, " ...")
  obj <- obj_list[[nm]]
  n_before <- ncol(obj)

  # dbr.sd = 0 fixes the expected doublet rate at the default for the cell count
  set.seed(SEED)
  sce <- scDblFinder(as.SingleCellExperiment(obj), dbr.sd = 0)
  obj <- as.Seurat(sce, data = NULL)
  obj <- subset(obj, subset = scDblFinder.class == "singlet")
  Idents(obj) <- "condition"

  singlet_list[[nm]] <- obj
  doublet_summary <- rbind(doublet_summary, data.frame(
    sample = nm, before = n_before, after = ncol(obj),
    removed = n_before - ncol(obj),
    pct = round(100 * (n_before - ncol(obj)) / n_before, 2)))
  message("    Kept ", ncol(obj), " / ", n_before, " cells")
}
rm(obj_list); gc()

# ---- Merge all singlets into one object ----
All_Data <- Reduce(function(a, b) merge(a, b, project = "Retina_v2",
                                         merge.data = TRUE), singlet_list)
rm(singlet_list); gc()
write.csv(doublet_summary, file.path(TABLE_DIR, "01_doublet_summary.csv"),
          row.names = FALSE)
message("  Merged: ", ncol(All_Data), " cells x ", nrow(All_Data), " genes")

# ============================================================
# Part 3: Quality control filtering
# ============================================================
message("\n=== Part 3: QC filtering ===")

DefaultAssay(All_Data) <- "RNA"
All_Data[["percent.mt"]] <- PercentageFeatureSet(All_Data, pattern = "^mt-")
All_Data[["percent.hb"]] <- PercentageFeatureSet(All_Data, pattern = "^Hb")
All_Data[["percent.rb"]] <- PercentageFeatureSet(All_Data, pattern = "^Rp[sl]")
gc()

# ---- QC figures BEFORE filtering ----
message("  Pre-QC figures ...")
meta_before <- All_Data@meta.data
p_vln_before <- meta_before %>%
  select(condition, nCount_RNA, nFeature_RNA, percent.mt, percent.hb) %>%
  pivot_longer(-condition, names_to = "metric", values_to = "value") %>%
  ggplot(aes(condition, value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "QC metrics BEFORE filtering", x = "", y = "")
save_fig(p_vln_before, "qc_violin_before", FIG_QC, width = 14, height = 5)

# ---- Apply filters (thresholds in config.R) ----
message("  Applying QC filters ...")
n_before <- ncol(All_Data)
meta <- All_Data@meta.data
keep <- (meta$nCount_RNA   > QC_MIN_COUNTS) &
        (meta$nCount_RNA   < QC_MAX_COUNTS) &
        (meta$nFeature_RNA > QC_MIN_FEATURES) &
        (meta$nFeature_RNA < QC_MAX_FEATURES) &
        (meta$percent.mt   < QC_MAX_PCT_MT) &
        (meta$percent.hb   < QC_MAX_PCT_HB)
All_Data <- subset(All_Data, cells = colnames(All_Data)[keep])
rm(meta, keep); gc()

# Keep genes detected in at least QC_GENE_MIN_CELLS cells
cts <- GetAssayData(All_Data, layer = "counts")
keep_genes <- Matrix::rowSums(cts > 0) >= QC_GENE_MIN_CELLS
All_Data <- subset(All_Data, features = rownames(All_Data)[keep_genes])
rm(cts, keep_genes); gc()

# Drop mitochondrial and haemoglobin genes so they cannot drive clustering
rm_genes <- grepl("^mt-", rownames(All_Data)) | grepl("^Hb", rownames(All_Data))
All_Data <- subset(All_Data, features = rownames(All_Data)[!rm_genes])
rm(rm_genes); gc()

n_after <- ncol(All_Data)
message("  Cells: ", n_before, " -> ", n_after,
        " (removed ", n_before - n_after, ")")

# Recompute QC metrics on the filtered object and log-normalise the RNA assay
All_Data[["percent.mt"]] <- PercentageFeatureSet(All_Data, pattern = "^mt-")
All_Data[["percent.hb"]] <- PercentageFeatureSet(All_Data, pattern = "^Hb")
All_Data <- NormalizeData(All_Data, verbose = FALSE)
gc()

filter_summary <- data.frame(
  metric        = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                     "percent.hb", "gene_min_cells"),
  threshold     = c(paste0(QC_MIN_COUNTS, "-", QC_MAX_COUNTS),
                     paste0(QC_MIN_FEATURES, "-", QC_MAX_FEATURES),
                     paste0("<", QC_MAX_PCT_MT),
                     paste0("<", QC_MAX_PCT_HB),
                     paste0(">=", QC_GENE_MIN_CELLS)),
  cells_before  = n_before,
  cells_after   = n_after)
write.csv(filter_summary, file.path(TABLE_DIR, "01_qc_filter_summary.csv"),
          row.names = FALSE)

# ---- Post-QC violin ----
meta_after <- All_Data@meta.data
p_vln_after <- meta_after %>%
  select(condition, nCount_RNA, nFeature_RNA, percent.mt, percent.hb) %>%
  pivot_longer(-condition, names_to = "metric", values_to = "value") %>%
  ggplot(aes(condition, value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "QC metrics AFTER filtering", x = "", y = "")
save_fig(p_vln_after, "qc_violin_after", FIG_QC, width = 14, height = 5)
rm(meta_before, meta_after); gc()

# ============================================================
# Part 4: SCTransform + PCA + Harmony integration
# ============================================================
message("\n=== Part 4: SCTransform + PCA + Harmony ===")

# Split the RNA assay by sample so SCTransform fits one model per sample
# (the standard Seurat v5 SCTransform plus Harmony workflow).
All_Data[["RNA"]] <- split(All_Data[["RNA"]], f = All_Data$sample)

message("  SCTransform ...")
All_Data <- SCTransform(All_Data, vars.to.regress = "percent.mt",
                         verbose = FALSE)

message("  PCA (npcs = ", N_PCS_INTEGRATION, ") ...")
All_Data <- RunPCA(All_Data, npcs = N_PCS_INTEGRATION, verbose = FALSE)

message("  RunHarmony on PCA (vars_use = 'sample', theta = ",
        HARMONY_THETA, ", lambda = ", HARMONY_LAMBDA,
        ", max_iter = ", HARMONY_MAX_ITER, ") ...")
set.seed(SEED)
All_Data <- harmony::RunHarmony(
  object         = All_Data,
  group.by.vars  = "sample",
  reduction.use  = "pca",
  reduction.save = "harmony",
  dims.use       = 1:N_PCS_INTEGRATION,
  theta          = HARMONY_THETA,
  lambda         = HARMONY_LAMBDA,
  max_iter       = HARMONY_MAX_ITER,
  verbose        = FALSE)

message("  FindNeighbors + FindClusters on harmony reduction ...")
All_Data <- FindNeighbors(All_Data, reduction = "harmony",
                           dims = 1:N_DIMS_NEIGHBORS, verbose = FALSE)
All_Data <- FindClusters(All_Data, resolution = CLUSTER_RESOLUTIONS,
                          verbose = FALSE)

message("  RunUMAP on harmony reduction ...")
All_Data <- RunUMAP(All_Data, reduction = "harmony",
                     dims = 1:N_DIMS_NEIGHBORS, verbose = FALSE)

# Default identity is the clustering at DEFAULT_RESOLUTION (0.6); these clusters are annotated in step 02
Idents(All_Data) <- All_Data[[paste0("SCT_snn_res.", DEFAULT_RESOLUTION)]][, 1]
All_Data$seurat_clusters <- Idents(All_Data)
gc()

# ---- Save integrated object ----
saveRDS(All_Data, file.path(RDS_DIR, "01_integrated.rds"))
message("\n  Saved: output/rds/01_integrated.rds")
message("  Cells: ", ncol(All_Data), " | Genes: ", nrow(All_Data),
        " | Clusters at res ", DEFAULT_RESOLUTION, ": ",
        length(unique(Idents(All_Data))))

# ---- Diagnostic UMAPs ----
message("\n  Diagnostic UMAP figures ...")

p_umap_clusters <- DimPlot(All_Data, reduction = "umap", label = TRUE,
                            label.size = 3, repel = TRUE) +
  pub_theme() +
  labs(title = paste0("Harmony UMAP: Seurat clusters (res ",
                       DEFAULT_RESOLUTION, ")"))
save_fig(p_umap_clusters, "umap_clusters", FIG_QC, width = 8, height = 7)

p_umap_sample <- DimPlot(All_Data, reduction = "umap", group.by = "sample") +
  pub_theme() + labs(title = "Harmony UMAP by sample (batch check)")
save_fig(p_umap_sample, "umap_by_sample", FIG_QC, width = 8, height = 7)

p_umap_time <- DimPlot(All_Data, reduction = "umap", group.by = "time_label") +
  pub_theme() + labs(title = "Harmony UMAP by time point")
save_fig(p_umap_time, "umap_by_timepoint", FIG_QC, width = 8, height = 7)

# ---- Marker dot plot by cluster (supports the cluster annotation used in step 02) ----
DefaultAssay(All_Data) <- "SCT"
markers_present <- intersect(unlist(CELLTYPE_PANELS), rownames(All_Data))
p_dot <- DotPlot(All_Data, features = unique(markers_present),
                  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  pub_theme() +
  labs(title = paste0("Cell-type marker expression by cluster (res ",
                       DEFAULT_RESOLUTION, ")"))
save_fig(p_dot, "dotplot_markers_by_cluster", FIG_QC,
         width = 14, height = 8)

message("\n=== 01_qc_harmony_integration.R complete ===")
