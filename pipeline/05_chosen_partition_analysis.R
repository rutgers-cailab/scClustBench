# ============================================================
# pipeline/05_chosen_partition_analysis.R
#
# Describes the illustrative rod partition (CHOSEN_METHOD at CHOSEN_K in config.R: scDHA + Louvain,
# k = 3) taken from the step 04 sweep. Clusters are renumbered by size (1 = largest); the script then
# tabulates cluster sizes and their composition by developmental timepoint, computes a UMAP of the
# Harmony-corrected rod scDHA latent, and finds cluster-defining markers (FindAllMarkers, positive only).
# Per-cluster marker tables are also exported in a format suitable for pathway tools such as IPA.
#
# Paper:    Figure 5 marker counts (via pipeline/06_publication_figures.R); marker audit in step 12b
# Inputs:   output/rds/03_benchmarking.rds  (step 04 bundle: assignments, embeddings, cell ids)
#           output/rds/02_rods.rds          (rod Seurat object, for expression)
# Outputs:  output/rds/04_rod_subclusters.rds
#           output/tables/04_partition_summary.csv
#           output/tables/04_per_timepoint_composition.csv
#           output/tables/04_harmony_scdha_umap.csv
#           output/tables/04_deg_counts.csv
#           output/tables/04_all_markers.csv
#           output/exports/ipa/04_DEGs_cluster<N>.csv
# Run:      Rscript run_pipeline.R 05
# ============================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(uwot)
  library(ggplot2)
  library(patchwork)
  library(dplyr); library(tidyr)
})
set.seed(SEED)

T0 <- Sys.time()
# Takes a message; prints it with the elapsed time
log_step <- function(m) message(sprintf("[%s] %s",
                                          format(Sys.time() - T0, digits = 3), m))

# ============================================================
# Part 1: load the step 04 bundle and the cone-cleaned rod cells
# ============================================================
log_step("=== Part 1: Loading benchmarking bundle + rod object ===")

bundle <- readRDS(file.path(RDS_DIR, "03_benchmarking.rds"))
log_step(sprintf("  Loaded benchmark: %d cells, %d (method, k) rows",
                  bundle$n_cells, nrow(bundle$scores_silrogue)))

rods_full <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
DefaultAssay(rods_full) <- "RNA"
rods_full[["RNA"]] <- JoinLayers(rods_full[["RNA"]])
rods_full <- NormalizeData(rods_full, verbose = FALSE)

# Keep the cone-cleaned rods used in step 04, in the same cell order
clean <- subset(rods_full, cells = bundle$cell_ids)
stopifnot(setequal(colnames(clean), bundle$cell_ids))
clean <- clean[, bundle$cell_ids]
rm(rods_full); gc()

# ============================================================
# Part 2: take the chosen partition and renumber clusters by size
# ============================================================
log_step(sprintf("=== Part 2: Applying %s @ k=%d ===",
                  CHOSEN_METHOD, CHOSEN_K))

partition_col <- paste0(CHOSEN_METHOD, "_k", CHOSEN_K)
stopifnot(partition_col %in% colnames(bundle$cluster_assignments))
raw_lab <- bundle$cluster_assignments[[partition_col]]
names(raw_lab) <- bundle$cluster_assignments$cell_id

# Standardize: cluster 1 = largest, cluster K = smallest
size_order <- order(table(raw_lab), decreasing = TRUE)
ordered_keys <- names(table(raw_lab))[size_order]
remap <- setNames(seq_along(ordered_keys), ordered_keys)
# Build a cell-named metadata frame so AddMetaData matches by barcode (factor()
# drops names, which otherwise triggers "No cell overlap" in Seurat v5).
cluster_md <- data.frame(
  cluster = factor(remap[as.character(raw_lab[colnames(clean)])],
                   levels = seq_along(ordered_keys)),
  row.names = colnames(clean))
clean <- AddMetaData(clean, cluster_md)
log_step("  Cluster sizes after size-ordering:")
print(table(clean$cluster))

sizes <- as.numeric(table(clean$cluster))
write.csv(data.frame(cluster = seq_along(sizes), n = sizes,
                      pct = round(100 * sizes / sum(sizes), 2)),
          file.path(TABLE_DIR, "04_partition_summary.csv"),
          row.names = FALSE)

# ============================================================
# Part 3: cluster composition (% of cells) per timepoint, plus an "All" row
# ============================================================
log_step("=== Part 3: Per-timepoint composition ===")

ptime <- clean@meta.data %>%
  count(time_label, cluster) %>%
  group_by(time_label) %>%
  mutate(pct = 100 * n / sum(n), n_total = sum(n)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(time_label, n_total),
              names_from = cluster,
              names_prefix = "cluster_",
              values_from = pct,
              values_fill = 0) %>%
  mutate(across(starts_with("cluster_"), ~ round(.x, 1)))

all_row <- data.frame(
  time_label = "All",
  n_total    = ncol(clean),
  setNames(as.list(round(100 * sizes / sum(sizes), 1)),
            paste0("cluster_", seq_along(sizes))))
ptime_full <- rbind(ptime, all_row[, colnames(ptime)])
write.csv(ptime_full,
          file.path(TABLE_DIR, "04_per_timepoint_composition.csv"),
          row.names = FALSE)
log_step("  Per-timepoint table:"); print(ptime_full)

# ============================================================
# Part 4: UMAP of the Harmony-corrected rod scDHA latent (for display only)
# ============================================================
log_step("=== Part 4: UMAP on Harmony-scDHA latent ===")

set.seed(SEED)
ump <- uwot::umap(bundle$har_scdha, n_neighbors = 30, min_dist = 0.3,
                   n_threads = N_CORES, verbose = FALSE)
rownames(ump) <- rownames(bundle$har_scdha)
colnames(ump) <- c("UMAP1", "UMAP2")
ump <- ump[colnames(clean), , drop = FALSE]
write.csv(data.frame(cell_id = rownames(ump), ump),
          file.path(TABLE_DIR, "04_harmony_scdha_umap.csv"),
          row.names = FALSE)

# ============================================================
# Part 5: cluster-defining markers (each cluster vs the rest, FindAllMarkers)
# ============================================================
log_step("=== Part 5: Cluster-defining markers (FindAllMarkers) ===")

Idents(clean) <- clean$cluster
markers <- FindAllMarkers(clean, only.pos = TRUE,
                           min.pct = 0.10, logfc.threshold = 0.25,
                           verbose = FALSE)
markers_sig <- markers %>% filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

deg_counts <- markers_sig %>%
  count(cluster, name = "n_sig_markers") %>%
  arrange(cluster)
write.csv(deg_counts, file.path(TABLE_DIR, "04_deg_counts.csv"),
          row.names = FALSE)
log_step("  DEG counts per cluster:"); print(deg_counts)

write.csv(markers_sig, file.path(TABLE_DIR, "04_all_markers.csv"),
          row.names = FALSE)

# One marker table per cluster, for pathway analysis tools such as IPA
ipa_dir <- file.path(EXPORT_DIR, "ipa")
dir.create(ipa_dir, recursive = TRUE, showWarnings = FALSE)
for (cl in unique(markers_sig$cluster)) {
  sub <- markers_sig %>% filter(cluster == cl)
  write.csv(sub, file.path(ipa_dir,
                            sprintf("04_DEGs_cluster%s.csv", cl)),
            row.names = FALSE)
}
log_step(sprintf("  IPA inputs saved to %s", ipa_dir))

# ============================================================
# Part 6: save the rod object with cluster labels and UMAP
# ============================================================
log_step("=== Part 6: Saving rod-subcluster bundle ===")

clean[["umap_har_scdha"]] <- CreateDimReducObject(
  embeddings = ump, key = "UMAPHARSC_", assay = "RNA")

saveRDS(list(clean      = clean,
              cluster    = clean$cluster,
              umap       = ump,
              ptime      = ptime_full,
              sizes      = sizes,
              deg_counts = deg_counts,
              markers    = markers_sig),
         file.path(RDS_DIR, "04_rod_subclusters.rds"))
log_step("  Saved: output/rds/04_rod_subclusters.rds")

log_step("=== Step 05 done ===")
log_step(sprintf("Total wall-time: %s",
                  format(Sys.time() - T0, digits = 3)))
