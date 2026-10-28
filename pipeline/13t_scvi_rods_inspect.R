# ============================================================
# pipeline/13t_scvi_rods_inspect.R
#
# Characterises the scVI rod k = 3 split (R = 0.84, but one of the five scVI
# embeddings does not reach k = 3, so the split fails the Stage B validity
# gates; the scDHA and PCA families reject the rods as well). Each of the five scVI rod
# embeddings from 13s is clustered at k = 3 (Louvain), and the seed-1 partition
# is profiled: cluster sizes, sample composition, sex-gene detection (Xist/Tsix
# and Y-linked genes), ribosomal fraction, median UMI, top markers, and overlap
# with the scDHA Louvain k = 3 rod partition from step 04. The assignments
# table is the input of 17b and 17c.
#
# Paper:    scVI rod split (Stage B rods, embedding-dependent result)
# Inputs:   output/rds/07_rods_positive_control.rds
#           output/rds/13s_scvi_rods_latents.rds (from 13s; if missing, the five models are retrained)
#           output/tables/03_all_cluster_assignments.csv (rod partitions written by 04_phase2_rod_benchmark.R)
# Outputs:  output/tables/13t_scvi_rods_assignments.csv (cell x seed, k = 3)
#           output/tables/13t_scvi_rods_profile.csv
#           output/tables/13t_scvi_rods_markers.csv
#           output/tables/13t_scvi_rods_overlap.csv
# Run:      Rscript run_pipeline.R 13t   or   Rscript pipeline/13t_scvi_rods_inspect.R
# Options:  RETICULATE_PYTHON  python with scvi-tools (needed only if the 13s latents file is missing)
# ============================================================
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages({ library(aricode); library(dplyr) })
pc  <- readRDS(file.path(RDS_DIR, "07_rods_positive_control.rds"))
obj <- sccb_prep(pc$bc, "LogNormalize", 3000)

# ---- k = 3 partition of each scVI embedding ----
lat_f <- file.path(RDS_DIR, "13s_scvi_rods_latents.rds")    # the five scVI embeddings saved by 13s
lats <- if (file.exists(lat_f)) readRDS(lat_f) else NULL
lab <- list()
for (s in SCDHA_SEED_LIST) {
  L <- if (!is.null(lats)) lats[[as.character(s)]] else sccb_embed(obj, "scvi", s, batch = "sample"); if (is.null(L)) stop("scVI unavailable")
  lab[[as.character(s)]] <- as.character(sccb_louvain_to_k(L, 3, SEED))
  message("seed ", s, ": ", paste(table(lab[[as.character(s)]]), collapse = "/"))
}
# A: cells x seeds matrix of cluster labels; column 1 is seed 1, the reference partition
A <- do.call(cbind, lab); rownames(A) <- colnames(obj)
write.csv(data.frame(cell_id = rownames(A), A, check.names = FALSE),
          file.path(TABLE_DIR, "13t_scvi_rods_assignments.csv"), row.names = FALSE)
pairs <- combn(names(lab), 2)
message("pairwise ARI: ", paste(round(apply(pairs, 2, function(p) ARI(lab[[p[1]]], lab[[p[2]]])), 3), collapse = " "))

# ---- Profile of the seed-1 clusters ----
cl <- A[, 1]; obj$scvi_k3 <- cl
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
# det(g): fraction of cells in each cluster with at least one count of gene g (NA if g is absent)
det <- function(g) if (g %in% rownames(cnt)) tapply(cnt[g, ] > 0, cl, mean) else rep(NA, length(unique(cl)))
rb  <- grepl("^Rp[sl]", rownames(cnt)); frac_rb <- Matrix::colSums(cnt[rb, ]) / Matrix::colSums(cnt)
prof <- data.frame(cluster = sort(unique(cl)), n = as.numeric(table(cl)),
                   Xist = det("Xist"), Tsix = det("Tsix"), Ddx3y = det("Ddx3y"), Eif2s3y = det("Eif2s3y"),
                   Uty = det("Uty"), Kdm5d = det("Kdm5d"), pct_ribo = round(100 * tapply(frac_rb, cl, mean), 2),
                   median_umi = as.numeric(tapply(Matrix::colSums(cnt), cl, median)))
# Share of each cluster coming from each sample (rows sum to 1)
comp <- as.data.frame.matrix(prop.table(table(cl, obj$sample), 1))
prof <- cbind(prof, round(comp[as.character(prof$cluster), , drop = FALSE], 3))
write.csv(prof, file.path(TABLE_DIR, "13t_scvi_rods_profile.csv"), row.names = FALSE); print(prof)

# ---- Top 20 positive markers per cluster ----
Idents(obj) <- "scvi_k3"
mk <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.2, verbose = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 20)
write.csv(top, file.path(TABLE_DIR, "13t_scvi_rods_markers.csv"), row.names = FALSE)
print(top %>% slice_max(avg_log2FC, n = 10) %>% summarise(genes = paste(gene, collapse = ", ")))

# ---- Overlap with the scDHA Louvain k = 3 rod partition (step 04) ----
# The file keeps its 03_ prefix but holds the rod partitions written by step 04
asg <- read.csv(file.path(TABLE_DIR, "03_all_cluster_assignments.csv"))
sh <- intersect(asg$cell_id, names(cl))
ov <- as.data.frame.matrix(table(scvi = cl[sh], scDHA_LL_k3 = asg$scDHA_LL_k3[match(sh, asg$cell_id)]))
write.csv(ov, file.path(TABLE_DIR, "13t_scvi_rods_overlap.csv")); print(ov)
message(sprintf("ARI scVI seed1 k3 vs scDHA_LL_k3 on %d shared cells: %.3f", length(sh),
                ARI(cl[sh], asg$scDHA_LL_k3[match(sh, asg$cell_id)])))
