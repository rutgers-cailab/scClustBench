# ============================================================
# pipeline/12d_rescore_phase1_ari.R
#
# Recomputes the Phase 1 ARI and NMI of every cached (method, k) partition against the reference
# annotation in output/rds/02_annotated.rds (obj$cell_type), without re-clustering. This keeps the
# step 03 tables consistent with the annotation in data/annotation/ even when step 03 is not rerun.
# The values found in the tables are kept as ari_old / nmi_old, and a per-partition comparison is written.
#
# Paper:    Phase 1 agreement with the reference annotation (ARI, NMI)
# Inputs:   output/rds/02_annotated.rds
#           output_full_dataset/tables/poc_all_cluster_assignments.csv
#           output_full_dataset/tables/poc_scores_matrix.csv, poc_scores_matrix_4metric.csv
# Outputs:  output_full_dataset/tables/poc_scores_matrix.csv, poc_scores_matrix_4metric.csv
#           (ari, nmi rewritten in place)
#           output_full_dataset/tables/12d_phase1_ari_rescore.csv
# Run:      Rscript run_pipeline.R 12d
# Runtime:  about 2 minutes; about 32 GB RAM (reads the annotated object once)
# ============================================================
source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({ library(Seurat); library(aricode) })
D <- "output_full_dataset/tables"
# Keep only the reference labels and free the object
obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds")); truth <- setNames(as.character(obj$cell_type), colnames(obj)); rm(obj); gc()
asg <- read.csv(file.path(D, "poc_all_cluster_assignments.csv")); stopifnot(all(asg$cell_id %in% names(truth)))
tr <- as.integer(factor(truth[asg$cell_id])); print(table(truth[asg$cell_id]))
# ARI and NMI per partition column, over cells with a cluster label
res <- do.call(rbind, lapply(setdiff(names(asg), "cell_id"), function(col) {
  ok <- !is.na(asg[[col]]); l <- as.integer(factor(asg[[col]][ok]))
  data.frame(col = col, ari = round(ARI(l, tr[ok]), 4), nmi = round(NMI(l, tr[ok]), 4)) }))
# Write the recomputed values into both score tables, keyed by <method>_k<target_k>
for (f in c("poc_scores_matrix.csv", "poc_scores_matrix_4metric.csv")) {
  t <- read.csv(file.path(D, f)); key <- paste0(t$method, "_k", t$target_k); m <- match(key, res$col)
  if (!"ari_old" %in% names(t)) { t$ari_old <- t$ari; t$nmi_old <- t$nmi }
  t$ari[!is.na(m)] <- res$ari[m[!is.na(m)]]; t$nmi[!is.na(m)] <- res$nmi[m[!is.na(m)]]
  write.csv(t, file.path(D, f), row.names = FALSE)
}
# Per-partition comparison of recomputed and stored (ari_old) ARI; print a few partitions discussed in the paper
out <- merge(res, transform(read.csv(file.path(D, "poc_scores_matrix.csv")), col = paste0(method, "_k", target_k))[, c("col", "ari_old", "nmi_old")])
out$d_ari <- round(out$ari - out$ari_old, 4); write.csv(out, file.path(D, "12d_phase1_ari_rescore.csv"), row.names = FALSE)
print(out[out$col %in% c("scDHA_LL_k6", "scDHA_KMeans_k6", "PCA_LL_k6", "PCA_LL_k8", "PCA_Leiden_k6", "PCA_KMeans_k5", "PCA_KMeans_k8"), ])
cat("partitions with |d_ari| > 0.01:", sum(abs(out$d_ari) > 0.01, na.rm = TRUE), "of", nrow(out), "\n")
