# ============================================================
# pipeline/04e_knn_robustness.R
#
# Sensitivity of the illustrative rod partition (scDHA + Louvain) to the number of neighbours (k.param)
# used to build the SNN graph. For each k in STABILITY_K_RANGE, the rod scDHA latent is clustered with
# Louvain at k.param = 10, 20, 30, 50 and 80, and each partition is compared with the one at
# k.param = KNN_K (20).
#
# Paper:    k.param sensitivity of the illustrative rod partition (supporting table)
# Inputs:   output/tables/03_harmony_scdha_latent.csv (step 04)
# Outputs:  output/tables/_knn_robustness_summary.csv. Columns: cluster sizes, smallest cluster, number of
#           clusters under 20 cells, and overlap_pct_vs_orig_kparam20 = % of cells in the cluster that
#           best overlaps a k.param = KNN_K cluster (max-overlap matching per cluster).
# Run:      Rscript run_pipeline.R 04e
# Runtime:  about 20 minutes (5 k.param values x 3 k, each a Louvain resolution search)
# ============================================================
source("config.R"); source("R/scClustBench.R")
df  <- read.csv(file.path(TABLE_DIR, "03_harmony_scdha_latent.csv"), check.names = FALSE)
lat <- as.matrix(df[, -1]); rownames(lat) <- as.character(df[[1]])
# Takes cluster sizes; returns them as one " / "-separated string with thousands separators
fmt <- function(sz) paste(formatC(sz, big.mark = ",", width = 6), collapse = " / ")
rows <- list()
for (k in STABILITY_K_RANGE) {
  # reference partition at the default neighbour count; k.param is passed through to FindNeighbors
  ref <- sccb_louvain_to_k(lat, k, SEED, k.param = KNN_K)
  for (kp in c(10, 20, 30, 50, 80)) {
    lab <- if (kp == KNN_K) ref else sccb_louvain_to_k(lat, k, SEED, k.param = kp)
    sz <- sort(as.integer(table(lab)), decreasing = TRUE)
    rows[[length(rows) + 1]] <- data.frame(k_param = kp, target_k = k, actual_k = length(sz), sizes = fmt(sz),
      smallest_cluster = min(sz), n_tiny_clusters = sum(sz < 20),
      overlap_pct_vs_orig_kparam20 = round(100 * sum(apply(table(lab, ref), 1, max)) / length(lab), 2))
    print(tail(rows, 1)[[1]])
  }
}
out <- do.call(rbind, rows); out <- out[order(out$k_param, out$target_k), ]
write.csv(out, file.path(TABLE_DIR, "_knn_robustness_summary.csv"), row.names = FALSE)
