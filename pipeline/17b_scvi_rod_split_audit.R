# ============================================================
# pipeline/17b_scvi_rod_split_audit.R
#
# Audits the scVI rod k = 3 split (from 13s/13t; R = 0.84, but one of the five
# embeddings does not reach k = 3, so it fails the validity gates) by sample, using
# tables only. It asks whether the split follows one sample (a batch effect)
# and whether it is the same axis that the single-embedding rod partitions see:
#   (1) share of each scVI cluster within each sample, and ARI/NMI against sample
#   (2) cross-seed ARI of the five scVI k = 3 partitions within each sample
#       (a split driven by batch cannot reproduce within a single sample)
#   (3) concordance with the rod partitions of step 04 (scDHA/PCA Louvain,
#       K-means and GMM, k = 2 and 3)
#
# Paper:    scVI rod split, audit by sample
# Inputs:   output/tables/13t_scvi_rods_assignments.csv (cell x seed, k = 3)
#           output/tables/03_all_cluster_assignments.csv (rod partitions written by 04_phase2_rod_benchmark.R)
# Outputs:  output/tables/17b_scvi_rod_split_by_sample.csv
#           output/tables/17b_scvi_rod_split_concordance.csv
# Run:      Rscript run_pipeline.R 17b   or   Rscript pipeline/17b_scvi_rod_split_audit.R
# Runtime:  under 10 s
# ============================================================
source("config.R"); source("R/utils.R")
# ari(x, y) / nmi(x, y): ARI() / NMI() from R/utils.R (labels converted to factors)
ari <- ARI; nmi <- NMI

# ---- scVI assignments; the sample is the cell ID prefix (e.g. S1_) ----
a <- read.csv(file.path(TABLE_DIR, "13t_scvi_rods_assignments.csv"), check.names = FALSE)
seeds <- setdiff(names(a), "cell_id"); cl <- a[[seeds[1]]]
smp <- sub("_.*", "", a$cell_id); pr <- combn(seeds, 2)
# xs(i): pairwise ARI between the five seed partitions, restricted to the cells selected by i
xs <- function(i) apply(pr, 2, function(p) ari(a[[p[1]]][i], a[[p[2]]][i]))

# ---- (1) and (2): cluster shares and within-sample reproducibility, per sample and overall ----
by <- do.call(rbind, lapply(c(sort(unique(smp)), "all"), function(s) {
  i <- if (s == "all") rep(TRUE, nrow(a)) else smp == s; r <- xs(i)
  f <- prop.table(table(factor(cl[i], sort(unique(cl)))))
  data.frame(sample = s, n = sum(i), frac_c1 = round(f[[1]], 3), frac_c2 = round(f[[2]], 3), frac_c3 = round(f[[3]], 3),
             within_cross_seed_ari_mean = round(mean(r), 3), within_cross_seed_ari_min = round(min(r), 3))
}))
write.csv(by, file.path(TABLE_DIR, "17b_scvi_rod_split_by_sample.csv"), row.names = FALSE); print(by, row.names = FALSE)
message(sprintf("scVI k3 vs sample: ARI %.3f, NMI %.3f", ari(cl, smp), nmi(cl, smp)))

# ---- (3) concordance with the step 04 rod partitions ----
# The file keeps its 03_ prefix but holds the rod partitions written by step 04
d <- read.csv(file.path(TABLE_DIR, "03_all_cluster_assignments.csv"))
m <- match(a$cell_id, d$cell_id); ok <- !is.na(m)
cols <- intersect(c("scDHA_LL_k2", "scDHA_LL_k3", "PCA_LL_k2", "PCA_LL_k3", "PCA_KMeans_k2", "PCA_KMeans_k3",
                    "scDHA_KMeans_k2", "scDHA_KMeans_k3", "PCA_GMM_k2"), names(d))
# Two-group version of the scVI split: the second-largest seed-1 cluster is the minority group,
# the other clusters the majority (Louvain label numbers are arbitrary, so sizes decide)
minor_lab <- names(sort(table(cl), decreasing = TRUE))[2]
axis <- ifelse(cl == minor_lab, "minor", "major")
conc <- data.frame(partition = cols, n_cells = sum(ok),
                   ari_vs_scvi_k3 = sapply(cols, function(p) round(ari(cl[ok], d[[p]][m[ok]]), 3)),
                   ari_vs_scvi_major_minor = sapply(cols, function(p) round(ari(axis[ok], d[[p]][m[ok]]), 3)))
write.csv(conc, file.path(TABLE_DIR, "17b_scvi_rod_split_concordance.csv"), row.names = FALSE); print(conc, row.names = FALSE)
