# ============================================================
# pipeline/12c_recompute_jaccard.R
#
# Computes the Stage A stability column (jacc_min: the minimum over clusters of the mean best-match
# bootstrap Jaccard) for every cached (method, k) partition of Phase 1 (step 03) and the rods
# (step 04), and writes it into their 4-metric tables. It uses .jacc_min() from R/scClustBench.R,
# which draws all 80% subsamples before clustering so that each bootstrap uses a different subsample
# (the clustering functions reseed the RNG internally). Nothing else is re-clustered; the value
# written by steps 03/04 is kept as jacc_min_old. Run pipeline/12_decision_rule.R and the figures
# afterwards.
#
# Paper:    Stage A admissibility (min bootstrap Jaccard >= 0.6) for Phase 1 and the rods
# Inputs:   output_full_dataset/tables/poc_harmony_{pca,scdha}_latent.csv, poc_all_cluster_assignments.csv,
#           poc_scores_matrix_4metric.csv                      (DATASET=phase1, from step 03)
#           output/tables/03_harmony_{pca,scdha}_latent.csv, 03_all_cluster_assignments.csv,
#           03_scores_matrix_4metric.csv                       (DATASET=rods, from step 04)
# Outputs:  the 4-metric table(s) above, rewritten in place (jacc_min updated, jacc_min_old added);
#           with ONLY set: <tables dir>/12c_parts/<method>_k<k>.csv
# Run:      DATASET=phase1 Rscript pipeline/12c_recompute_jaccard.R   (run_pipeline.R steps 12c_p1, 12c_rods)
# Options:  N_BOOT      bootstrap subsamples per partition (default 20)
#           DATASET     phase1 | rods | both (default both)
#           RESUME_LOG  log of an interrupted run; partitions it already recomputed are skipped and keep
#                       the value in the checkpointed table (optional)
#           ONLY        one partition, e.g. PCA_LL_k7: compute only it and write 12c_parts/<col>.csv
#                       without rewriting the table (for Slurm array jobs, hpc/slurm/jacc_rods.sbatch)
#           MERGE=1     apply every 12c_parts/*.csv to the 4-metric table(s) and exit
# Runtime:  dominated by the resolution search of the graph methods; roughly 3-5 h for the rods and
#           6-10 h for Phase 1 at N_BOOT=20 (N_BOOT=10 halves it). About 32 GB RAM.
# ============================================================
source("config.R"); source("R/utils.R"); source("R/decision_rule.R"); source("R/scClustBench.R")
N_BOOT <- as.integer(Sys.getenv("N_BOOT", "20")); DATASET <- Sys.getenv("DATASET", "both")
# Maps the method suffix in the table (e.g. "LL" in PCA_LL) to its SCCB_METHODS clustering function
METHOD_FN <- c(KMeans = "kmeans", GMM = "gmm", LL = "louvain", Leiden = "leiden")

# Takes a CSV path (cell_id column, then one column per dimension); returns a cells x dims matrix.
read_lat <- function(f) { d <- read.csv(f, check.names = FALSE); m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; m }
# File locations and prefixes of the two datasets
sets <- list(
  phase1 = list(dir = "output_full_dataset/tables", pre = "poc_", four = "poc_scores_matrix_4metric.csv"),
  rods   = list(dir = TABLE_DIR, pre = "03_", four = "03_scores_matrix_4metric.csv"))
if (DATASET != "both") sets <- sets[DATASET]

# Partitions already recomputed in the RESUME_LOG run (lines of the form "  <col>: x (old y)")
done <- if (nzchar(Sys.getenv("RESUME_LOG"))) sub("^\\s+(\\S+):.*", "\\1",
  grep("\\(old ", readLines(Sys.getenv("RESUME_LOG")), value = TRUE)) else character(0)
if (length(done)) message("Resuming: skipping ", length(done), " partitions already recomputed")

for (nm in names(sets)) {
  s <- sets[[nm]]; message("== ", nm)
  lat <- list(PCA = read_lat(file.path(s$dir, paste0(s$pre, "harmony_pca_latent.csv"))),
              scDHA = read_lat(file.path(s$dir, paste0(s$pre, "harmony_scdha_latent.csv"))))
  asg <- read.csv(file.path(s$dir, paste0(s$pre, "all_cluster_assignments.csv")))
  four_f <- file.path(s$dir, s$four); four <- read.csv(four_f)
  if (!"jacc_min_old" %in% names(four)) four$jacc_min_old <- four$jacc_min
  parts <- file.path(s$dir, "12c_parts")
  # MERGE mode: collect the per-partition results of an ONLY array run into the table
  if (nzchar(Sys.getenv("MERGE"))) {
    for (f in list.files(parts, "\\.csv$", full.names = TRUE)) {
      p <- read.csv(f); i <- which(paste0(four$method, "_k", four$target_k) == p$col)
      four$jacc_min[i] <- p$jacc_min; message(sprintf("  merged %s: %.3f", p$col, p$jacc_min))
    }
    write.csv(four, four_f, row.names = FALSE); next
  }
  # Skip partitions without assignments, failed partitions (actual_k NA) and resumed ones
  for (i in seq_len(nrow(four))) {
    col <- paste0(four$method[i], "_k", four$target_k[i])
    if (!col %in% names(asg) || is.na(four$actual_k[i]) || col %in% done) next
    if (nzchar(Sys.getenv("ONLY")) && col != Sys.getenv("ONLY")) next
    emb <- sub("_.*", "", four$method[i]); L <- lat[[emb]][asg$cell_id, , drop = FALSE]
    # Steps 03/04 build the PCA + Leiden graph on the first N_DIMS_NEIGHBORS PCs; re-cluster on the same dims
    if (four$method[i] == "PCA_Leiden") L <- L[, seq_len(N_DIMS_NEIGHBORS), drop = FALSE]
    # Bootstrap at actual_k (the k the partition reached), comparing against the cached partition
    fn <- SCCB_METHODS[[METHOD_FN[[sub(".*_", "", four$method[i])]]]]
    four$jacc_min[i] <- round(tryCatch(.jacc_min(fn, L, asg[[col]], four$actual_k[i], N_BOOT, SEED),
                                       error = function(e) NA_real_), 4)
    message(sprintf("  %s: %.3f (old %.3f)", col, four$jacc_min[i], four$jacc_min_old[i]))
    if (nzchar(Sys.getenv("ONLY"))) { dir.create(parts, showWarnings = FALSE)
      write.csv(data.frame(col = col, jacc_min = four$jacc_min[i]), file.path(parts, paste0(col, ".csv")), row.names = FALSE)
    } else write.csv(four, four_f, row.names = FALSE)                    # checkpoint
  }
}
message("Done. Next: Rscript run_pipeline.R 12 06")
