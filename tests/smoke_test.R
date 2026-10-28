# ============================================================
# tests/smoke_test.R
#
# End-to-end check of R/scClustBench.R on simulated counts, with no scDHA or
# ROGUE needed. Two datasets are run through sccb_run() and the verdicts checked:
#   Structured: 3 groups x 300 cells with group-specific DE genes, 2 batches -> expect ACCEPT, k* = 3
#   Null:       900 cells, one population, same 2 batches                 -> expect REJECT
#
# Outputs:  CSVs in tempdir()/structured and tempdir()/null (removed when R exits)
# Run:      Rscript tests/smoke_test.R
# Runtime:  about 10 minutes on a laptop
# ============================================================
source("R/scClustBench.R")
suppressPackageStartupMessages(library(Matrix))
# ---- Simulated data ----
# sim(structured, n_per, n_genes, seed): returns a Seurat object of negative binomial counts for 3 x n_per
# cells in two batches. With structured = TRUE each group has 60 up-regulated genes; otherwise all
# cells share one expression profile. Metadata: batch, and the true group (used only to score ARI).
sim <- function(structured, n_per = 300, n_genes = 2000, seed = 11) {
  set.seed(seed)
  n <- 3 * n_per; grp <- rep(1:3, each = n_per); bat <- rep(c("b1", "b2"), length.out = n)
  mu <- matrix(rexp(n_genes, 1 / 0.5), n_genes, n)                  # baseline gene means
  if (structured) for (g in 1:3) mu[(g - 1) * 60 + 1:60, grp == g] <- mu[(g - 1) * 60 + 1:60, grp == g] * 5
  mu[1500:1600, bat == "b2"] <- mu[1500:1600, bat == "b2"] * 1.5   # mild batch effect
  mu <- sweep(mu, 2, rlnorm(n, 0, 0.3), "*")                        # library-size variation
  x <- matrix(rnbinom(length(mu), mu = mu, size = 2), n_genes, n,
              dimnames = list(paste0("g", 1:n_genes), paste0("c", 1:n)))
  o <- CreateSeuratObject(Matrix(x, sparse = TRUE)); o$batch <- bat; o$group <- c("alpha", "beta", "gamma")[grp]; o   # character labels exercise the aricode factor conversion
}
# ---- Run the full rule on both datasets ----
# run(o, tag): sccb_run() with a small grid (k = 2..5, 2 bootstraps), pca_perturb as the Stage B family
run <- function(o, tag) sccb_run(o, batch = "batch", k_range = 2:5, embeddings = c("pca", "pca_perturb"),
                                 methods = c("kmeans", "gmm", "louvain", "leiden"), repro_embeddings = "pca_perturb",
                                 n_boot = 2, truth = "group", out_dir = file.path(tempdir(), tag))
s <- run(sim(TRUE), "structured"); n <- run(sim(FALSE), "null")

# ---- Print summaries ----
cat("\n--- structured: Stage A top 3 ---\n"); print(head(s$selection[, c("rank", "method", "target_k", "jacc_min", "silhouette", "score", "ari")], 3))
cat("--- structured: R(k) ---\n"); print(s$curve[, c("k", "actual_k", "mean_cross_seed_ari", "ari_vs_truth")])
cat("--- null: R(k) ---\n");       print(n$curve[, c("k", "actual_k", "mean_cross_seed_ari")])
cat("--- recommended ---\n"); print(rbind(s$recommended, n$recommended))
cat("--- verdicts ---\n"); print(rbind(cbind(data = "structured", s$verdict), cbind(data = "null", n$verdict)))

# ---- Assertions ----
# Output format: every metric column is present and the CSVs were written
stopifnot(all(c("method", "target_k", "actual_k", "balance", "jacc_min", "marker_spec", "rogue", "silhouette", "sizes_str")
              %in% names(s$scores)))
stopifnot(file.exists(file.path(tempdir(), "structured", "sccb_verdict.csv")))
# Structured data: ACCEPT at the true k, and Stage A ranks a k = 3 partition that matches the truth first
stopifnot(s$verdict$verdict == "ACCEPT", s$verdict$k_star == 3)
stopifnot(s$selection$target_k[1] == 3, s$selection$ari[1] > 0.95)
# Null data: no reproducible sub-structure
stopifnot(n$verdict$verdict == "REJECT")
# The reported answer: configuration chosen at k*, with labels; nothing reported for the null
stopifnot(s$recommended$k_star == 3, !is.na(s$recommended$method), nrow(s$partition) == 900,
          aricode::ARI(factor(s$partition$cluster), factor(sim(TRUE)$group)) > 0.95,
          "min_cluster_size" %in% names(s$curve))
stopifnot(n$recommended$verdict == "REJECT", is.null(n$partition))
cat("\nsmoke_test.R: ALL ASSERTIONS PASSED\n")
