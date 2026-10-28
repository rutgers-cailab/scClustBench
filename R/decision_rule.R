# ============================================================
# R/decision_rule.R
#
# The scClustBench decision rule. Both stages use label-free quantities only:
# Stage A screens and ranks candidate partitions from the metric benchmark, and
# Stage B decides from the across-embedding reproducibility curve R(k) whether a
# population has reproducible sub-structure, and at what granularity k*.
# Sourced by R/scClustBench.R and by the pipeline scripts (step 12 applies the
# rule to every population and embedding family).
#
#  Stage A  sccb_select(scores, jacc_min, min_cells, k)
#             A1 validity : achieved k == target k, and no cluster smaller
#                           than min_cells (drops degenerate singleton splits)
#             A2 stability: min per-cluster bootstrap Jaccard >= jacc_min
#                           (Hennig 2007: < 0.6 = unstable)
#             A3 ranking  : mean percentile rank over the metrics in SCCB_METRICS
#                           (ROGUE, balance, Jaccard, marker specificity,
#                           silhouette), equal weights; ties -> smaller k.
#           Ranking across k is biased toward small k (Jaccard, silhouette and
#           balance all rise as clusters merge), so Stage A does not choose k.
#           With `k = k*` from Stage B it ranks configurations at that k only,
#           which is how the reported partition is chosen.
#
#  Stage B  sccb_verdict(curve, tau, tol, min_cells)
#           R(k) = mean pairwise ARI between Louvain partitions of independently
#           trained or perturbed embeddings (from sccb_reproducibility()).
#             B1 validity : the Stage A validity gate applied at each k: achieved
#                           k == k for the reference (seed-1) partition and, when
#                           the curve stores n_embeddings_at_k, for every
#                           embedding; and, when it stores seed1_sizes, no
#                           reference cluster smaller than min_cells. Invalid k
#                           are excluded before the maximum is taken.
#             B2 ACCEPT if max over valid k of R(k) >= tau, else REJECT.
#             B3 k* = the largest valid k whose R(k) is within `tol` of that
#                maximum (merging reproducible clusters keeps them
#                reproducible, so the smallest such k is biased coarse).
#                k_star_smallest is reported alongside; k_star_at_edge flags a
#                k* on the last k of the sweep (k* is then a lower bound).
#                k_above_tau lists every valid k with R(k) >= tau; k* is a
#                granularity summary, not an estimate of the number of types.
#           Trend diagnostics (Spearman rho of R and silhouette vs k) are
#           reported but do not gate the decision.
#
#  Families sccb_tau(embedding, calib, tau) gives max(tau, the family's null 95th
#           percentile): calibration can only raise tau, never lower it.
#           sccb_consensus(verdicts) combines embedding families: ACCEPT or
#           REJECT only if every family agrees, otherwise EMBEDDING-DEPENDENT.
#
# To change a threshold, pass it as an argument (defaults: jacc_min = 0.6,
# min_cells = 20, tau = 0.8, tol = 0.01). To change the ranking metrics, edit
# SCCB_METRICS (each must be a column of the benchmark table).
#
# Run:      Rscript R/decision_rule.R   (runs the built-in self-check at the end of this file)
# ============================================================


# Metrics averaged in the Stage A ranking (higher is better for each)
SCCB_METRICS <- c("rogue", "balance", "jacc_min", "marker_spec", "silhouette")

# .min_size(s): smallest cluster size in a sizes string such as "250/120/30"
.min_size <- function(s) min(as.numeric(strsplit(gsub("[ ,]", "", s), "/")[[1]]))

# ---- Stage A ----
# sccb_select(scores, jacc_min, min_cells, k): takes the benchmark table (one row per method and target_k)
# and returns it with columns smallest, valid, stable, admissible, score (mean percentile, admissible
# rows only) and rank, sorted best first. With k given, only rows with target_k in k are kept.
sccb_select <- function(scores, jacc_min = 0.6, min_cells = 20, k = NULL) {
  s <- if (is.null(k)) scores else scores[scores$target_k %in% k, ]
  s$smallest <- if ("sizes_str" %in% names(s)) vapply(s$sizes_str, .min_size, numeric(1)) else Inf
  s$valid  <- !is.na(s$actual_k) & s$actual_k == s$target_k & s$smallest >= min_cells
  s$stable <- !is.na(s$jacc_min) & s$jacc_min >= jacc_min
  s$admissible <- s$valid & s$stable
  s$score <- NA_real_
  a <- s$admissible
  if (any(a)) {
    # Percentile rank of each metric among admissible rows, scaled to [0, 1]; NA metrics are skipped
    pr <- sapply(SCCB_METRICS, function(m) {
      v <- s[[m]][a]; r <- rank(v, na.last = "keep", ties.method = "average")
      if (sum(!is.na(v)) > 1) (r - 1) / (sum(!is.na(v)) - 1) else rep(1, length(v))
    })
    s$score[a] <- rowMeans(matrix(pr, nrow = sum(a)), na.rm = TRUE)
  }
  # Ties go to the smaller k; scores are compared after rounding so floating-point noise cannot break a tie
  s <- s[order(-s$admissible, -round(s$score, 10), s$target_k), ]
  s$rank <- ifelse(s$admissible, seq_len(nrow(s)), NA_integer_)
  rownames(s) <- NULL
  s
}

# sccb_valid_k(curve, min_cells): logical vector, TRUE for each row (k) of an R(k) curve that passes the
# Stage B validity gates: the reference partition reached k (actual_k), every embedding reached k
# (n_embeddings_at_k, or n_seeds_at_k in some saved curves) and the reference partition has no cluster
# smaller than min_cells (seed1_sizes). A gate whose column is absent is not applied.
sccb_valid_k <- function(curve, min_cells = 20) {
  ok <- is.na(curve$actual_k) | curve$actual_k == curve$k                      # drop k the reference never reached
  if ("seed1_sizes" %in% names(curve)) {                                        # degenerate clusters
    sz <- suppressWarnings(vapply(as.character(curve$seed1_sizes), .min_size, numeric(1)))
    ok <- ok & (is.na(sz) | sz >= min_cells) }
  ne_col <- intersect(c("n_embeddings_at_k", "n_seeds_at_k"), names(curve))
  # The family size is not stored in the curve, so it is taken as the largest number of embeddings
  # that reached any k
  if (length(ne_col)) { ne <- curve[[ne_col[1]]]; ok <- ok & ne == max(ne, na.rm = TRUE) }
  ok
}

# ---- Stage B ----
# sccb_verdict(curve, tau, tol, min_cells): takes an R(k) curve, data.frame(k, actual_k,
# mean_cross_seed_ari[, silhouette, seed1_sizes, n_embeddings_at_k]), and returns a one-row data.frame
# with the verdict, k*, the gated and ungated max R, the valid k above tau and the trend diagnostics.
sccb_verdict <- function(curve, tau = 0.8, tol = 0.01, min_cells = 20) {
  ok <- sccb_valid_k(curve, min_cells)
  cv <- curve[ok, ]
  R <- cv$mean_cross_seed_ari
  mR <- if (length(R)) max(R, na.rm = TRUE) else NA_real_
  hit <- if (!is.na(mR) && mR >= tau) cv$k[R >= mR - tol] else integer(0)
  rho <- function(y) suppressWarnings(cor(cv$k, y, method = "spearman", use = "complete.obs"))
  data.frame(
    verdict    = if (length(hit)) "ACCEPT" else "REJECT",
    k_star     = if (length(hit)) max(hit) else NA_integer_,
    k_star_smallest = if (length(hit)) min(hit) else NA_integer_,
    k_star_at_edge  = length(hit) > 0 && max(hit) == max(curve$k),
    max_R      = mR,
    max_R_ungated = max(curve$mean_cross_seed_ari, na.rm = TRUE),
    n_valid_k  = nrow(cv),
    k_first_tau = if (any(R >= tau)) min(cv$k[R >= tau]) else NA_integer_,
    k_above_tau = paste(cv$k[!is.na(R) & R >= tau], collapse = ","),     # every valid k with R >= tau
    rho_R_k    = if (nrow(cv) > 2) round(rho(R), 3) else NA_real_,
    rho_sil_k  = if ("silhouette" %in% names(cv) && nrow(cv) > 2) round(rho(cv$silhouette), 3) else NA_real_,
    tau        = tau)
}

# ---- Embedding families ----
# sccb_tau(embedding, calib, tau): returns tau for one embedding family: the fixed tau, raised (never
# lowered) to the family's null 95th percentile when a calibration table (embedding, tau_calibrated)
# is supplied, e.g. from step 16.
sccb_tau <- function(embedding, calib = NULL, tau = 0.8) {
  q <- if (is.null(calib)) NA else calib$tau_calibrated[match(embedding, calib$embedding)]
  max(tau, q, na.rm = TRUE)
}

# sccb_consensus(verdicts): combines per-family verdicts into ACCEPT, REJECT or EMBEDDING-DEPENDENT
sccb_consensus <- function(verdicts) {
  if (all(verdicts == "ACCEPT")) "ACCEPT" else if (all(verdicts == "REJECT")) "REJECT" else "EMBEDDING-DEPENDENT"
}

# ---- Self-check (runs only when the file is executed directly: Rscript R/decision_rule.R) ----
if (sys.nframe() == 0) {
  g <- data.frame(method = c("a","a","b"), target_k = c(2,3,3), actual_k = c(2,3,4),
                  rogue = c(.8,.9,.9), balance = c(.5,.6,.9), jacc_min = c(.9,.7,.9),
                  marker_spec = c(.1,.2,.3), silhouette = c(.3,.2,.5),
                  sizes_str = c("50/50","40/40/20","25/25/25/25"))
  # Stage A: "b" overshoots its target k, so it is not admissible; a lower Jaccard fails A2
  sel <- sccb_select(g)
  stopifnot(sel$method[1] == "a", sel$target_k[1] == 3, !sel$admissible[sel$method == "b"])
  stopifnot(!sccb_select(transform(g, jacc_min = .5))$admissible[1])
  # Stage B: ACCEPT at the largest k above tau; REJECT when no k reaches tau
  v <- sccb_verdict(data.frame(k = 2:4, actual_k = 2:4, mean_cross_seed_ari = c(.5,.85,.9), silhouette = 1:3))
  stopifnot(v$verdict == "ACCEPT", v$k_star == 4)
  stopifnot(sccb_verdict(data.frame(k = 2:3, actual_k = 2:3, mean_cross_seed_ari = c(.7,.6)))$verdict == "REJECT")
  # B1: a k whose reference partition has a 2-cell cluster cannot carry the verdict
  cz <- data.frame(k = 2:4, actual_k = 2:4, mean_cross_seed_ari = c(.5, .7, .85),
                   seed1_sizes = c("60/40", "50/30/20", "50/30/18/2"))
  vz <- sccb_verdict(cz)
  stopifnot(vz$verdict == "REJECT", vz$max_R == .7, vz$max_R_ungated == .85, vz$n_valid_k == 2)
  stopifnot(sccb_verdict(cz, min_cells = 1)$k_star == 4)
  stopifnot(sccb_verdict(cz[3, ])$verdict == "REJECT")                     # no valid k at all
  # B1: a k reached by only some embeddings is excluded
  cn <- data.frame(k = 2:4, actual_k = 2:4, mean_cross_seed_ari = c(.5, .82, .9), n_embeddings_at_k = c(5, 5, 4))
  stopifnot(sccb_verdict(cn)$k_star == 3, sccb_verdict(cn[, -4])$k_star == 4)
  # B3: largest k within tol; smallest reported
  v2 <- sccb_verdict(data.frame(k = 2:5, actual_k = 2:5, mean_cross_seed_ari = c(.95, .6, .945, .5)))
  stopifnot(v2$k_star == 4, v2$k_star_smallest == 2, !v2$k_star_at_edge, v2$k_above_tau == "2,4")
  # tau calibration only raises tau; consensus
  cal <- data.frame(embedding = c("scdha", "pca_perturb"), tau_calibrated = c(.85, .12))
  stopifnot(sccb_tau("scdha", cal) == .85, sccb_tau("pca_perturb", cal) == .8, sccb_tau("scvi", cal) == .8)
  stopifnot(sccb_consensus(c("ACCEPT", "ACCEPT")) == "ACCEPT",
            sccb_consensus(c("ACCEPT", "REJECT")) == "EMBEDDING-DEPENDENT")
  # Stage A ranking restricted to one k
  stopifnot(all(sccb_select(g, k = 2)$target_k == 2))
  cat("decision_rule.R self-check OK\n")
}
