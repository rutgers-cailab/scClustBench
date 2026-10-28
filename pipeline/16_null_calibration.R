# ============================================================
# pipeline/16_null_calibration.R
#
# Calibrates the Stage B threshold tau on Splatter simulations. Null datasets
# contain one population, so any accepted split is a false positive; structured
# datasets contain three groups, so they measure sensitivity and whether k* = 3.
# R(k) for k = 2 to 6 is computed with the scDHA family (if installed) and the
# PCA gene-resampling family, recording the achieved k. The calibrated tau per
# family is the 95th percentile of max_k R(k) over the null datasets (with 20
# nulls this is close to the second-largest null value). The script then
# reports sensitivity, specificity and correct k* of sccb_verdict() at tau in
# {0.75, 0.80, 0.85, calibrated}. It complements the synthetic validation of
# step 09 with more null datasets and the achieved-k record.
#
#   NULL        20 single-population datasets, nCells 2000 to 10000,
#               library-size location 10 to 12, distinct seeds
#   STRUCTURED  de.prob {0.02, 0.05, 0.10, 0.20} x groups
#               {balanced 3, or 3 with a 1% / 2% / 5% rare group}
#
# Paper:    Null calibration of tau (Splatter simulations)
# Inputs:   none (simulated); needs Bioconductor splatter
# Outputs:  output/tables/16_null_calibration_curves.csv      (also the resume checkpoint)
#           output/tables/16_null_calibration_tau.csv
#           output/tables/16_null_calibration_performance.csv
# Run:      Rscript run_pipeline.R 16   or   Rscript pipeline/16_null_calibration.R
# Options:  FRESH=1  ignore the curves checkpoint and recompute every dataset (default: resume)
# Runtime:  about 3 h with PCA gene resampling only; about 15 h more with scDHA
#           (36 datasets x 5 trainings). RAM about 16 GB. The curves table is also the
#           resume checkpoint: a rerun skips every dataset x family already in it. The
#           repository includes the paper's curves table, so without FRESH=1 the step only
#           recomputes the tau and performance tables from it.
# ============================================================
if (!requireNamespace("splatter", quietly = TRUE))
  stop('Needs Bioconductor splatter: BiocManager::install("splatter")')
source("config.R"); source("R/scClustBench.R")
EMB <- c(if (requireNamespace("scDHA", quietly = TRUE)) "scdha", "pca_perturb")
KS <- 2:6; NGENES <- 2000
ck <- file.path(TABLE_DIR, "16_null_calibration_curves.csv")

# ---- Simulation designs ----
# Each design is a list(id, type, true_k, n cells, lib location, de prob, group probs, seed)
set.seed(16)
nulls <- lapply(1:20, function(i) list(id = sprintf("null_%02d", i), type = "null", true_k = 1,
  n = sample(seq(2000, 10000, 500), 1), lib = runif(1, 10, 12), de = 0, gp = 1, seed = 1600 + i))
grid <- expand.grid(de = c(0.02, 0.05, 0.1, 0.2), rare = c(NA, 0.01, 0.02, 0.05))
# Rare designs: group proportions (0.85 - rare, 0.15, rare)
strs <- lapply(seq_len(nrow(grid)), function(i) { g <- grid[i, ]
  list(id = sprintf("de%.2f_%s", g$de, if (is.na(g$rare)) "bal" else sprintf("rare%.2f", g$rare)),
       type = "structured", true_k = 3, n = 5000, lib = 11, de = g$de, rare = g$rare,
       gp = if (is.na(g$rare)) c(.34, .33, .33) else c(0.85 - g$rare, 0.15, g$rare), seed = 1700 + i) })

# ---- Simulate and compute R(k), checkpointing after every dataset x family ----
done <- if (file.exists(ck) && Sys.getenv("FRESH") != "1") read.csv(ck) else NULL
# have(id, e): TRUE if dataset id has already been run with embedding family e
have <- function(id, e) !is.null(done) && any(done$dataset == id & done$embedding == e)
for (sc in c(nulls, strs)) {
  if (all(vapply(EMB, function(e) have(sc$id, e), TRUE))) next
  message("=== ", sc$id, " (n=", sc$n, ")")
  p <- splatter::newSplatParams(nGenes = NGENES, batchCells = sc$n, lib.loc = sc$lib, seed = sc$seed)
  sim <- if (sc$type == "null") splatter::splatSimulate(p, verbose = FALSE) else
    splatter::splatSimulateGroups(p, group.prob = sc$gp, de.prob = sc$de, verbose = FALSE)
  m <- SingleCellExperiment::counts(sim); colnames(m) <- paste0("c", seq_len(ncol(m)))
  o <- sccb_prep(CreateSeuratObject(as(m, "dgCMatrix")), "LogNormalize", 1000)
  truth <- if (sc$type == "null") rep("G1", ncol(m)) else as.character(sim$Group)
  for (e in EMB) {
    if (have(sc$id, e)) next
    cv <- sccb_reproducibility(o, e, SCDHA_SEED_LIST, KS, truth = truth)
    if (is.null(cv)) next
    row <- cbind(dataset = sc$id, type = sc$type, true_k = sc$true_k, n_cells = sc$n, lib_loc = round(sc$lib, 2),
                 de_prob = sc$de, rare_prop = if (is.null(sc$rare)) NA else sc$rare, embedding = e, cv)
    done <- rbind(done, row); write.csv(done, ck, row.names = FALSE)
  }
}

# ---- Calibrated tau: 95th percentile of the null max R(k), per family ----
keys <- split(done, list(done$dataset, done$embedding), drop = TRUE)
# null_q95(drop): 95th percentile and maximum of the null max R(k) per family, after removing the gate
# columns in `drop` (tau = Inf, so sccb_verdict() only reports max_R over the k it keeps)
null_q95 <- function(drop = character(0)) {
  m <- do.call(rbind, lapply(keys, function(d) { d[intersect(drop, names(d))] <- NULL
    cbind(d[1, c("type", "embedding")], max_R = sccb_verdict(d, tau = Inf)$max_R) }))
  m <- m[m$type == "null", ]
  merge(aggregate(max_R ~ embedding, m, function(x) quantile(x, 0.95, names = FALSE)),
        aggregate(max_R ~ embedding, m, max), by = "embedding", suffixes = c("_q95", "_max"))
}
# tau_calibrated uses the final rule (every Stage B validity gate); the second pair of columns is the same
# calibration without the every-embedding-reaches-k gate, reported in Supplementary Table S9
f <- null_q95(); u <- null_q95(c("n_embeddings_at_k", "n_seeds_at_k"))
tau <- data.frame(embedding = f$embedding, tau_calibrated = f$max_R_q95, null_max_R_max = f$max_R_max,
                  tau_without_embedding_gate = u$max_R_q95[match(f$embedding, u$embedding)],
                  null_max_R_max_without_embedding_gate = u$max_R_max[match(f$embedding, u$embedding)])
write.csv(tau, file.path(TABLE_DIR, "16_null_calibration_tau.csv"), row.names = FALSE); print(tau)

# ---- Operating characteristics of the rule at fixed and calibrated tau ----
# sensitivity = share of structured datasets accepted; specificity = share of nulls rejected;
# k_star_correct = share of accepted structured datasets with k* = 3
perf <- do.call(rbind, lapply(tau$embedding, function(e) {
  tc <- tau$tau_calibrated[tau$embedding == e]
  do.call(rbind, lapply(c(0.75, 0.8, 0.85, tc), function(t) {
    v <- do.call(rbind, lapply(keys[vapply(keys, function(d) d$embedding[1] == e, TRUE)], function(d)
      cbind(d[1, c("type", "true_k")], sccb_verdict(d, tau = t))))
    s <- v[v$type == "structured", ]
    data.frame(embedding = e, tau = round(t, 4), calibrated = t == tc,
               sensitivity = mean(s$verdict == "ACCEPT"), specificity = mean(v$verdict[v$type == "null"] == "REJECT"),
               k_star_correct = mean(s$k_star[s$verdict == "ACCEPT"] == s$true_k[s$verdict == "ACCEPT"]),
               n_null = sum(v$type == "null"), n_structured = nrow(s))
  }))
}))
write.csv(perf, file.path(TABLE_DIR, "16_null_calibration_performance.csv"), row.names = FALSE); print(perf)
