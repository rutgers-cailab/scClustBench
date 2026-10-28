# ============================================================
# pipeline/12_decision_rule.R
#
# Applies the two-stage decision rule in R/decision_rule.R to tables written by earlier steps;
# it does no clustering. Stage A (sccb_select) screens and ranks the Phase 1 and rod benchmark
# grids, with a sensitivity check over the Jaccard and minimum-size gates. Stage B
# (sccb_verdict) gives ACCEPT or REJECT and k* for every R(k) curve (retina populations per
# embedding family, Phase 1, synthetic, external and null-calibration datasets), with tau raised
# to the family's null 95th percentile when a calibration table exists (sccb_tau). It then
# combines families per population (sccb_consensus), reruns the rule under variants, and picks
# the Stage A partition at the Stage B k*.
#
# Paper:    Stage A and Stage B results and their supplementary tables
# Inputs:   output_full_dataset/tables/poc_scores_matrix{,_4metric}.csv   (Phase 1 grid, steps 03, 12d)
#           output/tables/03_scores_matrix{,_4metric}.csv                 (rod grid, step 04)
#           output/tables/13_embedding_generality_curves.csv              (step 13c; if absent, the
#             scDHA curves 07_bc_stability_by_k_extended.csv and 07_rods_stability_by_k.csv)
#           output/tables/07_all_stability_by_k.csv                       (Phase 1 curve, step 07_all)
#           output/tables/09_synthetic_validation.csv                     (step 09)
#           output/tables/15_external_curves.csv, 15_external_scores.csv  (step 15; optional)
#           output/tables/16_null_calibration_curves.csv, 16_null_calibration_tau.csv (step 16; optional)
# Outputs:  output/tables/12_selection_phase1.csv            ranked Phase 1 grid
#           output/tables/12_selection_rods.csv              ranked rod grid
#           output/tables/12_selection_sensitivity.csv       top pick under each Jaccard / size gate
#           output/tables/12_verdicts.csv                    Stage B verdict for every curve
#           output/tables/12_verdicts_consensus.csv          per population: do all families agree?
#           output/tables/12_verdicts_tau_sensitivity.csv    verdict and k* under rule variants
#           output/tables/12_stage16_performance.csv         sensitivity, specificity, k* on step 16 sims
#           output/tables/12_partition_at_kstar.csv          Stage A pick at the Stage B k*
# Run:      Rscript run_pipeline.R 12   (or: Rscript pipeline/12_decision_rule.R)
# ============================================================

source("config.R"); source("R/decision_rule.R")

# ---- Stage A: screen and rank the benchmark grids ----
# load_grid(four, silrogue): the metric grid in `four`, with the silhouette column taken from
# `silrogue`, joined on method and target k; returns one row per partition
load_grid <- function(four, silrogue) {
  a <- read.csv(four); b <- read.csv(silrogue)[, c("method", "target_k", "silhouette")]
  merge(a, b, by = c("method", "target_k"))
}
grids <- list(
  phase1 = load_grid("output_full_dataset/tables/poc_scores_matrix_4metric.csv",
                     "output_full_dataset/tables/poc_scores_matrix.csv"),
  rods   = load_grid(file.path(TABLE_DIR, "03_scores_matrix_4metric.csv"),
                     file.path(TABLE_DIR, "03_scores_matrix.csv")))

cols <- c("rank", "method", "target_k", "actual_k", "smallest", SCCB_METRICS, "score", "admissible")
sens <- list()
for (nm in names(grids)) {
  sel <- sccb_select(grids[[nm]])
  keep <- intersect(c(cols, "ari", "nmi"), names(sel))
  write.csv(sel[, keep], file.path(TABLE_DIR, sprintf("12_selection_%s.csv", nm)), row.names = FALSE)
  cat(sprintf("\n=== %s: top 8 admissible ===\n", nm)); print(head(sel[sel$admissible, keep], 8))
  # Sensitivity: top pick under other Jaccard and minimum-cluster-size thresholds
  for (jm in c(0.5, 0.6, 0.75)) for (mc in c(20, 100)) {
    t1 <- sccb_select(grids[[nm]], jacc_min = jm, min_cells = mc)[1, ]
    sens[[length(sens) + 1]] <- data.frame(dataset = nm, jacc_min = jm, min_cells = mc,
      method = t1$method, k = t1$target_k, score = round(t1$score, 3),
      ari = if ("ari" %in% names(t1)) t1$ari else NA)
  }
}
sens <- do.call(rbind, sens); print(sens)
write.csv(sens, file.path(TABLE_DIR, "12_selection_sensitivity.csv"), row.names = FALSE)

# ---- Stage B: collect every R(k) curve ----
# tf(f): path of table f; rd(f): read table f, or NULL if it does not exist
tf <- function(f) file.path(TABLE_DIR, f)
rd <- function(f) if (file.exists(tf(f))) read.csv(tf(f)) else NULL
calib <- rd("16_null_calibration_tau.csv")
cs <- list()                                   # one entry per curve: list(case, pop, emb, d)
# add(case, pop, emb, d): append one curve to cs
add <- function(case, pop, emb, d) cs[[length(cs) + 1]] <<- list(case = case, pop = pop, emb = emb, d = d)
# Retina populations: every embedding family from step 13c, else the scDHA curves from 07 / 07b
g13 <- rd("13_embedding_generality_curves.csv")
if (!is.null(g13)) for (d in split(g13, list(g13$celltype, g13$embedding_family), drop = TRUE))
  add("retina", d$celltype[1], d$embedding_family[1], d) else {
  add("retina", "BC", "scdha", rd("07_bc_stability_by_k_extended.csv"))
  add("retina", "Rods", "scdha", rd("07_rods_stability_by_k.csv")) }
if (!is.null(p1 <- rd("07_all_stability_by_k.csv"))) add("retina", "Phase1_all", "scdha", p1)
syn <- rd("09_synthetic_validation.csv")
for (d in split(syn, syn$scenario)) add("synthetic09", d$scenario[1], "scdha", d)
ext <- rd("15_external_curves.csv")
if (!is.null(ext)) for (d in split(ext, list(ext$dataset, ext$embedding), drop = TRUE))
  add("external", d$dataset[1], d$embedding[1], d)
n16 <- rd("16_null_calibration_curves.csv")
if (!is.null(n16)) for (d in split(n16, list(n16$dataset, n16$embedding), drop = TRUE))
  add("stage16", d$dataset[1], d$embedding[1], d)

# ---- Stage B: verdict and k* for every curve ----
# meta(x): identifying columns for curve x (true_k only for simulated data)
meta <- function(x) data.frame(case = x$case, population = x$pop, embedding = x$emb,
  true_k = if ("true_k" %in% names(x$d)) x$d$true_k[1] else NA)
# run_rule(x, tau, min_cells, kpick): sccb_verdict() row for curve x, with k* switched to the smallest
# qualifying k when kpick = "smallest", plus the ARI against the truth at k* when the curve has it
run_rule <- function(x, tau = sccb_tau(x$emb, calib), min_cells = 20, kpick = "largest") {
  v <- sccb_verdict(x$d, tau = tau, min_cells = min_cells)
  if (kpick == "smallest") v$k_star <- v$k_star_smallest
  v$ari_at_k_star <- if ("ari_vs_truth" %in% names(x$d)) x$d$ari_vs_truth[match(v$k_star, x$d$k)] else NA
  v
}
ver <- do.call(rbind, lapply(cs, function(x) cbind(meta(x), run_rule(x))))
print(ver[ver$case != "stage16", ])
write.csv(ver, tf("12_verdicts.csv"), row.names = FALSE)

# ---- Consensus across embedding families, per real-data population ----
nonsim <- ver[ver$case %in% c("retina", "external"), ]
con <- do.call(rbind, lapply(split(nonsim, nonsim$population), function(d) data.frame(
  population = d$population[1], families = paste(d$embedding, collapse = "/"),
  verdicts = paste(d$verdict, collapse = "/"), k_star = paste(d$k_star, collapse = "/"),
  consensus = sccb_consensus(d$verdict))))
print(con); write.csv(con, tf("12_verdicts_consensus.csv"), row.names = FALSE)

# ---- Rule variants: fixed tau grid, the calibrated rule ("final"), no size gate, smallest k* ----
variants <- list(tau0.75 = list(tau = 0.75), tau0.80_fixed = list(tau = 0.8), tau0.85 = list(tau = 0.85),
                 final = list(), no_size_gate = list(min_cells = 0), k_smallest = list(kpick = "smallest"))
sens <- do.call(rbind, lapply(names(variants), function(nv) do.call(rbind, lapply(cs, function(x) {
  v <- do.call(run_rule, c(list(x), variants[[nv]]))
  cbind(meta(x), variant = nv, verdict = v$verdict, k_star = v$k_star, tau = v$tau) }))))
write.csv(sens, tf("12_verdicts_tau_sensitivity.csv"), row.names = FALSE)

# ---- Step 16 simulations: operating characteristics under each variant (plus the raw calibrated tau) ----
if (!is.null(n16)) {
  s16 <- sens[sens$case == "stage16", ]
  if (!is.null(calib)) s16 <- rbind(s16, do.call(rbind, lapply(cs[vapply(cs, `[[`, "", "case") == "stage16"], function(x) {
    v <- run_rule(x, tau = calib$tau_calibrated[match(x$emb, calib$embedding)])
    cbind(meta(x), variant = "tau_calibrated_raw", verdict = v$verdict, k_star = v$k_star, tau = v$tau) })))
  # type ("structured" or "null") per simulated dataset
  n16t <- unique(n16[, c("dataset", "type")]); s16$type <- n16t$type[match(s16$population, n16t$dataset)]
  # Sensitivity: structured datasets accepted; specificity: null datasets rejected;
  # k_star_correct: accepted structured datasets whose k* equals the true k
  perf <- do.call(rbind, lapply(split(s16, list(s16$embedding, s16$variant), drop = TRUE), function(d) {
    st <- d[d$type == "structured", ]; acc <- st$verdict == "ACCEPT"
    data.frame(embedding = d$embedding[1], variant = d$variant[1], tau = round(max(d$tau), 3),
      sensitivity = sprintf("%d/%d", sum(acc), nrow(st)),
      specificity = sprintf("%d/%d", sum(d$verdict[d$type == "null"] == "REJECT"), sum(d$type == "null")),
      k_star_correct = sprintf("%d/%d", sum(st$k_star[acc] == st$true_k[acc]), sum(acc))) }))
  print(perf); write.csv(perf, tf("12_stage16_performance.csv"), row.names = FALSE)
}

# ---- Partition at k*: Stage A ranking restricted to the Stage B k* ----
# pick(grid, k, label): top admissible Stage A row of grid at k, or NULL if k is NA or none is admissible
pick <- function(grid, k, label) {
  if (is.na(k)) return(NULL)
  s <- sccb_select(grid, k = k)
  s <- s[s$admissible, intersect(c("method", "target_k", "jacc_min", "smallest", "score", "ari"), names(s))]
  if (nrow(s)) cbind(population = label, k_star = k, s[1, ]) else NULL
}
at_k <- list(pick(grids$phase1, ver$k_star[ver$population == "Phase1_all"][1], "Phase1_all"))
# External datasets: one pick per (dataset, embedding family) verdict
sc15 <- rd("15_external_scores.csv")
if (!is.null(sc15)) for (d in split(ver[ver$case == "external", ], seq_len(sum(ver$case == "external"))))
  at_k[[length(at_k) + 1]] <- pick(sc15[sc15$dataset == d$population, ], d$k_star, paste(d$population, d$embedding))
at_k <- do.call(rbind, at_k); print(at_k)
write.csv(at_k, tf("12_partition_at_kstar.csv"), row.names = FALSE)
