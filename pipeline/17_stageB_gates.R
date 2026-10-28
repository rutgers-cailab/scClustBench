# ============================================================
# pipeline/17_stageB_gates.R
#
# Audits the Stage B validity gates and tau on the saved R(k) curves (tables
# only; no data or embeddings are loaded). Each curve is passed through
# sccb_verdict() after one of three explicit k filters:
#   as_is      no extra filter (sccb_verdict() applies its own gates)
#   consensus  keep k only if every embedding of the family reached k
#              (n_embeddings_at_k equals the number of embeddings), so R(k) is
#              never a mean over partitions with different k
#   size_gate  consensus plus the minimum cluster size gate
#              (smallest seed-1 cluster >= MIN_CELLS)
# sccb_verdict() in R/decision_rule.R applies both gates itself, so the three
# rows should agree; the table documents that the gates are in force. Each
# curve is judged at tau = 0.8 and at tau_cal, the 95th percentile of the null
# max R(k) under the same filter (step 16 nulls; scDHA and PCA gene resampling only).
#
# Paper:    Stage B validity gates and calibrated tau (supplementary tables)
# Inputs:   output/tables/13_embedding_generality_curves.csv, output/tables/13s_scvi_curves.csv,
#           output/tables/15_external_curves.csv, output/tables/16_null_calibration_curves.csv
# Outputs:  output/tables/17_stageB_gated_verdicts.csv
#           output/tables/17_stageB_gated_calibration.csv
# Run:      Rscript run_pipeline.R 17   or   Rscript pipeline/17_stageB_gates.R
# Runtime:  under 10 s
# ============================================================
source("config.R"); source("R/decision_rule.R")
MIN_CELLS <- 20
# rd(f): read table f from TABLE_DIR
rd <- function(f) read.csv(file.path(TABLE_DIR, f))

# gate(d, g): R(k) curve d -> the rows (k values) that pass filter g. sccb_verdict() applies a gate
# whenever its column is present, so the columns of the gates a filter leaves out are removed:
#   as_is      no gates
#   consensus  every embedding reaches k
#   size_gate  every embedding reaches k, and no reference cluster under MIN_CELLS cells
gate <- function(d, g) {
  n <- max(d$n_embeddings_at_k)
  if (g != "as_is") d <- d[d$n_embeddings_at_k == n, ]
  if (g == "size_gate") d <- d[vapply(d$seed1_sizes, .min_size, 1) >= MIN_CELLS, ]
  if (g == "as_is") d$n_embeddings_at_k <- NULL
  if (g != "size_gate") d$seed1_sizes <- NULL
  d
}
# verdict(d, tau, kmax): sccb_verdict() row for the kept k, or a NO_ADMISSIBLE_K row (same columns,
# all NA) if none are left. k_star_at_edge is re-referenced to the original sweep (kmax), since
# sccb_verdict() sees only the kept k
verdict <- function(d, tau, kmax = max(d$k)) {
  if (nrow(d)) return(transform(sccb_verdict(d, tau), k_star_at_edge = k_star %in% kmax))
  v <- suppressWarnings(sccb_verdict(data.frame(k = 1, actual_k = 1, mean_cross_seed_ari = NA_real_), tau))
  v[] <- NA; v$verdict <- "NO_ADMISSIBLE_K"; v$tau <- tau; v
}
GATES <- c("as_is", "consensus", "size_gate")

# ---- tau calibrated per embedding family and per filter on the step 16 nulls ----
nc <- rd("16_null_calibration_curves.csv")
nk <- split(nc, list(nc$dataset, nc$embedding), drop = TRUE)
# max R(k) of every simulated dataset under every filter (tau = Inf: nothing is accepted)
cal <- do.call(rbind, lapply(GATES, function(g) do.call(rbind, lapply(nk, function(d) {
  v <- verdict(gate(d, g), Inf)
  data.frame(gate = g, dataset = d$dataset[1], type = d$type[1], true_k = d$true_k[1],
             embedding = d$embedding[1], max_R = v$max_R)
}))))
tau_cal <- aggregate(max_R ~ gate + embedding, cal[cal$type == "null", ],
                     function(x) quantile(x, 0.95, names = FALSE))
names(tau_cal)[3] <- "tau_cal"
# Sensitivity / specificity / correct k* on the simulations at tau = 0.8 and at tau_cal.
# Here specificity counts any non-ACCEPT verdict on a null (including NO_ADMISSIBLE_K) as correct
perf <- do.call(rbind, lapply(seq_len(nrow(tau_cal)), function(i) {
  g <- tau_cal$gate[i]; e <- tau_cal$embedding[i]
  do.call(rbind, lapply(c(0.8, tau_cal$tau_cal[i]), function(t) {
    v <- do.call(rbind, lapply(nk[vapply(nk, function(d) d$embedding[1] == e, TRUE)], function(d)
      cbind(d[1, c("type", "true_k")], verdict(gate(d, g), t, max(d$k)))))
    s <- v[v$type == "structured", ]; a <- s$verdict == "ACCEPT"
    data.frame(gate = g, embedding = e, tau = round(t, 3), sensitivity = mean(a),
               specificity = mean(v$verdict[v$type == "null"] != "ACCEPT"),
               k_star_correct = mean(s$k_star[a] == s$true_k[a]),
               k_star_smallest_correct = mean(s$k_star_smallest[a] == s$true_k[a]))
  }))
}))
write.csv(perf, file.path(TABLE_DIR, "17_stageB_gated_calibration.csv"), row.names = FALSE); print(perf)

# ---- Real-data curves: retina (stage 13) and external datasets (stage 15) ----
g13 <- rd("13_embedding_generality_curves.csv")
g13 <- rbind(g13[g13$embedding_family != "scvi", ], rd("13s_scvi_curves.csv"))   # 13s is the scVI source
ext <- rd("15_external_curves.csv")
# Bring both sources to one layout: source, case, embedding, then the R(k) curve columns
cur <- rbind(data.frame(source = "13", case = g13$celltype, embedding = g13$embedding_family, g13[, -(1:2)]),
             data.frame(source = "15", case = ext$dataset, embedding = ext$embedding,
                        ext[, setdiff(names(ext), c("dataset", "n_cells", "n_true", "embedding", "ari_vs_truth"))]))
ck <- split(cur, list(cur$case, cur$embedding), drop = TRUE)
# One row per case x embedding x filter x tau; tau_cal exists only for families calibrated in step 16
out <- do.call(rbind, lapply(GATES, function(g) do.call(rbind, lapply(ck, function(d) {
  tc <- tau_cal$tau_cal[tau_cal$gate == g & tau_cal$embedding == d$embedding[1]]
  do.call(rbind, lapply(c(0.8, if (length(tc)) tc), function(t)
    cbind(source = d$source[1], case = d$case[1], embedding = d$embedding[1], gate = g,
          tau_type = if (t == 0.8) "fixed" else "calibrated", n_k_kept = nrow(gate(d, g)), verdict(gate(d, g), t, max(d$k)))))
}))))
out$tau <- round(out$tau, 3)
out <- out[order(out$source, out$case, out$embedding, out$tau_type, match(out$gate, GATES)), ]
write.csv(out, file.path(TABLE_DIR, "17_stageB_gated_verdicts.csv"), row.names = FALSE)
print(out[, c("case", "embedding", "gate", "tau_type", "tau", "n_k_kept", "verdict", "k_star", "k_star_smallest", "max_R")], row.names = FALSE)
