# ============================================================
# pipeline/02b_annotation_resolution_sensitivity.R
#
# Tests how much the cell-level annotation depends on the clustering resolution used to define
# the annotated clusters (reference: 0.6). It reuses the cluster columns already stored by step 01
# (SCT_snn_res.{0.2,0.4,0.6,0.8,1}) and does not re-cluster. For each resolution r, clusters are
# collapsed to cell types in two ways:
#   auto     : the automated rule of step 02 (panel with the highest cluster-mean module score,
#              "Unassigned" below CELLTYPE_MIN_SCORE)
#   majority : each cluster takes the majority reference label (the hand annotation at res 0.6).
#              This is the best any cluster-level labelling at resolution r can do, so its agreement
#              with the reference measures how well the clusters at r respect the reference classes.
# Agreement is reported against (i) the same collapse at r = 0.6 and (ii) the reference
# `cell_type` column, as percent of cells and ARI.
#
# Paper:    Supplementary analysis of annotation robustness (Methods, reference annotation)
# Inputs:   output/rds/02_annotated.rds
# Outputs:  output/tables/02b_annotation_resolution_sensitivity.csv
#           output/tables/02b_annotation_resolution_crosstab_<r>.csv  (majority vs reference, per r)
# Run:      Rscript run_pipeline.R 02b
# ============================================================
source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(aricode) })

# ---- Inputs: reference labels, panel module scores (from step 02), stored cluster columns ----
obj  <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
md   <- obj@meta.data
canon <- as.character(md$cell_type)
score_cols <- grep("^score_", colnames(md), value = TRUE)
panel_to_ct <- setNames(names(CELLTYPE_PANELS),
                        paste0("score_", gsub("[^A-Za-z0-9]", "_", names(CELLTYPE_PANELS))))
res_cols <- grep("^SCT_snn_res", colnames(md), value = TRUE)

# Takes a cluster vector cl (one per cell); returns per-cell labels from the automated module-score rule.
collapse_auto <- function(cl) {
  cm <- aggregate(md[, score_cols], by = list(cl = cl), FUN = mean)
  top <- apply(cm[, score_cols], 1, which.max); top_s <- apply(cm[, score_cols], 1, max)
  lab <- ifelse(top_s >= CELLTYPE_MIN_SCORE, panel_to_ct[score_cols[top]], "Unassigned")
  unname(as.character(unlist(setNames(lab, cm$cl)[as.character(cl)])))
}
# Takes a cluster vector cl; returns per-cell labels where each cluster takes its majority reference label.
collapse_major <- function(cl) {
  maj <- tapply(canon, cl, function(v) names(which.max(table(v))))
  unname(as.character(unlist(maj[as.character(cl)])))
}

# Collapses at the reference resolution, used as the baseline for the "vs 0.6" columns
ref_auto  <- collapse_auto(md$SCT_snn_res.0.6)
ref_major <- collapse_major(md$SCT_snn_res.0.6)
# Takes two label vectors; returns the percent of cells with identical labels.
pct <- function(a, b) round(100 * mean(a == b), 2)

# ---- One summary row per resolution, plus a majority-vs-reference crosstab ----
rows <- lapply(res_cols, function(rc) {
  cl <- md[[rc]]; r <- sub("SCT_snn_res\\.", "", rc)
  a <- collapse_auto(cl); m <- collapse_major(cl)
  ct <- as.data.frame.matrix(table(canonical = canon, majority = m))
  write.csv(ct, file.path(TABLE_DIR, sprintf("02b_annotation_resolution_crosstab_%s.csv", r)))
  data.frame(resolution = as.numeric(r), n_clusters = length(unique(cl)),
             auto_pct_vs_auto0.6 = pct(a, ref_auto),  auto_ari_vs_auto0.6 = round(ARI(a, ref_auto), 4),
             auto_pct_vs_canonical = pct(a, canon),   auto_ari_vs_canonical = round(ARI(a, canon), 4),
             major_pct_vs_major0.6 = pct(m, ref_major), major_ari_vs_major0.6 = round(ARI(m, ref_major), 4),
             major_pct_vs_canonical = pct(m, canon),  major_ari_vs_canonical = round(ARI(m, canon), 4),
             n_types_recovered_major = length(unique(m)))
})
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "02b_annotation_resolution_sensitivity.csv"), row.names = FALSE)
print(out)
message("canonical types: ", paste(sort(unique(canon)), collapse = ", "))
