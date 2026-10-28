# ============================================================
# run_pipeline.R
#
# Runs the scClustBench retina analysis end to end on one machine, or any
# subset of its steps. Each step is a script in pipeline/ run in its own
# Rscript process (so memory is released between steps), with the
# environment variables listed in STEPS below.
#
# Usage (from the repository root):
#   Rscript run_pipeline.R                  # every step, in the order below
#   Rscript run_pipeline.R 01 02 06         # only these steps
#   Rscript run_pipeline.R --from 07_bc     # this step and everything after it
#   Rscript run_pipeline.R --list           # print the steps and exit
#
# Steps run in dependency order; a step only needs the outputs of steps
# listed above it. The controls that need no retina data (15, 16, 09, 11)
# run right after the annotation, because steps 13c, 17 and 12 read their
# tables. The full run takes several days on a 32-core
# workstation; docs/HPC.md describes the Slurm version
# (hpc/slurm/submit_all.sh), which runs independent steps in parallel.
# ============================================================

# id, script in pipeline/, environment variables for that run
STEPS <- list(
  # ---- Preprocessing and reference annotation ----
  list("01",        "01_qc_harmony_integration.R"),
  list("02",        "02_cell_type_annotation.R"),
  list("02b",       "02b_annotation_resolution_sensitivity.R"),
  list("02c",       "02c_harmony_annotation_evidence.R"),
  # ---- Controls that need no retina data: external datasets, null calibration, synthetic data ----
  list("15",        "15_external_datasets.R"),
  list("16",        "16_null_calibration.R"),
  list("09",        "09_synthetic_validation.R"),
  list("11",        "11_synthetic_umaps.R"),
  # ---- Phase 1: whole retina ----
  list("03",        "03_phase1_full_dataset_benchmark.R"),
  list("12c_p1",    "12c_recompute_jaccard.R",          c(DATASET = "phase1")),
  list("12d",       "12d_rescore_phase1_ari.R"),
  # ---- Phase 2: rods (Stage A screen and the illustrative partition) ----
  list("04",        "04_phase2_rod_benchmark.R"),
  list("12c_rods",  "12c_recompute_jaccard.R",          c(DATASET = "rods")),
  list("04b",       "04b_chosen_partition_stability.R", c(POP = "Rods")),
  list("05",        "05_chosen_partition_analysis.R"),
  list("12b",       "12b_marker_technical_audit.R"),
  # ---- Stage B: R(k) across embedding families ----
  list("07_bc",     "07_positive_control_bc.R",         c(CELLTYPE = "BC")),
  list("07_rods",   "07_positive_control_bc.R",         c(CELLTYPE = "Rods")),
  list("07_all",    "07_positive_control_bc.R",         c(CELLTYPE = "ALL")),
  list("07b",       "07b_bc_extended_ksweep.R"),
  list("04b_bc",    "04b_chosen_partition_stability.R", c(POP = "BC")),
  list("04c",       "04c_seed_robustness.R"),
  list("13_bc",     "13_embedding_generality.R",        c(POP = "BC")),
  list("13_rods",   "13_embedding_generality.R",        c(POP = "Rods")),
  list("13s_bc",    "13s_scvi_family.R",                c(POP = "BC")),
  list("13s_rods",  "13s_scvi_family.R",                c(POP = "Rods")),
  list("13c",       "13c_combine.R"),
  list("13t",       "13t_scvi_rods_inspect.R"),
  list("17c",       "17c_scvi_rods_technical.R"),
  list("17b",       "17b_scvi_rod_split_audit.R"),
  list("17",        "17_stageB_gates.R"),
  list("12",        "12_decision_rule.R"),
  # ---- Sensitivity analyses ----
  list("04e",       "04e_knn_robustness.R"),
  list("14",        "14_preprocessing_sensitivity.R",   c(PART = "ab")),
  list("14_tables", "14_preprocessing_sensitivity.R",   c(PART = "c")),
  # ---- Paper figures (written to output/paper/) ----
  list("fig1",      "00_pipeline_schematic.py"),
  list("06",        "06_publication_figures.R"),        # Figures 2, 3, 5; Supplementary Figures S1, S2 and S5
  list("08b",       "08b_bc_positive_control_figure.R"),# Figure 4
  list("08",        "08_discovery_comparison_figure.R"),# Figure 6
  list("10",        "10_synthetic_figure.R"),           # Figure 7
  list("04d",       "04d_stability_figures.R"),         # Supplementary Figure S6
  list("13f",       "13f_supp_figures.R")               # Supplementary Figures S3, S4
)
ids <- vapply(STEPS, `[[`, "", 1)

args <- commandArgs(trailingOnly = TRUE)
if ("--list" %in% args) {
  for (s in STEPS) cat(sprintf("%-10s %-40s %s\n", s[[1]], s[[2]],
    if (length(s) > 2) paste(names(s[[3]]), s[[3]], sep = "=", collapse = " ") else ""))
  quit(save = "no")
}
if (length(args) && args[1] == "--from" && (length(args) < 2 || !args[2] %in% ids))
  stop("--from needs a known step id. See Rscript run_pipeline.R --list")
run <- if (!length(args)) ids else if (args[1] == "--from") ids[match(args[2], ids):length(ids)] else args
bad <- setdiff(run, ids)
if (length(bad)) stop("unknown step(s): ", paste(bad, collapse = ", "), ". See Rscript run_pipeline.R --list")

T0 <- Sys.time()
for (s in STEPS[ids %in% run]) {
  script <- file.path("pipeline", s[[2]])
  env <- if (length(s) > 2) paste0(names(s[[3]]), "=", s[[3]]) else character(0)
  message(sprintf("\n### [%s] step %s: %s %s", format(Sys.time(), "%H:%M:%S"), s[[1]], script, paste(env, collapse = " ")))
  cmd <- if (grepl("\\.py$", script)) "python3" else file.path(R.home("bin"), "Rscript")
  status <- system2(cmd, script, env = env)
  if (status != 0) stop(sprintf("step %s failed (exit status %d)", s[[1]], status))
}
message(sprintf("\nDone: %d step(s) in %s", sum(ids %in% run), format(Sys.time() - T0, digits = 3)))
