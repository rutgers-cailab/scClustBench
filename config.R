# ============================================================
# config.R
#
# Shared settings for every script in pipeline/: output folders, the six GEO
# samples and their locations, QC thresholds, Harmony integration and
# clustering parameters, marker panels, colour palettes and core counts.
# Scripts run from the repository root and start with source("config.R")
# followed by source("R/utils.R").
#
# To change an analysis setting, edit the value here; every step reads it
# from this file. To point the pipeline at the raw data, set SCCB_RAW_DATA
# (see RAW_DATA_BASE below). To add or remove a sample, edit SAMPLES.
#
# Options:  SCCB_RAW_DATA        folder holding the GEO sample folders
#                                (default: the parent of the repository root)
#           SLURM_CPUS_PER_TASK  cores to use (default: all detected cores minus one)
# ============================================================

# ---- Paths (all relative to the repository root) ----
PIPELINE_DIR   <- "."
RAW_DATA_BASE  <- Sys.getenv("SCCB_RAW_DATA", file.path(PIPELINE_DIR, ".."))   # folder holding the six GEO sample folders

RDS_DIR    <- file.path(PIPELINE_DIR, "output", "rds")
TABLE_DIR  <- file.path(PIPELINE_DIR, "output", "tables")
EXPORT_DIR <- file.path(PIPELINE_DIR, "output", "exports")
# Working figures written by individual steps; the paper figures go to PAPER_DIR
FIG_QC     <- file.path(PIPELINE_DIR, "output", "figures", "fig1_qc")
FIG_CT     <- file.path(PIPELINE_DIR, "output", "figures", "fig2_celltypes")
FIG_PUB    <- file.path(PIPELINE_DIR, "output", "figures", "publication")
PAPER_DIR  <- file.path(PIPELINE_DIR, "output", "paper")      # final paper figures (Figure_1-7, Supplementary_Figure_S1-S6) and source CSVs

# Create every output folder up front so scripts can write without checking
for (d in c(RDS_DIR, TABLE_DIR, EXPORT_DIR,
            FIG_QC, FIG_CT, FIG_PUB, PAPER_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Sample metadata ----
# One row per 10x sample folder under RAW_DATA_BASE: the six public GEO samples
# analysed in the paper.
SAMPLES <- tibble::tribble(
  ~sample, ~dir,                                                         ~time_label, ~animal,  ~condition,           ~source,
  "S1",    file.path(RAW_DATA_BASE, "P15_GSM8133387"),                   "P15",       "Mouse",  "P015 B6C3 WT1",      "YuChen Lab",
  "S2",    file.path(RAW_DATA_BASE, "P15_GSM8133389"),                   "P15",       "Mouse",  "P015 B6C3 WT2",      "YuChen Lab",
  "S3",    file.path(RAW_DATA_BASE, "P49_GSM7734028"),                   "P49",       "Mouse",  "P049 C57BL/6J",      "University Oxford",
  "S4",    file.path(RAW_DATA_BASE, "P75_GSM7720305"),                   "P75",       "Mouse",  "P075 C57BL/6J",      "Johns Hopkins University",
  "S5",    file.path(RAW_DATA_BASE, "P75_GSM7720306"),                   "P75",       "Mouse",  "P075 C57BL/6J",      "Johns Hopkins University",
  "S6",    file.path(RAW_DATA_BASE, "P91_GSM6205478_genes_no_features"), "P91",       "Mouse",  "P091 C57BL/6J",      "Peking Univ First Hospital"
)

# ---- QC thresholds (step 01) ----
# Cells are kept when counts, features, % mitochondrial and % haemoglobin reads
# fall inside these bounds; genes must be detected in at least QC_GENE_MIN_CELLS cells.
QC_MIN_COUNTS     <- 400
QC_MAX_COUNTS     <- 80000
QC_MIN_FEATURES   <- 400
QC_MAX_FEATURES   <- 8000
QC_MAX_PCT_MT     <- 18
QC_MAX_PCT_HB     <- 0.30
QC_GENE_MIN_CELLS <- 5

# ---- Integration / Harmony parameters (step 01) ----
# DEFAULT_RESOLUTION gives the 34 Harmony clusters that carry the reference annotation.
N_PCS_INTEGRATION   <- 100
N_DIMS_NEIGHBORS    <- 15
CLUSTER_RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0)
DEFAULT_RESOLUTION  <- 0.6
HARMONY_THETA       <- 2     # Harmony default
HARMONY_LAMBDA      <- 1     # Harmony default
HARMONY_MAX_ITER    <- 10    # Harmony default

# ---- Clustering benchmark settings (steps 03 and 04) ----
K_RANGE          <- 2:9       # Phase 2 (rod) k sweep; k = 1 is trivial and not benchmarked
N_PCS_CLUSTERING <- 50
SCDHA_N_GENES    <- 3000
KNN_K            <- 20
MAX_SIL_CELLS    <- 10000
JACCARD_BOOTSTRAPS <- 50      # default bootstrap re-clusterings for the dedicated stability
                              # analysis of the chosen partition (04b, shown by 04d)
BENCH_JACCARD_BOOTSTRAPS <- 3 # bootstraps inside the full benchmark grid (steps 03, 04, 14, 15);
                              # kept small because graph methods run a 60-step resolution
                              # search inside every bootstrap
SEED             <- 1

# ---- Stability and robustness ----
STABILITY_K_RANGE <- c(2, 3, 4)        # k values profiled for bootstrap stability (04b, 04e)
SCDHA_SEED_LIST   <- c(1, 2, 3, 7, 42) # embedding training seeds for Stage B R(k) (04c, 07, 13, ...)

# Embedding x clustering configurations benchmarked in steps 03 and 04
# (LL = Louvain on the SNN graph, Leiden uses the CPM objective)
METHODS <- c("PCA_KMeans", "PCA_GMM", "scDHA_KMeans",
             "PCA_LL",     "scDHA_LL",
             "PCA_Leiden", "scDHA_Leiden")

# Illustrative rod partition (Phase 2): the configuration and k examined in
# detail by steps 04b, 04c, 05 and 06. Its reproducibility across scDHA training
# seeds is measured in 04b/04c, and it is reported together with that result
# rather than as a set of rod subtypes.
CHOSEN_METHOD <- "scDHA_LL"
CHOSEN_K      <- 3

# ---- Phase 1 (positive control) benchmark settings ----
# Phase 1 (whole retina) and the bipolar control sweep a wider k range so the
# eight known cell classes fall inside it.
POC_K_RANGE <- 2:12

# ---- Cone contamination filter (rod population; see load_population() in R/utils.R) ----
# A cell is dropped from the rods if its cone score is in the top
# (1 - CONE_QUANTILE_CUT) and exceeds its rod score, or exceeds it by CONE_DELTA_CUT.
CONE_PANEL <- c("Opn1sw","Opn1mw","Pde6h","Pde6c","Gngt2","Gnat2",
                "Arr3","Pcp2","Gng13","Cngb3","Cnga3","Crx","Tulp1",
                "Slc24a2","Guca1a","Guca1b")
ROD_PANEL  <- c("Rho","Nrl","Nr2e3","Gnat1","Pde6a","Pde6b","Pde6g",
                "Gngt1","Rcvrn","Sag","Cnga1","Cngb1","Slc24a1","Reep6",
                "Rom1","Prph2","Aipl1")
CONE_QUANTILE_CUT <- 0.95
CONE_DELTA_CUT    <- 0.10

# ---- Cell-type marker panels (steps 01, 02, 02b, 02c) ----
# Each entry: cell-type name -> canonical marker genes. Step 02 scores every
# panel with AddModuleScore and assigns each cluster the panel with the highest
# mean score, provided it reaches CELLTYPE_MIN_SCORE (else "Unassigned"). This
# automated label is kept as cell_type_auto; the reference cell_type comes from
# the hand annotation in data/annotation/harmony_cluster_annotation.csv, whose
# class names must match the names below.
CELLTYPE_PANELS <- list(
  "Rods"        = c("Rho","Nrl","Nr2e3","Gnat1","Pde6a","Pde6b","Pde6g",
                    "Sag","Cnga1","Cngb1","Reep6"),
  "Cones"       = c("Opn1sw","Opn1mw","Pde6h","Pde6c","Arr3","Gngt2",
                    "Gnat2","Cngb3","Cnga3","Pcp2"),
  "BC"          = c("Vsx2","Otx2","Grik1","Scgn","Prkca","Grm6","Vsx1"),
  "Muller Glia" = c("Vim","Apoe","Glul","Rlbp1","Aqp4","Slc1a3"),
  "RGC/Neuron"  = c("Slc17a6","Pou4f1","Pou4f2","Rbpms","Nefl","Snhg11",
                    "Thy1"),
  "AC"          = c("Slc6a9","Pax6","Gad1","Gad2","Slc32a1","C1ql2"),
  "RPE"         = c("Rpe65","Mlana"),
  "Microglia"   = c("C1qa","Cx3cr1","Tmem119","Csf1r")
)
CELLTYPE_MIN_SCORE <- 0.05    # min module score above which assignment is trusted

# ---- Colour palettes (names must match cell-type, method, cluster and time-point labels) ----
CELLTYPE_COLORS <- c(
  "Rods"        = "#E41A1C", "Cones"       = "#377EB8",
  "BC"          = "#4DAF4A", "Muller Glia" = "#984EA3",
  "RGC/Neuron"  = "#FF7F00", "AC"          = "#A65628",
  "RPE"         = "#F781BF", "Microglia"   = "#999999",
  "Unassigned"  = "#CCCCCC"
)

# ---- Figure settings (used by save_fig() and pub_theme() in R/utils.R) ----
FIG_DPI       <- 600
PUB_BASE_SIZE <- 10

# ---- Parallel processing ----
# Cores for scDHA, BLAS and parallel loops (R/scClustBench.R reads SLURM_CPUS_PER_TASK directly)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L)))  # respect Slurm allocation
