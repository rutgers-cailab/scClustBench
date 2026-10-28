# ============================================================
# pipeline/15_external_datasets.R
#
# Applies the full scClustBench workflow (sccb_run: benchmark, Stage A, R(k),
# Stage B) to public datasets outside the mouse retina, chosen where possible
# so that the reference labels were obtained independently of RNA clustering.
# The labels are used only post hoc (ARI/NMI); the rule never sees them. The
# cases include positive controls with known groups and negative or
# within-population controls with no expected sub-structure.
#
#   Case                 Data (platform)                         Reference labels             Expected
#   stoeckius_lines      4 human cell lines, HEK293/K562/KG-1/   hashtag (HTO) demultiplexing ACCEPT, k* = 4
#                        THP-1 (10x, cell hashing)               (experimental, not RNA)
#   kotliarov_pbmc       human PBMC, 20 donors (10x CITE-seq)    surface-protein gates:       high ARI
#                                                                T / B / NK / myeloid
#   baron_pancreas       human pancreas, 4 donors (inDrop)       published cell types (14)    high ARI
#   shekhar_bc           mouse bipolar cells (Drop-seq),         published BC types           ACCEPT, k* near
#                        Shekhar et al. 2016                     (14 labels; BC8/9 merged)    published count
#   stoeckius_hek        one cell line (HEK293 singlets from     none; truth = cell-cycle     REJECT, or ACCEPT only
#                        stoeckius_lines): negative control      phase (post hoc audit)       k* = cell cycle
#   baron_beta           beta cells only (4 donors): a within-   none; truth = donor          rod analogue:
#                        type population, like the rods          (post hoc batch audit)       any split audited
#
#   Protein gates (kotliarov_pbmc): per marker, a two-component Gaussian mixture on
#   log1p(ADT counts) plus seeded jitter in [-0.1, 0.1]; a cell is called positive or
#   negative only when the posterior of that component is > 0.9, and cells that
#   satisfy no lineage gate (or more than one) are excluded. CD14 staining is
#   unimodal in this panel and is not used.
#   Cell calling (stoeckius_lines): barcodes with >= 500 RNA UMIs; HTOs summed per
#   cell line, CLR-normalised, Seurat::HTODemux (positive.quantile 0.99); singlets only.
#
# Paper:    External datasets (positive and negative controls outside the retina)
# Inputs:   Bioconductor scRNAseq datasets (downloaded through ExperimentHub on first use)
# Outputs:  output/tables/15_external/<case>/sccb_{scores,selection,curve,verdict,...}.csv, meta.rds
#           output/tables/15_external_labels_<case>.csv
#           output/tables/15_external_{scores,selection,curves,verdicts}.csv (all finished cases)
# Run:      Rscript run_pipeline.R 15   or   CASE=baron_beta Rscript pipeline/15_external_datasets.R
# Options:  CASE  comma list of case names to run (default: all, one after another);
#                 a name that matches no case (e.g. CASE=none) only rebuilds the combined tables
# Runtime:  about 1 to 3 h per case with scDHA and PCA gene resampling; shekhar_bc is the largest.
# ============================================================
if (!requireNamespace("scRNAseq", quietly = TRUE))
  stop('Needs Bioconductor scRNAseq: BiocManager::install("scRNAseq")')
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages({ library(SingleCellExperiment); library(mclust) })
# Stage B families: scDHA when installed, and always PCA gene resampling
REPRO <- c(if (requireNamespace("scDHA", quietly = TRUE)) "scdha", "pca_perturb")
OUT15 <- file.path(TABLE_DIR, "15_external")
N_MAX <- 12000                                   # subsample cap per case (runtime)

# ---- Helpers ----
# to_seurat(m, meta): count matrix (genes x cells) and per-cell metadata -> Seurat object
to_seurat <- function(m, meta) {
  rownames(m) <- make.unique(rownames(m)); colnames(m) <- rownames(meta)
  CreateSeuratObject(as(m, "dgCMatrix"), meta.data = meta)
}
# cap(o, n): random subsample of at most n cells (seeded), or o unchanged if it is smaller
cap <- function(o, n = N_MAX) { if (ncol(o) <= n) return(o); set.seed(SEED); subset(o, cells = sample(colnames(o), n)) }

# gate(x): ADT counts of one marker -> +1 (positive), -1 (negative) or NA (ambiguous) per cell
gate <- function(x) {
  # small fixed-seed jitter: low ADT counts are discrete and make a mixture component degenerate
  set.seed(SEED); y <- log1p(x) + runif(length(x), -0.1, 0.1)
  g <- Mclust(y, G = 2, modelNames = "V", verbose = FALSE)
  if (is.null(g)) stop("mixture fit failed")
  hi <- which.max(g$parameters$mean); lo <- 3 - hi
  ifelse(g$z[, hi] > 0.9, 1L, ifelse(g$z[, lo] > 0.9, -1L, NA_integer_))
}

# ---- Case builders ----
# Each returns list(obj = Seurat object with a `truth` column, batch = metadata column
# for Harmony or NULL, k = k range to sweep)
build <- list(
  stoeckius_lines = function() {
    h <- scRNAseq::StoeckiusHashingData(type = "mixed"); keep <- colSums(counts(h)) >= 500; h <- h[, keep]
    # Sum the three hashtags of each cell line (suffixes _A, _B, _C)
    H <- as.matrix(counts(altExp(h, "hto"))); L <- rowsum(H, sub("_[ABC]$", "", rownames(H)))
    # HTODemux needs a Seurat object; a 10-gene RNA assay is enough to carry the HTO assay
    d <- CreateSeuratObject(counts(h)[1:10, ], min.cells = 0, min.features = 0)
    d[["HTO"]] <- CreateAssayObject(counts = L); d <- NormalizeData(d, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
    d <- HTODemux(d, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
    sing <- d$HTO_classification.global == "Singlet"
    meta <- data.frame(truth = as.character(d$hash.ID[sing]), row.names = colnames(h)[sing])
    list(obj = to_seurat(counts(h)[, sing], meta), batch = NULL, k = 2:8)
  },
  kotliarov_pbmc = function() {
    s <- kotliarov(); keep <- !is.na(s$lineage)
    meta <- data.frame(truth = s$lineage[keep], batch = as.character(s$sce$batch[keep]), row.names = colnames(s$sce)[keep])
    list(obj = cap(to_seurat(counts(s$sce)[, keep], meta)), batch = "batch", k = 2:10)
  },
  baron_pancreas = function() {
    b <- scRNAseq::BaronPancreasData("human")
    meta <- data.frame(truth = b$label, donor = b$donor, row.names = colnames(b))
    list(obj = to_seurat(counts(b), meta), batch = "donor", k = 2:16)
  },
  shekhar_bc = function() {
    # Bipolar clusters only (names starting BC or RBC); the label is the first word of CLUSTER
    s <- scRNAseq::ShekharRetinaData(); bc <- grepl("^BC|^RBC", s$CLUSTER)
    meta <- data.frame(truth = sub(" .*", "", s$CLUSTER[bc]), row.names = colnames(s)[bc])
    list(obj = cap(to_seurat(counts(s)[, bc], meta), 15000), batch = NULL, k = 2:20)
  },
  # Negative and within-population controls: `truth` holds an audit covariate, not a cell type
  stoeckius_hek = function() {
    o <- build$stoeckius_lines()$obj; o <- NormalizeData(o[, o$truth == "HEK"], verbose = FALSE)
    o$truth <- tryCatch(CellCycleScoring(o, s.features = cc.genes.updated.2019$s.genes,
                                         g2m.features = cc.genes.updated.2019$g2m.genes)$Phase,
                        error = function(e) { message("cell-cycle scoring failed: ", conditionMessage(e)); "unknown" })
    list(obj = o, batch = NULL, k = 2:8)
  },
  baron_beta = function() {
    o <- build$baron_pancreas()$obj; o <- o[, o$truth == "beta"]; o$truth <- o$donor
    list(obj = o, batch = "donor", k = 2:8)
  })

# kotliarov(): CITE-seq PBMC data plus a protein-gated lineage per cell
# (T, B, NK, Myeloid, or NA when no single lineage gate is satisfied)
kotliarov <- function() {
  sce <- scRNAseq::KotliarovPBMCData(); A <- as.matrix(counts(altExp(sce, "ADT")))
  rownames(A) <- sub("_PROT$", "", trimws(rownames(A)))
  G <- sapply(c("CD3", "CD19", "CD20", "CD56", "CD16", "HLA-DR", "CD11c"), function(g) gate(A[g, ]))
  # p(g) / n(g): cell is confidently positive / negative for marker g
  p <- function(g) !is.na(G[, g]) & G[, g] == 1; n <- function(g) !is.na(G[, g]) & G[, g] == -1
  cand <- cbind(T = p("CD3") & n("CD19"),
                B = (p("CD19") | p("CD20")) & n("CD3"),
                NK = n("CD3") & n("CD19") & n("HLA-DR") & (p("CD56") | p("CD16")),
                Myeloid = n("CD3") & n("CD19") & p("HLA-DR") & p("CD11c"))
  lineage <- ifelse(rowSums(cand) == 1, colnames(cand)[max.col(cand)], NA)
  list(sce = sce, lineage = lineage)
}

# ---- Run each selected case ----
cases <- if (nzchar(Sys.getenv("CASE"))) intersect(strsplit(Sys.getenv("CASE"), ",")[[1]], names(build)) else names(build)   # CASE=none: combine only
for (nm in cases) {
  message("\n##### ", nm)
  b <- build[[nm]](); o <- b$obj
  lab <- as.data.frame(table(truth = o$truth)); write.csv(lab, file.path(TABLE_DIR, sprintf("15_external_labels_%s.csv", nm)), row.names = FALSE)
  message(ncol(o), " cells; reference labels: ", paste(sprintf("%s=%d", lab$truth, lab$Freq), collapse = ", "))
  r <- sccb_run(o, batch = b$batch, k_range = b$k, embeddings = "pca",
                methods = c("kmeans", "gmm", "louvain", "leiden"), repro_embeddings = REPRO,
                n_boot = BENCH_JACCARD_BOOTSTRAPS, truth = "truth", out_dir = file.path(OUT15, nm))
  # Add the post hoc ARI of the reference partition at k* to each verdict row
  v <- r$verdict
  v$ari_at_k_star <- mapply(function(e, k) { d <- r$curve[r$curve$embedding == e & r$curve$k == k, ]
                                             if (nrow(d)) d$ari_vs_truth[1] else NA }, v$embedding, v$k_star)
  write.csv(v, file.path(OUT15, nm, "sccb_verdict.csv"), row.names = FALSE)
  # meta.rds marks the case as finished and stores the counts used in the combined tables
  saveRDS(list(n_cells = ncol(o), n_true = nrow(lab)), file.path(OUT15, nm, "meta.rds"))
}

# ---- Combine whatever cases have finished ----
done <- intersect(list.dirs(OUT15, recursive = FALSE, full.names = FALSE), names(build))
# verdicts are always recomputed from the saved curves, so a change to sccb_verdict() needs no rerun
for (nm in done) {
  cv <- read.csv(file.path(OUT15, nm, "sccb_curve.csv"))
  v <- do.call(rbind, lapply(split(cv, cv$embedding), function(d) cbind(embedding = d$embedding[1], sccb_verdict(d, tau = 0.8))))
  v$ari_at_k_star <- mapply(function(e, k) { d <- cv[cv$embedding == e & cv$k == k, ]
                                             if (nrow(d)) d$ari_vs_truth[1] else NA }, v$embedding, v$k_star)
  write.csv(v, file.path(OUT15, nm, "sccb_verdict.csv"), row.names = FALSE)
}
done <- done[file.exists(file.path(OUT15, done, "meta.rds"))]
# rbind_fill(l): row-bind data frames whose columns differ, filling missing columns with NA
# (tables from different cases need not share every column)
rbind_fill <- function(l) { cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] })) }
# One combined table per result type; for selection only the top-ranked row of each case is kept
for (t in c("scores", "selection", "curve", "verdict")) {
  tabs <- lapply(done, function(nm) { f <- file.path(OUT15, nm, paste0("sccb_", t, ".csv")); if (!file.exists(f)) return(NULL)
    d <- read.csv(f); if (t == "selection") d <- d[1, , drop = FALSE]; m <- readRDS(file.path(OUT15, nm, "meta.rds"))
    cbind(dataset = nm, n_cells = m$n_cells, n_true = m$n_true, d) })
  tabs <- Filter(Negate(is.null), tabs)
  if (length(tabs)) write.csv(rbind_fill(tabs), file.path(TABLE_DIR, sprintf("15_external_%s.csv", sub("curve$", "curves", sub("verdict$", "verdicts", t)))), row.names = FALSE)
}
v <- file.path(TABLE_DIR, "15_external_verdicts.csv"); if (file.exists(v)) print(read.csv(v))
