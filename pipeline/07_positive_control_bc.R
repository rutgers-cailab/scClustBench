# ============================================================
# pipeline/07_positive_control_bc.R
#
# Stage B reproducibility curve R(k) for one cell population, using the scDHA embedding family.
# Five scDHA embeddings are trained with different seeds and each is Harmony-corrected on sample.
# For every k in the sweep, each embedding is clustered with Louvain at that k, and R(k) is the
# mean pairwise ARI between the five partitions. Quality metrics of the seed-1 partition
# (silhouette, ROGUE, size balance) are recorded alongside. Run for bipolar cells (positive
# control), rods (Phase 2) and the whole retina (Phase 1).
#
# Paper:    Stage B curves behind Figures 4 and 6 and the verdicts of pipeline/12_decision_rule.R
# Inputs:   output/rds/02_annotated.rds (BC, ALL); output/rds/02_rods.rds (Rods), via load_population()
# Outputs:  output/tables/07_<tag>_stability_by_k.csv   (one row per k: R(k), min ARI, seed-1 metrics)
#           output/tables/07_<tag>_seed_sizes.csv       (seed-1 cluster sizes per k)
#           output/rds/07_<tag>_positive_control.rds    (curve, the five embeddings, the Seurat object)
#           output/rds/07_<tag>_seed<s>.rds             (one embedding, written when SEED_ONLY is set)
#           <tag> is CELLTYPE in lower case with non-letters removed: bc, rods, all
# Run:      CELLTYPE=BC Rscript pipeline/07_positive_control_bc.R   (or: Rscript run_pipeline.R 07_bc)
#           CELLTYPE=Rods ... (07_rods) and CELLTYPE=ALL ... (07_all) likewise
# Options:  CELLTYPE   "BC" (default), "Rods" (cone-filtered rods), "ALL" (whole dataset), or any
#                      cell_type label in 02_annotated.rds
#           K_SWEEP    R expression for the k values, e.g. "2:10" (default 2:10; POC_K_RANGE for ALL)
#           SEED_ONLY  train only this scDHA seed, save it to 07_<tag>_seed<s>.rds and exit (used by
#                      hpc/slurm/submit_all.sh to train seeds as parallel jobs; a later run without
#                      SEED_ONLY loads the saved seeds instead of retraining them)
# Runtime:  CELLTYPE=ALL needs about 64 GB RAM because counts are densified for scDHA and ROGUE.
# ============================================================

source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix)
  library(igraph); library(aricode); library(cluster); library(ROGUE)
  library(tibble); library(dplyr)
})
SEED_LIST <- SCDHA_SEED_LIST          # scDHA training seeds (config.R); seed "1" is the reference
K_SWEEP   <- 2:10
# log_step(m): print message m with the elapsed time since start
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0,digits=3), m))

# ---- Load the target cell population ----
CELLTYPE <- Sys.getenv("CELLTYPE", "BC")
TAG <- tolower(gsub("[^A-Za-z]", "", CELLTYPE))
log_step(sprintf("Target cell population: %s", CELLTYPE))
bc <- load_population(CELLTYPE)   # R/utils.R: "Rods" (cone-filtered), "ALL" (whole dataset), or a cell_type
if (CELLTYPE == "ALL") K_SWEEP <- POC_K_RANGE   # Phase 1 curve uses the Phase 1 k range (config.R)
if (nzchar(Sys.getenv("K_SWEEP"))) K_SWEEP <- eval(parse(text = Sys.getenv("K_SWEEP")))
gc(verbose=FALSE)
log_step(sprintf("%s cells: %d", CELLTYPE, ncol(bc)))

# Standard Seurat preprocessing; the processed object is saved with the embeddings for later steps
bc <- NormalizeData(bc, verbose=FALSE) |> FindVariableFeatures(nfeatures=3000, verbose=FALSE) |>
      ScaleData(verbose=FALSE) |> RunPCA(npcs=N_PCS_CLUSTERING, verbose=FALSE)
har_meta <- data.frame(sample=bc$sample, row.names=colnames(bc))
# scDHA input: cells x genes, log2(counts + 1), dense
counts_log2 <- log2(t(as.matrix(GetAssayData(bc, assay="RNA", layer="counts")))+1)


# ---- Train the scDHA embeddings, each Harmony-corrected on sample ----
suppressPackageStartupMessages({ library(torch); library(scDHA) })
if (!torch_is_installed()) install_torch()
# SEED_ONLY=<s>: train that one seed, save it to 07_<tag>_seed<s>.rds and stop (seeds run as parallel jobs);
# without SEED_ONLY, seeds already saved that way are loaded instead of retrained.
# seed_file(s): path of the saved embedding for seed s
SEED_ONLY <- Sys.getenv("SEED_ONLY")
seed_file <- function(s) file.path(RDS_DIR, sprintf("07_%s_seed%s.rds", TAG, s))
seed_lat <- list()
for (s in SEED_LIST) {
  if (nzchar(SEED_ONLY) && as.character(s) != SEED_ONLY) next
  # Reuse a saved seed; the cell order must match the current object
  if (!nzchar(SEED_ONLY) && file.exists(seed_file(s))) {
    seed_lat[[as.character(s)]] <- readRDS(seed_file(s)); stopifnot(identical(rownames(seed_lat[[as.character(s)]]), rownames(counts_log2)))
    log_step(sprintf("scDHA seed %d loaded from %s", s, basename(seed_file(s)))); next }
  log_step(sprintf("scDHA seed %d", s)); set.seed(s)
  sc <- scDHA(data=counts_log2, k=NULL, method="scDHA", sparse=FALSE, n=SCDHA_N_GENES,
              ncores=N_CORES, gen_fil=TRUE, do.clus=TRUE, seed=s)
  lat <- as.matrix(sc$latent); rownames(lat)<-rownames(counts_log2)
  set.seed(s)
  seed_lat[[as.character(s)]] <- harmony::RunHarmony(lat, har_meta, vars_use="sample",
      theta=HARMONY_THETA, lambda=HARMONY_LAMBDA, max_iter=HARMONY_MAX_ITER, verbose=FALSE)
  if (nzchar(SEED_ONLY)) { saveRDS(seed_lat[[as.character(s)]], seed_file(s)); log_step("seed saved"); quit(save = "no") }
}
stopifnot(length(seed_lat) == length(SEED_LIST))

# ---- ROGUE input: filtered dense counts, used for per-k cluster purity ----
counts_dense <- as.matrix(GetAssayData(bc, assay="RNA", layer="counts"))
rogue_expr <- ROGUE::matr.filter(counts_dense, min.cells=10, min.genes=10); ROGUE::SE_fun(rogue_expr)
# gini(x): Gini coefficient of cluster sizes x (0 = equal sizes); balance = 1 - gini
gini <- function(x){x<-sort(x);n<-length(x); if(n<=1||sum(x)==0)return(0); sum((2*seq_len(n)-n-1)*x)/(n*sum(x))}

# ---- Sweep k: R(k) across the seed embeddings, plus seed-1 quality metrics ----
rows<-list(); size_rows<-list()
ref_ids <- Reduce(intersect, lapply(seed_lat, rownames))
# Louvain for every (seed, k): one cached resolution path per seed, seeds in parallel;
# (one Louvain resolution search per seed, shared across k)
source("R/scClustBench.R")
lab_all <- .par(seed_lat, function(L) { p <- sccb_louvain_path(L, SEED); lapply(K_SWEEP, function(k) setNames(p(k), rownames(L))) })
res_k <- .par(K_SWEEP, function(k) {   # k values in parallel (ROGUE per k dominates)
  labs <- lapply(lab_all, `[[`, match(k, K_SWEEP))
  # Pairwise ARI across seeds, on the cells shared by all embeddings; R(k) is their mean
  aris<-c()
  sn<-names(labs)
  for (i in 1:(length(sn)-1)) for (j in (i+1):length(sn))
    aris<-c(aris, ARI(labs[[i]][ref_ids], labs[[j]][ref_ids]))
  # Quality metrics on the seed-1 (reference) partition; NA if a metric fails
  l1 <- labs[["1"]]
  sil <- tryCatch(compute_silhouette(seed_lat[["1"]], l1, max_cells=MAX_SIL_CELLS, seed=SEED), error=function(e)NA)
  rg  <- tryCatch({ r<-ROGUE::rogue(rogue_expr, labels=as.character(l1[colnames(rogue_expr)]),
                    samples=rep(1,ncol(rogue_expr)), platform="UMI"); mean(unlist(r),na.rm=TRUE)}, error=function(e)NA)
  bal <- 1 - gini(as.numeric(table(l1)))
  # n_embeddings_at_k and seed1_sizes feed the Stage B validity gates in sccb_verdict()
  row <- data.frame(k=k, actual_k=length(unique(l1)),
      n_embeddings_at_k=sum(vapply(labs, function(l) length(unique(l)) == k, TRUE)),
      mean_cross_seed_ari=round(mean(aris),4), min_cross_seed_ari=round(min(aris),4),
      silhouette=round(sil,4), rogue=round(rg,4), balance=round(bal,4),
      seed1_sizes=paste(sort(as.integer(table(l1)),decreasing=TRUE),collapse="/"))
  srow <- data.frame(k=k, seed1_sizes=paste(sort(as.integer(table(l1)),decreasing=TRUE),collapse="/"))
  log_step(sprintf("  k=%d: mean cross-seed ARI=%.3f (sil=%.3f rogue=%.3f bal=%.3f)", k, mean(aris), sil, rg, bal))
  list(row, srow)
}, cores = min(.ncores(), 6))
# ---- Write the curve, the seed-1 cluster sizes and the embeddings ----
rows <- lapply(res_k, `[[`, 1); size_rows <- lapply(res_k, `[[`, 2)
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, sprintf("07_%s_stability_by_k.csv", TAG)), row.names=FALSE)
write.csv(do.call(rbind,size_rows), file.path(TABLE_DIR, sprintf("07_%s_seed_sizes.csv", TAG)), row.names=FALSE)
saveRDS(list(stability=out, seed_lat=seed_lat, bc=bc), file.path(RDS_DIR, sprintf("07_%s_positive_control.rds", TAG)))
log_step(sprintf("=== %s DISCOVERY-STABILITY ===", CELLTYPE)); print(out)
log_step("DONE")
