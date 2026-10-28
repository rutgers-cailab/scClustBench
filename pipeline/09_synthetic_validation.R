# ============================================================
# pipeline/09_synthetic_validation.R
#
# Tests the Stage B reproducibility criterion on Splatter simulations with known truth. Null
# scenarios have one homogeneous group; structured scenarios have three groups at several
# levels of differential expression, including one with a rare group. For each scenario five
# scDHA embeddings are trained and, at each k, the script records R(k) (mean pairwise ARI across
# the five Louvain partitions), the ARI of the seed-1 partition against the true groups, the
# achieved k and how many embeddings reached k (so sccb_verdict() can apply its validity gate).
#
# Paper:    Figure 7 (panels B and C); synthetic verdicts in pipeline/12_decision_rule.R
# Inputs:   none (data are simulated)
# Outputs:  output/tables/09_synthetic_validation.csv
#           output/rds/09_synthetic_validation.rds
# Run:      Rscript run_pipeline.R 09   (or: Rscript pipeline/09_synthetic_validation.R)
# ============================================================

source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({
  library(splatter); library(SingleCellExperiment); library(Seurat)
  library(aricode); library(cluster)
})
# Cells and genes per simulation, scDHA seeds, and the k values swept
NG <- 5000; NGENES <- 2000; SEEDS <- SCDHA_SEED_LIST; KS <- 2:6
# log_step(m): print message m with the elapsed time since start
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0,digits=3), m))

# ---- Simulation scenarios: group proportions, DE probability and Splatter seed ----
scenarios <- list(
  list(id="null_A",     type="null",       true_k=1, group.prob=1,               de.prob=0,    seed=101),
  list(id="null_B",     type="null",       true_k=1, group.prob=1,               de.prob=0,    seed=202),
  list(id="strong_3",   type="structured", true_k=3, group.prob=c(.34,.33,.33),  de.prob=0.30, seed=303),
  list(id="medium_3",   type="structured", true_k=3, group.prob=c(.34,.33,.33),  de.prob=0.10, seed=404),
  list(id="weak_3",     type="structured", true_k=3, group.prob=c(.34,.33,.33),  de.prob=0.03, seed=505),
  list(id="rare_3",     type="structured", true_k=3, group.prob=c(.80,.15,.05),  de.prob=0.15, seed=606)
)

# louvain_to_k(latent, target_k): Louvain labels (named by cell) whose cluster count is closest to
# target_k over a fixed 50-point resolution grid
louvain_to_k <- function(latent, target_k) {
  o <- CreateSeuratObject(matrix(0, 2, nrow(latent), dimnames = list(c("d1", "d2"), rownames(latent))))
  o[["lat"]] <- CreateDimReducObject(embeddings=latent, key="L_", assay="RNA")
  o <- FindNeighbors(o, reduction="lat", dims=1:ncol(latent), verbose=FALSE)
  bl<-NULL; bd<-Inf
  for (r in exp(seq(log(1e-4), log(10), length.out=50))) {
    set.seed(SEED); o<-FindClusters(o, resolution=r, algorithm=1, verbose=FALSE)
    lab<-as.integer(Idents(o)); nk<-length(unique(lab))
    if (abs(nk-target_k)<bd){bd<-abs(nk-target_k); bl<-lab}; if(bd==0)break
  }
  names(bl)<-rownames(latent); bl
}

suppressPackageStartupMessages({ library(torch); library(scDHA) })
if (!torch_is_installed()) install_torch()

# ---- Simulate each scenario and compute R(k) and ARI against the truth ----
rows <- list()
for (sc in scenarios) {
  log_step(sprintf("=== simulate %s (%s, true_k=%d, de=%.2f) ===", sc$id, sc$type, sc$true_k, sc$de.prob))
  params <- newSplatParams(nGenes=NGENES, batchCells=NG, seed=sc$seed)
  if (sc$type=="null") {
    sim <- splatSimulate(params, verbose=FALSE); truth <- rep("G1", ncol(sim))
  } else {
    sim <- splatSimulateGroups(params, group.prob=sc$group.prob, de.prob=sc$de.prob, verbose=FALSE)
    truth <- as.character(colData(sim)$Group)
  }
  cts <- as.matrix(counts(sim)); colnames(cts) <- paste0("c", seq_len(ncol(cts)))
  names(truth) <- colnames(cts)
  data_scdha <- log2(t(cts) + 1)   # scDHA input: cells x genes, log2(counts + 1)

  # One scDHA embedding per seed (no Harmony: simulations have a single batch)
  lat <- list()
  for (s in SEEDS) { set.seed(s)
    r <- scDHA(data=data_scdha, k=NULL, method="scDHA", sparse=FALSE, n=min(NGENES,3000),
               ncores=N_CORES, gen_fil=TRUE, do.clus=TRUE, seed=s)
    m <- as.matrix(r$latent); rownames(m) <- rownames(data_scdha); lat[[as.character(s)]] <- m
  }
  ids <- Reduce(intersect, lapply(lat, rownames))
  for (k in KS) {
    labs <- lapply(lat, function(L) louvain_to_k(L, k))
    # R(k): mean pairwise ARI across seeds; truth and silhouette use the seed-1 partition
    aris <- c(); sn <- names(labs)
    for (i in 1:(length(sn)-1)) for (j in (i+1):length(sn)) aris <- c(aris, ARI(labs[[i]][ids], labs[[j]][ids]))
    ari_truth <- ARI(labs[["1"]][ids], truth[ids])
    sil <- tryCatch(compute_silhouette(lat[["1"]][ids,,drop=FALSE], labs[["1"]][ids],
                    max_cells=MAX_SIL_CELLS, seed=SEED), error=function(e) NA_real_)
    rows[[length(rows)+1]] <- data.frame(scenario=sc$id, type=sc$type, true_k=sc$true_k, de_prob=sc$de.prob,
        k=k, actual_k=length(unique(labs[["1"]])), n_seeds_at_k=sum(sapply(labs, function(l) length(unique(l)))==k),
        mean_cross_seed_ari=round(mean(aris),4), ari_vs_truth=round(ari_truth,4), silhouette=round(sil,4),
        seed1_sizes=paste(sort(as.integer(table(labs[["1"]])), decreasing=TRUE), collapse="/"))
    log_step(sprintf("  k=%d: cross-seed ARI=%.3f | ARI vs truth=%.3f | sil=%.3f", k, mean(aris), ari_truth, sil))
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "09_synthetic_validation.csv"), row.names=FALSE)
saveRDS(out, file.path(RDS_DIR, "09_synthetic_validation.rds"))

# ---- Console summary: peak R(k) per scenario ----
# Verdicts (tau, validity gates, k*) are computed from the CSV by pipeline/12_decision_rule.R.
peak <- aggregate(mean_cross_seed_ari ~ scenario + type + true_k, data = out, FUN = max)
log_step("=== SUMMARY: peak cross-seed ARI per scenario ==="); print(peak)
log_step("DONE")
