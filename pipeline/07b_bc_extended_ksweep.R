# ============================================================
# pipeline/07b_bc_extended_ksweep.R
#
# Extends the bipolar cell R(k) curve of pipeline/07_positive_control_bc.R to a wider k range
# (default k = 2 to 18), which spans the number of transcriptomic bipolar types reported by
# Shekhar et al. (2016). It reuses the five scDHA embeddings saved by step 07 (no retraining)
# and computes the same quantities as step 07: R(k) as the mean pairwise ARI across embeddings,
# and seed-1 silhouette, ROGUE and size balance.
#
# Paper:    bipolar R(k) curve in Figures 4 and 6 (both prefer this table when it exists)
# Inputs:   output/rds/07_bc_positive_control.rds   (seed_lat and bc, from CELLTYPE=BC step 07)
# Outputs:  output/tables/07_bc_stability_by_k_extended.csv
#           output/tables/07_bc_seed_sizes_extended.csv
# Run:      Rscript run_pipeline.R 07b   (or: Rscript pipeline/07b_bc_extended_ksweep.R)
# Options:  K_EXT   k range as "from:to" (default "2:18")
# ============================================================

source("config.R"); source("R/utils.R")
suppressPackageStartupMessages({ library(Seurat); library(aricode); library(ROGUE); library(tibble); library(dplyr) })

# Parse K_EXT ("from:to") into an integer sequence
K_EXT <- as.integer(strsplit(Sys.getenv("K_EXT", "2:18"), ":")[[1]])
K_EXT <- seq(K_EXT[1], K_EXT[2])
# log_step(m): print message m with the elapsed time since start
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0, digits=3), m))

pc <- readRDS(file.path(RDS_DIR, "07_bc_positive_control.rds"))
seed_lat <- pc$seed_lat
bc <- pc$bc
log_step(sprintf("Loaded %d cached embeddings (seeds %s); %d cells",
                 length(seed_lat), paste(names(seed_lat), collapse=","), ncol(bc)))

# ---- ROGUE input: filtered dense counts (same as step 07) ----
counts_dense <- as.matrix(GetAssayData(bc, assay="RNA", layer="counts"))
rogue_expr <- ROGUE::matr.filter(counts_dense, min.cells=10, min.genes=10); ROGUE::SE_fun(rogue_expr)
# gini(x): Gini coefficient of cluster sizes x (0 = equal sizes); balance = 1 - gini
gini <- function(x){x<-sort(x);n<-length(x); if(n<=1||sum(x)==0)return(0); sum((2*seq_len(n)-n-1)*x)/(n*sum(x))}

# ---- Sweep k: R(k) across the seed embeddings, plus seed-1 quality metrics ----
rows<-list(); size_rows<-list()
ref_ids <- Reduce(intersect, lapply(seed_lat, rownames))
source("R/scClustBench.R")   # cached Louvain resolution path per seed, seeds in parallel
lab_all <- .par(seed_lat, function(L) { p <- sccb_louvain_path(L, SEED); lapply(K_EXT, function(k) setNames(p(k), rownames(L))) })
res_k <- .par(K_EXT, function(k) {   # k values in parallel (ROGUE per k dominates)
  labs <- lapply(lab_all, `[[`, match(k, K_EXT))
  # Pairwise ARI across seeds, on the cells shared by all embeddings; R(k) is their mean
  aris<-c(); sn<-names(labs)
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
  log_step(sprintf("  k=%d: mean cross-seed ARI=%.3f (min=%.3f sil=%.3f rogue=%.3f bal=%.3f)",
                   k, mean(aris), min(aris), sil, rg, bal))
  list(row, srow)
}, cores = min(.ncores(), 6))
# ---- Write the curve and the seed-1 cluster sizes ----
rows <- lapply(res_k, `[[`, 1); size_rows <- lapply(res_k, `[[`, 2)
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "07_bc_stability_by_k_extended.csv"), row.names=FALSE)
write.csv(do.call(rbind,size_rows), file.path(TABLE_DIR, "07_bc_seed_sizes_extended.csv"), row.names=FALSE)
log_step(sprintf("Done. Peak mean cross-seed ARI = %.3f at k=%d (range k=%d-%d).",
                 max(out$mean_cross_seed_ari), out$k[which.max(out$mean_cross_seed_ari)], min(K_EXT), max(K_EXT)))
