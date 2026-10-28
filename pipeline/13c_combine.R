# ============================================================
# pipeline/13c_combine.R
#
# Merges the per-population R(k) curves of step 13 (scDHA, PCA gene resampling, PCA
# hyperparameters) with those of step 13s (scVI), applies the Stage B rule (sccb_verdict, with
# tau from sccb_tau: 0.8, raised to the family's null calibration from step 16 when that is
# higher) to each population and embedding family, and draws a summary plot of all curves.
#
# Paper:    Stage B embedding-family curves used by pipeline/12_decision_rule.R and Supplementary
#           Figures S5, S6 (pipeline/13f)
# Inputs:   output/tables/13_embedding_generality_curves_{bc,rods}.csv  (step 13 with POP=BC / Rods)
#           output/tables/13s_scvi_curves_{bc,rods}.csv                  (step 13s with POP=BC / Rods)
#           output/tables/16_null_calibration_tau.csv                    (step 16; optional)
# Outputs:  output/tables/13_embedding_generality_curves.csv
#           output/tables/13_embedding_generality_verdicts.csv
#           output/tables/13s_scvi_curves.csv
#           output/figures/publication/13_embedding_generality.{pdf,png}
# Run:      Rscript run_pipeline.R 13c   (or: Rscript pipeline/13c_combine.R)
# ============================================================
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages(library(ggplot2))
# rd(f): row-bind the bc and rods versions of table pattern f (sprintf with the population tag)
rd <- function(f) do.call(rbind, lapply(c("bc", "rods"), function(p) read.csv(file.path(TABLE_DIR, sprintf(f, p)))))
sc <- rd("13s_scvi_curves_%s.csv"); write.csv(sc, file.path(TABLE_DIR, "13s_scvi_curves.csv"), row.names = FALSE)
# scVI curves come from 13s; any scVI rows written by step 13 are replaced
cur <- rd("13_embedding_generality_curves_%s.csv"); cur <- rbind(cur[cur$embedding_family != "scvi", ], sc[, names(cur)])
# Stage B verdict per (population, embedding family)
# tau: 0.8, raised to the step 16 null calibration of the family when that is higher (as in step 12)
calib_f <- file.path(TABLE_DIR, "16_null_calibration_tau.csv"); calib <- if (file.exists(calib_f)) read.csv(calib_f) else NULL
ver <- do.call(rbind, lapply(split(cur, list(cur$celltype, cur$embedding_family), drop = TRUE), function(d)
  cbind(celltype = d$celltype[1], embedding_family = d$embedding_family[1], sccb_verdict(d, tau = sccb_tau(d$embedding_family[1], calib)))))
write.csv(cur, file.path(TABLE_DIR, "13_embedding_generality_curves.csv"), row.names = FALSE)
write.csv(ver, file.path(TABLE_DIR, "13_embedding_generality_verdicts.csv"), row.names = FALSE)
print(ver)
# Summary plot of R(k), only at k values the reference partition actually reached
p <- ggplot(cur[cur$actual_k == cur$k, ], aes(k, mean_cross_seed_ari, colour = embedding_family)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1.5) + facet_wrap(~celltype, scales = "free_x") +
  scale_y_continuous(limits = c(NA, 1)) +
  labs(x = "k", y = "R(k): mean pairwise ARI across embeddings", colour = "Embedding family",
       title = "Across-embedding reproducibility is not scDHA-specific") + pub_theme()
save_fig(p, "13_embedding_generality", FIG_PUB, width = 8, height = 3.8)
