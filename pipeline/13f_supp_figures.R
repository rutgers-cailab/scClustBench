# ============================================================
# pipeline/13f_supp_figures.R
#
# Draws Supplementary Figures S3 and S4 from the combined Stage B curve table.
# S5 shows R(k) and the minimum pairwise ARI for bipolar cells in the scDHA
# family over k = 2 to 18, with k* and the k values excluded by the Stage B
# validity gates marked. S6 shows R(k) for bipolar cells and rods in all four
# embedding families. Only the CSV table is needed (no .rds files), so the
# script runs on any machine with ggplot2.
#
# Paper:    Supplementary Figures S3 and S4
# Inputs:   output/tables/13_embedding_generality_curves.csv (written by 13c_combine.R)
# Outputs:  output/paper/Supplementary_Figure_S3.{pdf,png}
#           output/paper/Supplementary_Figure_S4.{pdf,png}
#           (or in the folder given as the first argument)
# Run:      Rscript run_pipeline.R 13f   or   Rscript pipeline/13f_supp_figures.R [output_dir]
# ============================================================
source("config.R"); source("R/utils.R"); source("R/decision_rule.R")
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args)) args[1] else PAPER_DIR
TAU  <- 0.8

# ---- Curves and per-k validity flags ----
# One group per population x embedding family
cur <- read.csv(file.path(TABLE_DIR, "13_embedding_generality_curves.csv"))
grp <- interaction(cur$celltype, cur$embedding_family, drop = TRUE)

# Stage B validity gates, as in sccb_verdict(): achieved k equals target k, every
# embedding of the family reached k, smallest reference cluster >= 20 cells.
cur$valid <- cur$actual_k == cur$k &
  cur$n_embeddings_at_k == ave(cur$n_embeddings_at_k, grp, FUN = max) &
  vapply(cur$seed1_sizes, .min_size, numeric(1)) >= 20
ver <- do.call(rbind, lapply(split(cur, grp), function(d)
  cbind(celltype = d$celltype[1], embedding_family = d$embedding_family[1], sccb_verdict(d, tau = TAU))))
# Self-check: the per-k flags reproduce the rule's own count of valid k, and k* is a valid k
stopifnot(all(tapply(cur$valid, grp, sum)[rownames(ver)] == ver$n_valid_k))
bc <- cur[cur$celltype == "BC" & cur$embedding_family == "scdha", ]
k_star <- ver$k_star[ver$celltype == "BC" & ver$embedding_family == "scdha"]
message("bipolar scDHA valid k: ", paste(bc$k[bc$valid], collapse = ", "), "; k* = ", k_star)
stopifnot(k_star %in% bc$k[bc$valid])

# Point shape encodes the gate status in both figures (filled = valid, open = excluded)
cur$status <- factor(ifelse(cur$valid, "Valid k", "Excluded by Stage-B gates"),
                     levels = c("Valid k", "Excluded by Stage-B gates"))
bc <- cur[cur$celltype == "BC" & cur$embedding_family == "scdha", ]
shape_scale <- scale_shape_manual(values = c("Valid k" = 16, "Excluded by Stage-B gates" = 1), name = NULL)
tau_line <- geom_hline(yintercept = TAU, linetype = "dashed", colour = "black", linewidth = 0.3)

# ---- Supplementary Figure S3: bipolar cells, scDHA family, k = 2 to 18 ----
# Long format so the mean (R(k)) and minimum pairwise ARI share one colour legend
long <- rbind(data.frame(k = bc$k, ari = bc$mean_cross_seed_ari, stat = "Mean pairwise ARI, R(k)"),
              data.frame(k = bc$k, ari = bc$min_cross_seed_ari,  stat = "Minimum pairwise ARI"))
s5 <- ggplot(bc, aes(k, mean_cross_seed_ari)) +
  tau_line +
  geom_vline(xintercept = k_star, linetype = "dotted", linewidth = 0.3) +
  annotate("text", x = k_star - 0.2, y = 0.02, label = paste0("k* = ", k_star), hjust = 1, size = 3) +
  annotate("text", x = 2, y = TAU + 0.03, label = "tau == 0.8", parse = TRUE, hjust = 0, size = 3) +   # plotmath: the pdf device cannot encode a UTF-8 tau
  geom_line(data = long, aes(k, ari, colour = stat), linewidth = 0.5) +
  geom_point(aes(shape = status), colour = "#0072B2", fill = "white", size = 2.2, stroke = 0.8) +
  scale_colour_manual(values = c("Mean pairwise ARI, R(k)" = "#0072B2", "Minimum pairwise ARI" = "#9EC9E2"), name = NULL) +
  shape_scale +
  scale_x_continuous(breaks = 2:18) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Cluster number k", y = "Pairwise ARI across five scDHA embeddings") +
  pub_theme() + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"))
save_fig(s5, "Supplementary_Figure_S3", OUT, width = 6.5, height = 4.6)

# ---- Supplementary Figure S4: bipolar cells and rods, four embedding families ----
# Display names for the embedding families; the stopifnot catches a family missing from this map
fam <- c(scdha = "scDHA", pca_perturb = "PCA gene resampling", pca_hyper = "PCA hyperparameters", scvi = "scVI")
stopifnot(all(cur$embedding_family %in% names(fam)))
cur$family <- factor(fam[cur$embedding_family], levels = fam)
cur$panel  <- factor(c(BC = "Bipolar cells", Rods = "Rods")[cur$celltype], levels = c("Bipolar cells", "Rods"))
# Families are dodged along k so points that coincide (e.g. bipolar k = 13 to 18) stay readable
dodge <- position_dodge(width = 0.5)
s6 <- ggplot(cur, aes(k, mean_cross_seed_ari, colour = family, group = family)) +
  tau_line +
  geom_line(linewidth = 0.5, position = dodge) +
  geom_point(aes(shape = status), fill = "white", size = 1.9, stroke = 0.7, position = dodge) +
  facet_wrap(~panel, scales = "free_x") +
  scale_colour_manual(values = setNames(c("#0072B2", "#E69F00", "#D55E00", "#009E73"), fam), name = "Embedding family") +
  shape_scale +
  scale_x_continuous(breaks = function(l) seq(2, floor(l[2]), 2)) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Cluster number k", y = "R(k): mean pairwise ARI across embeddings") +
  pub_theme() + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"))
save_fig(s6, "Supplementary_Figure_S4", OUT, width = 8, height = 4.8)
