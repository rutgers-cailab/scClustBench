# ============================================================
# pipeline/12b_marker_technical_audit.R
#
# Checks whether the markers of the illustrative rod partition are dominated by genes that commonly
# reflect technical rather than biological variation: sex-linked genes (confounded with sample),
# dissociation and immediate-early stress genes, ribosomal protein genes, and nuclear lncRNAs (Malat1 and
# similar). For the top 50 and top 200 markers of each cluster (ranked by specificity score), reports the
# percentage in each gene class and the share of total specificity score they carry.
#
# Paper:    Technical audit of the illustrative rod partition markers (supporting table)
# Inputs:   output/tables/04_all_markers.csv (step 05)
# Outputs:  output/tables/12b_marker_technical_audit.csv
# Run:      Rscript run_pipeline.R 12b
# ============================================================
source("config.R")
# Gene-symbol regular expressions for each technical gene class
TECH <- list(
  sex        = "^(Xist|Tsix|Ddx3y|Eif2s3y|Kdm5d|Uty)$",
  stress_IEG = "^(Fos|Fosb|Jun|Junb|Jund|Atf3|Egr[1-3]|Ier[235]|Hspa1[ab]|Hsp90aa1|Dusp1|Klf[246]|Nr4a1|Zfp36|Socs3|Btg2|Cebpd|Gadd45[bg])$",
  # ribosomal protein genes, paralogs and pseudogenes; not Rps6k* kinases or Rps19bp1
  ribosomal  = "^(Rp[ls][0-9]+[a-z]?[0-9]*(l[0-9]*)?(-ps[0-9]+)?|Rps27rt|Rplp[0-2]|Rpsa|Fau|Uba52)$",
  nuclear_lnc = "^(Malat1|Neat1|Meg3|Gm[0-9]+)$")
m <- read.csv(file.path(TABLE_DIR, "04_all_markers.csv"))
# significant positive markers only
m <- m[m$p_val_adj < 0.05 & m$avg_log2FC > 0.25, ]
m$score <- pmax(0, (m$pct.1 - m$pct.2) * abs(m$avg_log2FC))   # marker-specificity contribution (clipped at 0)
rows <- list()
for (cl in sort(unique(m$cluster))) for (n in c(50, 200)) {
  top <- head(m[m$cluster == cl, ][order(-m$score[m$cluster == cl]), ], n)
  # logical matrix: marker gene (rows) x technical class (columns)
  cls <- sapply(TECH, function(p) grepl(p, top$gene))
  rows[[length(rows) + 1]] <- data.frame(cluster = cl, top_n = n,
    t(round(100 * colMeans(cls), 1)), any_technical_pct = round(100 * mean(rowSums(cls) > 0), 1),
    tech_share_of_specificity = round(100 * sum(top$score[rowSums(cls) > 0]) / sum(top$score), 1),
    top10 = paste(head(top$gene, 10), collapse = ";"))
}
out <- do.call(rbind, rows); print(out[, -ncol(out)])
write.csv(out, file.path(TABLE_DIR, "12b_marker_technical_audit.csv"), row.names = FALSE)
