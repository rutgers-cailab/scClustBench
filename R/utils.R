# ============================================================
# R/utils.R
#
# Helper functions shared by the pipeline scripts and by R/scClustBench.R.
# Source it after config.R: several defaults (SEED, MAX_SIL_CELLS, FIG_DPI,
# PUB_BASE_SIZE, RDS_DIR and the cone/rod panels) are read
# from there. R/scClustBench.R sources this file itself and supplies fallback
# values for the few settings it needs, so the API also works without config.R.
#
# Main functions:
#   pub_theme(), save_fig()                figure theme and PDF/PNG export
#   ARI(), NMI()                           partition agreement (labels converted to factors)
#   read10x_robust()                       read a 10x folder whatever its gene column
#   compute_silhouette()                   mean silhouette (subsampled for large data)
#   build_seurat_snn_igraph()              Seurat SNN graph as a weighted igraph
#   find_resolution_for_k()                Louvain/Leiden resolution search to a target k
#   load_population()                      load the whole retina, the rods or one cell class
# ============================================================

# ---- Partition agreement ----
# ARI(a, b), NMI(a, b): adjusted Rand index / normalized mutual information between two label vectors.
# Labels are converted to factors first, because aricode 1.1.0 fails on character labels.
ARI <- function(a, b) aricode::ARI(factor(a), factor(b))
NMI <- function(a, b) aricode::NMI(factor(a), factor(b))

# ---- Figure theme and export ----
# pub_theme(base_size): returns a white, gridless ggplot2 theme used by every figure
pub_theme <- function(base_size = PUB_BASE_SIZE) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text          = element_text(color = "black"),
      plot.title    = element_text(size = base_size + 2, face = "bold", hjust = 0,
                                   margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, hjust = 0, color = "grey30"),
      axis.text     = element_text(size = base_size - 1, color = "black"),
      axis.title    = element_text(size = base_size, color = "black"),
      axis.line     = element_line(color = "black", linewidth = 0.4),
      axis.ticks    = element_line(color = "black", linewidth = 0.3),
      legend.text   = element_text(size = base_size - 1),
      legend.title  = element_text(size = base_size, face = "bold"),
      legend.key.size = unit(0.4, "cm"),
      strip.text    = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      panel.grid    = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin   = margin(8, 8, 8, 8, "mm")
    )
}

# save_fig(plot, filename, dir, width, height): writes <dir>/<filename>.pdf and .png (FIG_DPI)
save_fig <- function(plot, filename, dir, width = 7, height = 5) {
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  png_path <- file.path(dir, paste0(filename, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height,
         device = "pdf")
  ggsave(png_path, plot = plot, width = width, height = height,
         dpi = FIG_DPI, device = "png", bg = "white")
  message("  Saved: ", filename, " (.pdf + .png)")
}


# ---- Data input ----
# read10x_robust(path): reads a 10x matrix folder and returns a genes x cells sparse matrix.
# Gene symbols come from column 2 of features.tsv; folders with one gene column fall back to column 1.
read10x_robust <- function(path) {
  tryCatch(
    Seurat::Read10X(data.dir = path, gene.column = 2),
    error = function(e1) {
      tryCatch(
        Seurat::Read10X(data.dir = path, gene.column = 1),
        error = function(e2) stop("Read10X failed for ", path, "\n", e1$message)
      )
    }
  )
}

# ---- Internal metrics ----
# compute_silhouette(embedding, labels, max_cells, seed): returns the mean silhouette width
# (NA if fewer than two clusters). Above max_cells a random subsample is used, since dist() is O(n^2).
compute_silhouette <- function(embedding, labels, max_cells = MAX_SIL_CELLS,
                                seed = SEED) {
  k <- length(unique(labels))
  if (k < 2) return(NA_real_)

  n <- nrow(embedding)
  if (n > max_cells) {
    set.seed(seed)
    idx <- sample(n, max_cells)
    embedding <- embedding[idx, , drop = FALSE]
    labels    <- labels[idx]
    if (length(unique(labels)) < 2) return(NA_real_)   # a tiny cluster can be absent from the subsample
  }

  d   <- dist(embedding)
  sil <- cluster::silhouette(as.integer(factor(labels)), d)
  mean(sil[, "sil_width"])
}


# ---- Graph construction ----
# build_seurat_snn_igraph(embedding, cell_ids, counts_for_dummy, reduction_name, key_prefix, dims_use):
# returns an undirected weighted igraph of Seurat's SNN graph on `embedding` (cells x dims).
# Using Seurat's SNN construction means Leiden clusters the same graph that the Louvain
# configurations use.
# counts_for_dummy only provides a Seurat object to hold the embedding (its values are unused);
# key_prefix must match the column-name prefix of `embedding`.
build_seurat_snn_igraph <- function(embedding, cell_ids, counts_for_dummy,
                                     reduction_name, key_prefix, dims_use) {
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("igraph", quietly = TRUE)
  requireNamespace("Matrix", quietly = TRUE)

  tmp_obj <- Seurat::CreateSeuratObject(counts = counts_for_dummy[, cell_ids])
  tmp_obj[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = embedding, key = key_prefix, assay = "RNA")
  tmp_obj <- Seurat::FindNeighbors(tmp_obj, reduction = reduction_name,
                                    dims = dims_use, verbose = FALSE)
  snn_mat <- tmp_obj@graphs[["RNA_snn"]]
  rm(tmp_obj); gc(verbose = FALSE)

  # Triplet form exposes (i, j, x); keep the upper triangle so each undirected edge appears once
  snn_mat <- methods::as(snn_mat, "TsparseMatrix")
  keep    <- snn_mat@i < snn_mat@j
  edges   <- cbind(snn_mat@i[keep] + 1L, snn_mat@j[keep] + 1L)
  weights <- snn_mat@x[keep]

  g <- igraph::make_graph(t(edges), n = ncol(snn_mat), directed = FALSE)
  igraph::E(g)$weight <- weights
  igraph::V(g)$name   <- rownames(snn_mat)
  g
}


# ---- Resolution search ----
# find_resolution_for_k(g, target_k, algorithm, res_min, res_max, max_iter, seed): clusters graph g
# with Louvain (modularity) or Leiden (CPM) and searches the resolution for one giving target_k
# clusters. Returns list(labels, resolution, actual_k, exact); if no resolution hits target_k exactly,
# the closest k found is returned with exact = FALSE.
find_resolution_for_k <- function(g, target_k, algorithm = "louvain",
                                  res_min = 1e-6, res_max = 50,
                                  max_iter = 20, seed = SEED) {
  requireNamespace("igraph", quietly = TRUE)

  if (target_k == 1) {
    labs <- rep(1L, igraph::vcount(g))
    names(labs) <- igraph::V(g)$name
    return(list(labels = labs, resolution = 0, actual_k = 1, exact = TRUE))
  }

  cluster_at_res <- function(res) {
    set.seed(seed)
    w <- igraph::E(g)$weight
    if (algorithm == "louvain") {
      comm <- igraph::cluster_louvain(g, weights = w, resolution = res)
    } else {
      # CPM (Constant Potts Model) is Leiden's default and avoids the
      # modularity resolution limit (Fortunato & Barthélemy 2007). On
      # graphs with one dominant community, modularity-Leiden peels off
      # tiny outliers instead of splitting the bulk; CPM gives balanced
      # partitions across the resolution range.
      comm <- igraph::cluster_leiden(g, weights = w,
                                     resolution_parameter = res,
                                     objective_function = "CPM",
                                     n_iterations = 15L)
    }
    memb <- as.integer(factor(igraph::membership(comm)))
    names(memb) <- igraph::V(g)$name
    list(labels = memb, n_clusters = length(unique(memb)))
  }

  # Bisection on a log scale: more clusters than target lowers the upper bound, fewer raises the lower
  best      <- NULL
  best_diff <- Inf
  lo <- res_min; hi <- res_max

  for (iter in seq_len(max_iter)) {
    mid <- sqrt(lo * hi)
    res <- cluster_at_res(mid)
    diff <- abs(res$n_clusters - target_k)
    if (diff < best_diff) {
      best <- list(labels = res$labels, resolution = mid,
                   actual_k = res$n_clusters, exact = (diff == 0))
      best_diff <- diff
    }
    if (diff == 0) break
    if (res$n_clusters < target_k) lo <- mid else hi <- mid
  }

  # The number of clusters is not monotone in resolution, so bisection can miss target_k;
  # fall back to a 60-point log grid over the full range
  if (best_diff > 0) {
    grid <- exp(seq(log(res_min), log(res_max), length.out = 60))
    for (r in grid) {
      res <- cluster_at_res(r)
      diff <- abs(res$n_clusters - target_k)
      if (diff < best_diff) {
        best <- list(labels = res$labels, resolution = r,
                     actual_k = res$n_clusters, exact = (diff == 0))
        best_diff <- diff
      }
      if (diff == 0) break
    }
  }
  best
}

# ---- Populations ----
# load_population(celltype): returns a Seurat object (RNA layers joined) for one population:
#   "Rods": output/rds/02_rods.rds after the cone-contamination filter (the Phase 2 rod input)
#   "ALL" : the whole annotated dataset, output/rds/02_annotated.rds (Phase 1)
#   other : cells of 02_annotated.rds whose reference cell_type equals `celltype` (e.g. "BC")
load_population <- function(celltype) {
  if (celltype == "Rods") {
    x <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
    DefaultAssay(x) <- "RNA"; x[["RNA"]] <- JoinLayers(x[["RNA"]]); x <- NormalizeData(x, verbose = FALSE)
    cone_g <- intersect(CONE_PANEL, rownames(x)); rod_g <- intersect(ROD_PANEL, rownames(x))
    x <- AddModuleScore(x, features = list(cone_g), name = "cone_score", seed = SEED)
    x <- AddModuleScore(x, features = list(rod_g),  name = "rod_score",  seed = SEED)
    x$cone_score <- x$cone_score1; x$rod_score <- x$rod_score1
    # Drop likely cone doublets/contaminants: high cone score that beats the rod score, or a
    # cone score clearly above the rod score (thresholds in config.R)
    cq <- quantile(x$cone_score, CONE_QUANTILE_CUT, na.rm = TRUE)
    contam <- (x$cone_score >= cq & x$cone_score > x$rod_score) | ((x$cone_score - x$rod_score) > CONE_DELTA_CUT)
    return(subset(x, cells = colnames(x)[!contam]))
  }
  obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
  DefaultAssay(obj) <- "RNA"; obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  if (celltype == "ALL") return(obj)
  # Accept either name for the cell-type column ("cell_type" is written by step 02)
  ct_col <- if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "celltype"
  subset(obj, cells = colnames(obj)[obj@meta.data[[ct_col]] == celltype])
}
