# Using scClustBench on your own data

Everything below runs from the repository root. The tool itself is two files,
`R/scClustBench.R` and `R/decision_rule.R`, and `source("R/scClustBench.R")` loads both
(plus the helpers in `R/utils.R`).

## Input

A Seurat object (v5) with raw counts in the `RNA` assay. Subset it to the population you want to
test before calling the tool:

- the whole dataset, to choose the number of clusters and the clustering configuration; or
- one cell type, to ask whether it contains reproducible subtypes.

If the cells come from several samples or donors, pass the metadata column as `batch`. The PCA
embedding is then corrected with Harmony.

## One call: `sccb_run()`

```r
source("R/scClustBench.R")
res <- sccb_run(obj,
                batch            = "sample",          # metadata column, or NULL
                k_range          = 2:10,              # numbers of clusters to test
                embeddings       = "pca",             # Stage A embeddings
                methods          = c("kmeans", "gmm", "louvain", "leiden"),
                repro_embeddings = c("scdha", "pca_perturb"),   # Stage B families
                repro_method     = "louvain",         # clustering used to compute R(k)
                seeds            = c(1, 2, 3, 7, 42), # one embedding per seed in Stage B
                tau              = 0.8,               # accept threshold on R(k)
                jacc_min         = 0.6,               # Stage A bootstrap stability floor
                min_cells        = 20,                # smallest allowed cluster
                n_boot           = 20,                # bootstrap resamples for the Jaccard
                truth            = NULL,              # optional labels, scored as ARI/NMI only
                out_dir          = "sccb_out")        # writes sccb_*.csv when set
```

| Argument | Default | Notes |
|---|---|---|
| `norm` | `"LogNormalize"` | also `"SCT"` or `"scran"` |
| `nfeatures` | 3000 | highly variable genes |
| `rogue` | `TRUE` | compute ROGUE purity (needs the ROGUE package) |
| `repro_embeddings` | `"scdha"` if installed, else `"pca_perturb"` | list several to compare families |

### Output

| Element | File in `out_dir` | Contents |
|---|---|---|
| `res$scores` | `sccb_scores.csv` | one row per (embedding_method, k): balance, jacc_min, marker_spec, rogue, silhouette, achieved k, cluster sizes |
| `res$selection` | `sccb_selection.csv` | Stage A: admissible flag and mean-percentile score |
| `res$curve` | `sccb_curve.csv` | Stage B: R(k) per family, number of embeddings reaching k, smallest cluster |
| `res$verdict` | `sccb_verdict.csv` | per family: ACCEPT or REJECT, `k_star`, `k_star_smallest`, `k_star_at_edge`, `max_R` |
| `res$recommended` | `sccb_recommended.csv` | the answer: verdict, k\*, chosen configuration, and a note when families disagree |
| `res$partition` | `sccb_partition.csv` | cell to cluster for the chosen configuration (absent on REJECT) |

## Step by step

`sccb_run()` is a thin wrapper around these functions, which can be called separately (for
example to cache the benchmark and try several thresholds):

| Function | Does |
|---|---|
| `sccb_prep(obj, norm, nfeatures)` | normalise, find HVGs, scale (once; every embedding reuses it) |
| `sccb_embed(obj, embedding, seed, batch)` | one embedding, cells x dimensions |
| `sccb_benchmark(obj, embeddings, methods, k_range, ...)` | every partition and its five metrics |
| `sccb_select(scores, jacc_min, min_cells)` | Stage A screen and ranking |
| `sccb_reproducibility(obj, embedding, seeds, k_range, ...)` | R(k) for one family |
| `sccb_verdict(curve, tau)` | Stage B verdict and k\* |
| `sccb_recommend(scores, verdict)` | the reported configuration at k\* |

`R/decision_rule.R` has a self-check: `Rscript R/decision_rule.R` prints `self-check OK`.

## Runtime

Cost grows with cells x methods x k x `n_boot`. Graph methods (Louvain, Leiden) repeat their
resolution search on every bootstrap resample, so they dominate. Suggestions:

- `n_boot = 5` for a first look; 20 (the default) for reported results.
- Parallel loops use `SLURM_CPUS_PER_TASK` workers, or `N_CORES` from `config.R` when it is
  sourced. Set `SCCB_BENCH_WORKERS` to cap the number of workers in the benchmark if memory is
  tight.
- Above about 20,000 cells, run on a cluster ([HPC.md](HPC.md)).

## Interpreting the verdict

- **ACCEPT** in a family: partitions from independent embeddings agree (max R(k) >= tau). k\* is
  the largest k within 0.01 of the best R.
- **REJECT**: no k reaches tau. The population has no sub-structure that is reproducible across
  embeddings, even if one partition looks stable under bootstrap.
- Families disagree: report the structure as embedding-dependent.
- `k_star_at_edge = TRUE`: extend `k_range`.
- A family whose embeddings do not all reach a given k, or whose reference partition has a
  cluster under `min_cells`, is not counted at that k (validity gates).
