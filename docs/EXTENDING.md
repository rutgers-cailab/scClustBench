# Extending and editing scClustBench

## Add a clustering method

A clustering method is a function `(latent, k, seed)` that returns one integer label per row of
`latent` (cells x dimensions). Register it in `SCCB_METHODS` after sourcing the tool:

```r
source("R/scClustBench.R")
SCCB_METHODS$hclust_ward <- function(latent, k, seed) cutree(hclust(dist(latent), "ward.D2"), k)
res <- sccb_run(obj, methods = c("louvain", "hclust_ward"))
```

The same function re-clusters each bootstrap subsample, so the Jaccard stability metric works for
the new method with no further code. To use it for Stage B as well, pass
`repro_method = "hclust_ward"`.

If the method's expensive work can be shared across k (for example a graph built once and cut
at several resolutions), also register a path function in `SCCB_METHOD_PATHS`:
`function(latent, seed)` returning `function(k)` labels. `sccb_louvain_path()` in
`R/scClustBench.R` is the reference implementation. The results must be identical to calling the
method once per k; `tests/test_speedups.R` checks this for the built-in methods.

## Add an embedding

An embedding is a function `(obj, seed, batch = NULL, ...)` that returns a cells x dimensions
matrix with cell names as row names. `obj` has already been through `sccb_prep()` (normalised,
variable features set, variable genes scaled).

```r
SCCB_EMBEDDINGS$my_vae <- function(obj, seed, batch = NULL, ...) my_vae_latent(obj, seed)
res <- sccb_run(obj, embeddings = c("pca", "my_vae"), repro_embeddings = "my_vae")
```

For Stage B the function is called once per seed, so it should give a genuinely different
embedding for each seed (retraining, resampling or perturbation). An embedding that ignores the
seed gives R(k) = 1 and a meaningless ACCEPT.

Built in: `pca` (Harmony-corrected when `batch` is given), `pca_perturb` (each seed resamples 80%
of the variable genes and runs a seeded randomized PCA), `scdha`, and `scvi` (needs
`RETICULATE_PYTHON`). The PCA-hyperparameter family used in the paper is built in
`pipeline/13_embedding_generality.R`.

## Add a metric

Stage A metrics are computed in `sccb_metrics()` (`R/scClustBench.R`) and ranked by
`sccb_select()` over the columns named in `SCCB_METRICS` (`R/decision_rule.R`). Add the column in
`sccb_metrics()`, oriented so that higher is better, and add its name to `SCCB_METRICS`. Run
`Rscript R/decision_rule.R` and `Rscript tests/smoke_test.R` afterwards.

## Change thresholds

The rule's thresholds are function arguments with the paper's values as defaults:

| Threshold | Argument | Default |
|---|---|---|
| R(k) accept threshold | `tau` in `sccb_run()` / `sccb_verdict()` | 0.8 |
| Tolerance for k\* | `tol` in `sccb_verdict()` | 0.01 |
| Bootstrap stability floor | `jacc_min` | 0.6 |
| Smallest allowed cluster | `min_cells` | 20 |

For the retina analysis every parameter (QC cut-offs, Harmony settings, k ranges, seeds, number
of bootstraps, colours) is set in [config.R](../config.R), with a comment on each.

## Edit the retina analysis

- Each step is one script in `pipeline/`. The header of every script states its inputs, outputs,
  options and the command that runs it. Run a single step with
  `Rscript run_pipeline.R <step>` (see `--list`).
- Steps read and write files only through the paths in `config.R` (`RDS_DIR`, `TABLE_DIR`,
  `PAPER_DIR`), so a step can be rerun on its own once its inputs exist.
- To use other data: change `SAMPLES` in `config.R` and the cell-type marker panels used for
  annotation, and replace `data/annotation/harmony_cluster_annotation.csv` with your own
  cluster-to-cell-type table.
- Figures are made with ggplot2 and saved with `save_fig()` from `R/utils.R` (PNG and PDF).

## Code conventions

- Every script starts with a header block: purpose, paper element, inputs, outputs, run command,
  options.
- Scripts run from the repository root and `source("config.R")` first.
- Long loops use `.par()` (fork-based `parallel::mclapply`; nested calls run serially).
  On Windows set `SLURM_CPUS_PER_TASK=1`.
- Seeds come from `SEED` in `config.R` or from the `seeds` argument.

## Tests

| Test | Command | Checks |
|---|---|---|
| Smoke test | `Rscript tests/smoke_test.R` | the full API on simulated data, about 10 minutes |
| Decision rule | `Rscript R/decision_rule.R` | Stage A and Stage B on small hand-made tables |
| Speed-ups | `Rscript tests/test_speedups.R` | shared-path clustering gives the same partitions as the per-k reference implementation |
| Comment-only edits | `bash tests/check_comments_only.sh <old copy>` | two copies of the code differ only in comments |
