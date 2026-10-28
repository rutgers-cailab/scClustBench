# Software environment

The analysis used R 4.5.3 with these package versions (`hpc/slurm/setup_env.sbatch`):

| Package | Version | Package | Version |
|---|---|---|---|
| Seurat | 5.5.1 | SeuratObject | 5.4.0 |
| harmony | 2.0.5 | scDHA | 1.2.3 |
| torch (R) | 0.17.0 | ROGUE | 1.0 |
| mclust | 6.1.3 | aricode | 1.0.3 |
| igraph | 2.3.3 | leidenbase | 0.1.36 |
| uwot | 0.2.5 | FNN | 1.1.4.1 |
| scDblFinder | 1.24.10 | splatter | 1.34.0 |
| scran | 1.38.1 | scRNAseq | 2.24.0 |
| reticulate | 1.47.0 | Matrix | 1.7.6 |
| cluster | 2.1.8.3 | ggplot2 | 4.0.3 |
| patchwork | 1.3.2 | dplyr / tidyr | 1.2.1 / 1.3.2 |

Python (scVI): Python 3.11, scvi-tools 1.4.2, torch 2.11.0.

GitHub packages are pinned to commits in `setup_env.sbatch`: scDHA `duct317/scDHA@ac9060b70768`
and ROGUE `PaulingLiu/ROGUE@6e1c8f9ee419`.

The paper used aricode 1.0.3. Every ARI and NMI call, in `R/` and in `pipeline/`, goes through
wrappers that convert the labels to factors first (`ARI()` and `NMI()` in `R/utils.R`,
`.ari()` and `.nmi()` in `R/scClustBench.R`), so aricode 1.1.0, which fails on character
labels, works as well.

Figure 1 is drawn with Python matplotlib (`pipeline/00_pipeline_schematic.py`).
