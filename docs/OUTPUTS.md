# Outputs

| Folder | Contents | In the repository |
|---|---|---|
| `output/paper/` | paper figures (`Figure_1` to `Figure_7`, `Supplementary_Figure_S1` to `S6`, PNG and PDF) and the `source_*.csv` data behind them | PNG and CSV (PDFs are regenerated) |
| `output/tables/` | result tables for rods, bipolar cells, sensitivity analyses and controls | yes |
| `output_full_dataset/tables/` | Phase 1 (whole retina) benchmark tables | yes |
| `output/rds/`, `output_full_dataset/rds/` | Seurat objects, benchmark bundles and cached embeddings | no (large; see below) |
| `*_latent.csv`, `poc_harmony_*.csv`, `*partitions*.csv` | embedding coordinates, the Phase 1 UMAP and partition caches | no (large) |

Steps that write RDS files (none of them tracked):

| Step | RDS files |
|---|---|
| 01 | `output/rds/01_integrated.rds` |
| 02 | `output/rds/02_annotated.rds`, `02_rods.rds` |
| 03 | `output_full_dataset/rds/poc_benchmarking.rds` |
| 04 | `output/rds/03_benchmarking.rds` |
| 04c | `output/rds/04c_seed_robustness.rds` |
| 05 | `output/rds/04_rod_subclusters.rds` |
| 07 | `output/rds/07_<pop>_positive_control.rds`, `07_<pop>_seed<s>.rds` |
| 09 | `output/rds/09_synthetic_validation.rds` |
| 13s | `output/rds/13s_scvi_<pop>_seed<s>.rds`, `13s_scvi_<pop>_latents.rds` |
| 15 | `output/tables/15_external/<case>/meta.rds` |
| 17c | `output/tables/17c_parts/<condition>_seed<s>.rds` |

Most table names start with the id of the step that writes them (`12_*` from step 12, `13t_*`
from 13t, and so on). The exceptions:

| Table names | Written by |
|---|---|
| `output_full_dataset/tables/poc_*` | step 03 (Phase 1 benchmark); `poc_scores_matrix_4metric.csv` is rewritten by 12c (`jacc_min`) and 12d (ARI, NMI); `poc_harmony_pca_umap.csv` by step 06 |
| `output/tables/03_*` | step 04 (rod benchmark); `03_scores_matrix_4metric.csv` is rewritten by 12c (`jacc_min`) |
| `output/tables/04_*` | step 05 (illustrative rod partition); `04_harmony_scdha_umap.csv` and `04_deg_counts.csv` are recomputed by step 06 if absent |
| `output/tables/07_bc_*_extended.csv` | step 07b |
| `output/tables/13_*` | steps 13 (`*_bc.csv`, `*_rods.csv`) and 13c (combined tables); 13c also writes `13s_scvi_curves.csv` |
| `output_full_dataset/tables/12d_*` | step 12d |
| `output/tables/_knn_robustness_summary.csv` | step 04e |
| `output/tables/12c_parts/`, `output_full_dataset/tables/12c_parts/` | step 12c (per-partition results of the Slurm array) |

## Paper figures

| Figure | Made by | Main inputs |
|---|---|---|
| Figure 1 | `00_pipeline_schematic.py` | none (drawing) |
| Figure 2 | `06` | `output/rds/02_annotated.rds` (QC, clusters, markers, cell types) |
| Figure 3 | `06` | `output_full_dataset/tables/poc_scores_matrix_4metric.csv`, `poc_all_cluster_assignments.csv`, `poc_harmony_pca_umap.csv` (or `poc_harmony_pca_latent.csv`), `output/tables/07_all_stability_by_k.csv` |
| Figure 4 | `08b` | `output/rds/07_bc_positive_control.rds`, `12_verdicts.csv`, `07_bc_stability_by_k_extended.csv` (or `07_bc_stability_by_k.csv`) |
| Figure 5 | `06` | `03_scores_matrix_4metric.csv`, `03_all_cluster_assignments.csv`, `03_harmony_scdha_latent.csv`, `output/rds/02_rods.rds`, `04_harmony_scdha_umap.csv`, `04_deg_counts.csv` |
| Figure 6 | `08` | `07_bc_stability_by_k.csv`, `07_bc_stability_by_k_extended.csv`, `07_rods_stability_by_k.csv` |
| Figure 7 | `10` | `09_synthetic_validation.csv`, `11_synthetic_umap_coords.csv` |
| Supplementary Figure S1 | `06` | `output/rds/02_annotated.rds` (QC and batch mixing) |
| Supplementary Figure S2 | `06` | `output_full_dataset/tables/poc_scores_matrix_4metric.csv` (Phase 1 ARI across k, partition sizes) |
| Supplementary Figure S3 | `13f` | `13_embedding_generality_curves.csv` (bipolar k sweep) |
| Supplementary Figure S4 | `13f` | `13_embedding_generality_curves.csv` (four embedding families) |
| Supplementary Figure S5 | `06` | `output/rds/02_rods.rds`, `03_all_cluster_assignments.csv` (rod partition by sample and time point) |
| Supplementary Figure S6 | `04d` | `04b_chosen_partition_stability.csv`, `04c_seed_ari_matrix.csv`, `04c_seed_sizes.csv` |

Step 06 draws all of its figures in one run and reads every input listed for Figures 2, 3, 5 and
S1 to S3, so it needs `output/rds/02_annotated.rds`, `output/rds/02_rods.rds`,
`output/tables/03_harmony_scdha_latent.csv` and `output_full_dataset/tables/poc_harmony_pca_umap.csv`
(or `poc_harmony_pca_latent.csv`), none of which are tracked. Step 08b needs
`output/rds/07_bc_positive_control.rds`. Steps `fig1`, `08`, `10`, `04d` and `13f` run from the
tracked tables alone.

## Key result tables

| Question | Table |
|---|---|
| Final verdict and k\* for each population and embedding family | `output/tables/12_verdicts.csv`, `13_embedding_generality_verdicts.csv` |
| R(k) curves for all four families | `output/tables/13_embedding_generality_curves.csv` |
| Stage A selection (Phase 1, rods) | `output/tables/12_selection_phase1.csv`, `12_selection_rods.csv` |
| Partition reported at k\* | `output/tables/12_partition_at_kstar.csv` |
| Verdicts under alternative tau, k\* rule and gates | `output/tables/12_verdicts_tau_sensitivity.csv`, `12_selection_sensitivity.csv`, `17_stageB_gated_verdicts.csv` |
| Phase 1 scores for every partition | `output_full_dataset/tables/poc_scores_matrix_4metric.csv` |
| Rod scores for every partition | `output/tables/03_scores_matrix_4metric.csv` |
| Reference annotation (cell to type) | `output/tables/02_reference_annotation.csv` |
| Annotation evidence per Harmony cluster | `output/tables/02c_*.csv` |
| Rod partition markers and technical audit | `output/tables/04_all_markers.csv`, `12b_marker_technical_audit.csv` |
| scVI rod split and technical tests | `output/tables/13t_*.csv`, `17b_*.csv`, `17c_*.csv` |
| Normalization x distance | `output/tables/14_preprocessing_verdicts.csv`, `14_preprocessing_phase1_selection.csv` |
| External datasets | `output/tables/15_external_verdicts.csv`, `15_external/<case>/sccb_*.csv` |
| Null calibration | `output/tables/16_null_calibration_*.csv` |
| Synthetic recovery | `output/tables/09_synthetic_validation.csv` |

`output/tables/_sex_artifact_summary.csv` (Xist and Y-gene detection per sample and per rod
cluster) was computed on a broader rod selection and has no step in `pipeline/`; it is kept
because the paper's technical audit table quotes its per-sample values.
