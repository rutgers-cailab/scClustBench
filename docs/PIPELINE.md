# The retina analysis, step by step

`run_pipeline.R` runs these steps in this order; `Rscript run_pipeline.R --list` prints them with
their environment variables. Each step is one script in `pipeline/`. The header of each script
lists its exact inputs, outputs and options.

Phase 1 is the whole retina, where the cell classes are known (positive control). Bipolar cells
are the subtype-level positive control. Rods are the population tested for unknown sub-structure.

## Preprocessing and reference labels

| Step | Script | What it does |
|---|---|---|
| 01 | `01_qc_harmony_integration.R` | Reads the six 10x samples; QC filters, scDblFinder doublet removal, SCTransform, PCA, Harmony integration, clustering at resolution 0.6 |
| 02 | `02_cell_type_annotation.R` | Applies the hand annotation of the 34 Harmony clusters (`data/annotation/`) as the reference labels; also records the automated module-score label; writes the rod subset |
| 02b | `02b_annotation_resolution_sensitivity.R` | Agreement of the annotation across clustering resolutions 0.2 to 1.0 |
| 02c | `02c_harmony_annotation_evidence.R` | Per-cluster evidence behind the hand annotation: marker panels, top markers, sample mix, QC |

## Controls that need no retina data

These steps run early because steps 13c, 17 and 12 read their tables.

| Step | Script | What it does |
|---|---|---|
| 15 | `15_external_datasets.R` | The rule on external public datasets with independent labels (`CASE=` selects one) |
| 16 | `16_null_calibration.R` | Splatter null and structured simulations; operating characteristics of tau. `16_null_calibration_curves.csv` is also its resume checkpoint, so on a fresh clone (where that table is present) only the summary tables are recomputed; `FRESH=1` recomputes every simulation |
| 09, 11 | `09_synthetic_validation.R`, `11_synthetic_umaps.R` | Splatter ground-truth recovery and UMAPs |

## Phase 1: whole retina

| Step | Script | What it does |
|---|---|---|
| 03 | `03_phase1_full_dataset_benchmark.R` | Seven configurations x k = 2 to 12, five metrics per partition, ARI/NMI against the reference (post hoc) |
| 12c_p1 | `12c_recompute_jaccard.R` (`DATASET=phase1`) | Bootstrap Jaccard of every Phase 1 partition (20 resamples of 80% of cells) |
| 12d | `12d_rescore_phase1_ari.R` | Rescores Phase 1 ARI/NMI against the reference labels |

The Phase 1 R(k) curve is step 07_all, run with the other scDHA curves below.

## Rods: Stage A and the illustrative partition

| Step | Script | What it does |
|---|---|---|
| 04 | `04_phase2_rod_benchmark.R` | Rod benchmark: seven configurations x k = 2 to 9, five metrics |
| 12c_rods | `12c_recompute_jaccard.R` (`DATASET=rods`) | Bootstrap Jaccard of every rod partition |
| 04b | `04b_chosen_partition_stability.R` (`POP=Rods`) | 50-resample bootstrap stability of the scDHA + Louvain rod partition |
| 05 | `05_chosen_partition_analysis.R` | Applies the illustrative k = 3 rod partition; markers, composition per sample and time point, UMAP |
| 12b | `12b_marker_technical_audit.R` | Classifies the rod cluster markers (ribosomal, sex-linked, stress, lncRNA) |

Its seed robustness (04c) and neighbour-count robustness (04e) run later, after 07_rods and in
the sensitivity analyses.

## Stage B: reproducibility across embeddings

| Step | Script | What it does |
|---|---|---|
| 07_bc, 07_rods, 07_all | `07_positive_control_bc.R` | R(k) across five scDHA embeddings for bipolar cells, rods and the whole retina (Phase 1) |
| 07b | `07b_bc_extended_ksweep.R` | Bipolar R(k) extended to k = 18 |
| 04b_bc | `04b_chosen_partition_stability.R` (`POP=BC`) | Bootstrap stability of the bipolar partitions at k = 9 and 16 |
| 04c | `04c_seed_robustness.R` | The illustrative rod partition under five scDHA seeds (reuses the 07_rods embeddings) |
| 13_bc, 13_rods | `13_embedding_generality.R` | R(k) for the scDHA family (the step 07 embeddings) and the two scDHA-independent PCA families: gene resampling and hyperparameters. With `POP` set the scVI family is left to 13s |
| 13s_bc, 13s_rods | `13s_scvi_family.R` | R(k) under five scVI seeds (needs `RETICULATE_PYTHON`) |
| 13c | `13c_combine.R` | Merges the four families into one curve and verdict table (tau from the step 16 calibration when higher than 0.8) |
| 13t | `13t_scvi_rods_inspect.R` | Describes the scVI rod k = 3 partition (depth, ribosomal fraction, markers, samples) |
| 17c | `17c_scvi_rods_technical.R` | scVI rod R(k) after depth thinning, with technical covariates, and without technical genes |
| 17b | `17b_scvi_rod_split_audit.R` | The 13t partition by sample, and its agreement across seeds and with Stage A partitions |
| 17 | `17_stageB_gates.R` | Stage B verdicts with and without the validity gates |
| 12 | `12_decision_rule.R` | Applies the full rule to every population and writes the final verdicts and selections |

## Sensitivity analyses

| Step | Script | What it does |
|---|---|---|
| 04e | `04e_knn_robustness.R` | The illustrative rod partition under different neighbour counts |
| 14, 14_tables | `14_preprocessing_sensitivity.R` | Normalization (LogNormalize, SCT, scran) x graph distance (Euclidean, cosine) for Phase 1, bipolar cells and rods; `PART=c` combines the parts |

## Figures

| Step | Script | Paper figures |
|---|---|---|
| fig1 | `00_pipeline_schematic.py` | Figure 1 |
| 06 | `06_publication_figures.R` | Figures 2, 3, 5; Supplementary Figures S1, S2 and S5 |
| 08b | `08b_bc_positive_control_figure.R` | Figure 4 |
| 08 | `08_discovery_comparison_figure.R` | Figure 6 |
| 10 | `10_synthetic_figure.R` | Figure 7 |
| 04d | `04d_stability_figures.R` | Supplementary Figure S6 |
| 13f | `13f_supp_figures.R` | Supplementary Figures S3, S4 |

Steps `fig1`, `08`, `10`, `04d` and `13f` read only tracked tables, so they can be rerun quickly
after any change. Step 06 also needs `output/rds/02_annotated.rds`, `output/rds/02_rods.rds`,
`output/tables/03_harmony_scdha_latent.csv` and `output_full_dataset/tables/poc_harmony_pca_umap.csv`
(computed from `poc_harmony_pca_latent.csv` if absent), and step 08b needs
`output/rds/07_bc_positive_control.rds`; none of these files are tracked.

## Runtime (one 32-core node)

| Steps | Time |
|---|---|
| 01, 02 | about 1 hour |
| 03, 04 | half a day to two days each |
| 12c | several hours (one array task per partition on Slurm) |
| 07, 13, 13s, 17c | several hours each; scVI is faster on a GPU |
| 14, 15, 16 | one to two days each |
| figures | minutes |
