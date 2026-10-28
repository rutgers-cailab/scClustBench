#!/bin/bash
# ============================================================
# submit_all.sh
#
# Submits the whole analysis to Slurm as a dependency graph, so steps that do
# not depend on each other run at the same time. Same steps as run_pipeline.R.
#
# Usage (from the repository root, after hpc/slurm/setup_env.sbatch and
# setup_scvi.sbatch have built the two environments):
#   bash hpc/slurm/submit_all.sh
# Optional environment variables:
#   SCCB_RAW_DATA  folder with the six GEO sample folders (see config.R)
#   SCCB_ENV       R environment prefix   (default ~/micromamba/envs/sccb)
#   SCVI_PYTHON    python with scvi-tools (default ~/micromamba/envs/scvi/bin/python)
#   MAIL           address for Slurm failure e-mails (default: none)
#   SCVI_GPU=1     train scVI on a GPU node instead of CPU
# Job logs go to logs/; job IDs are written to logs/submitted.txt.
#
# Set the Slurm partition with SBATCH_PARTITION (and SCVI_PARTITION for GPU jobs) if your
# cluster has no default. Memory and time requests are in the R() calls below.
# ============================================================
set -euo pipefail
cd "$(dirname "$0")/../.."; mkdir -p logs
H=hpc/slurm
M=${MAIL:+--mail-type=FAIL --mail-user=$MAIL}
SCVI=$H/run_scvi_cpu.sbatch; [ "${SCVI_GPU:-0}" = 1 ] && SCVI="--partition=${SCVI_PARTITION:-gpu} $H/run_gpu.sbatch"

sb() { sbatch --parsable --export=ALL $M "$@"; }
# R <job name> <cpus> <memory> <time> [sbatch options] <sbatch script> <pipeline script> [VAR=value ...]
R() { local n=$1 c=$2 m=$3 t=$4; shift 4; sb --job-name=$n --cpus-per-task=$c --mem=$m --time=$t "$@"; }
after() { local IFS=:; echo "--dependency=afterok:$*"; }

# ---- Preprocessing and reference annotation ----
J01=$(R 01_qc 16 180G 12:00:00 $H/run.sbatch 01_qc_harmony_integration.R)
J02=$(R 02_annot 8 160G 06:00:00 $(after $J01) $H/run.sbatch 02_cell_type_annotation.R)
J02b=$(R 02b 4 96G 03:00:00 $(after $J02) $H/run.sbatch 02b_annotation_resolution_sensitivity.R)
J02c=$(R 02c 4 96G 03:00:00 $(after $J02) $H/run.sbatch 02c_harmony_annotation_evidence.R)

# ---- Controls that need no retina data (read later by 13c, 17 and 12) ----
J15=$(R 15 16 96G 48:00:00 $H/run.sbatch 15_external_datasets.R)
J16=$(R 16 16 96G 48:00:00 $H/run.sbatch 16_null_calibration.R)
J09=$(R 09 16 64G 24:00:00 $H/run.sbatch 09_synthetic_validation.R)
J11=$(R 11 4 32G 02:00:00 $(after $J09) $H/run.sbatch 11_synthetic_umaps.R)

# ---- Phase 1: whole retina ----
J03=$(R 03_phase1 32 180G 48:00:00 $(after $J02) $H/run.sbatch 03_phase1_full_dataset_benchmark.R)
J12cp=$(R 12c_phase1 16 96G 72:00:00 $(after $J03) $H/run.sbatch 12c_recompute_jaccard.R DATASET=phase1)
# 12d runs after 12c because both rewrite poc_scores_matrix_4metric.csv
J12d=$(R 12d 4 96G 02:00:00 $(after $J02 $J12cp) $H/run.sbatch 12d_rescore_phase1_ari.R)

# ---- Phase 2: rods. Stage-A bootstraps run as a 56-task array, one partition per task ----
J04=$(R 04_rods 32 180G 36:00:00 $(after $J02) $H/run.sbatch 04_phase2_rod_benchmark.R BENCH_JACC=0)
J12c=$(sb $(after $J04) --array=1-56 $H/jacc_rods.sbatch)
J12cm=$(R 12c_merge 2 16G 01:00:00 --dependency=afterany:$J12c $H/run.sbatch 12c_recompute_jaccard.R DATASET=rods MERGE=1)
J04b=$(R 04b_rods 16 64G 12:00:00 $(after $J04) $H/run.sbatch 04b_chosen_partition_stability.R POP=Rods)
J05=$(R 05 8 160G 12:00:00 $(after $J04) $H/run.sbatch 05_chosen_partition_analysis.R)
J12b=$(R 12b 2 16G 01:00:00 $(after $J05) $H/run.sbatch 12b_marker_technical_audit.R)

# ---- Stage B: the reference scDHA family, one job per seed, then a combine job ----
declare -A J07
for pop in BC Rods ALL; do
  S=()
  for s in 1 2 3 7 42; do S+=($(R 07s_${pop}_$s 16 160G 24:00:00 $(after $J02) $H/run.sbatch 07_positive_control_bc.R CELLTYPE=$pop SEED_ONLY=$s)); done
  J07[$pop]=$(R 07_$pop 32 200G 24:00:00 $(after ${S[@]}) $H/run.sbatch 07_positive_control_bc.R CELLTYPE=$pop)
done
J07b=$(R 07b_bc 16 128G 12:00:00 $(after ${J07[BC]}) $H/run.sbatch 07b_bc_extended_ksweep.R)
J04bb=$(R 04b_bc 16 64G 12:00:00 $(after ${J07[BC]}) $H/run.sbatch 04b_chosen_partition_stability.R POP=BC)
J04c=$(R 04c 8 128G 06:00:00 $(after ${J07[Rods]}) $H/run.sbatch 04c_seed_robustness.R)

# ---- Stage B: the other three embedding families ----
J13bc=$(R 13_bc 32 160G 24:00:00 $(after ${J07[BC]}) $H/run.sbatch 13_embedding_generality.R POP=BC)
J13r=$(R 13_rods 32 180G 24:00:00 $(after ${J07[Rods]}) $H/run.sbatch 13_embedding_generality.R POP=Rods)
J13sbc=$(R 13s_bc 8 96G 24:00:00 $(after ${J07[BC]}) $SCVI 13s_scvi_family.R POP=BC)
J13sr=$(R 13s_rods 8 128G 24:00:00 $(after ${J07[Rods]}) $SCVI 13s_scvi_family.R POP=Rods)
J13c=$(R 13c 2 16G 01:00:00 $(after $J13bc $J13r $J13sbc $J13sr $J16) $H/run.sbatch 13c_combine.R)
J13t=$(R 13t 8 96G 06:00:00 $(after $J13sr $J04) $H/run.sbatch 13t_scvi_rods_inspect.R)
J17c=$(R 17c 8 128G 48:00:00 $(after $J13t) $SCVI 17c_scvi_rods_technical.R)
J17b=$(R 17b 4 64G 02:00:00 $(after $J13t) $H/run.sbatch 17b_scvi_rod_split_audit.R)
J17=$(R 17 2 32G 02:00:00 $(after $J13c $J15 $J16) $H/run.sbatch 17_stageB_gates.R)
J12=$(R 12 2 8G 01:00:00 $(after $J12d $J12cm $J13c ${J07[ALL]} $J09 $J15 $J16) $H/run.sbatch 12_decision_rule.R)

# ---- Sensitivity analyses. Step 14 is split by normalization x distance x population ----
J04e=$(R 04e 8 96G 06:00:00 $(after $J04) $H/run.sbatch 04e_knn_robustness.R)
P14=()
for n in LogNormalize SCT scran; do
  for d in euclidean cosine; do
    P14+=($(R 14a_${n}_$d 32 180G 48:00:00 $(after $J02) $H/run.sbatch 14_preprocessing_sensitivity.R SCCB_BENCH_WORKERS=3 PART=a NORMS=$n DISTS=$d TAG=a_${n}_$d))
  done
  for c in BC Rods; do
    P14+=($(R 14b_${c}_$n 16 160G 48:00:00 $(after $J02) $H/run.sbatch 14_preprocessing_sensitivity.R PART=b CELLTYPES=$c NORMS=$n TAG=b_${c}_$n))
  done
done
J14=$(R 14c 2 16G 01:00:00 $(after ${P14[@]}) $H/run.sbatch 14_preprocessing_sensitivity.R PART=c)

# ---- Paper figures (output/paper/) once every table exists ----
ALL_TABLES="$J02b $J02c $J12d $J04b $J04bb $J05 $J12b $J07b $J04c $J17c $J17b $J17 $J12 $J04e $J14 $J15 $J16 $J11"
JF=$(R figures 4 96G 04:00:00 $(after $ALL_TABLES) $H/run.sbatch FIGURES)

echo "01=$J01 02=$J02 03=$J03 04=$J04 07=${J07[*]} 13c=$J13c 12=$J12 14=$J14 figures=$JF" | tee logs/submitted.txt
