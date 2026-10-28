# Running on a Slurm cluster

The scripts in `hpc/slurm/` run the analysis on a Slurm cluster. If your cluster has no default
partition, set it with `export SBATCH_PARTITION=<name>` before submitting (and `SCVI_PARTITION`
for GPU jobs). Memory and time requests are in `submit_all.sh`.

## 1. Build the environments (once)

```bash
sbatch hpc/slurm/setup_env.sbatch     # R 4.5 + all R packages, in a micromamba env "sccb"
sbatch hpc/slurm/setup_scvi.sbatch    # Python 3.11 + scvi-tools 1.4.2, env "scvi"
```

Both install into `~/micromamba` without administrator rights. `setup_env.sbatch` ends with
`ENV_OK` and writes `~/sccb_env_sessionInfo.txt`.

## 2. Put the data in place

Copy the six GEO sample folders to the cluster and set `SCCB_RAW_DATA` (see the README).

## 3. Submit

From the repository root:

```bash
export SCCB_RAW_DATA=/path/to/retina_samples
bash hpc/slurm/submit_all.sh
```

This submits every step as its own job with `afterok` dependencies, so independent steps run at
the same time:

- 07 runs one job per scDHA seed (`SEED_ONLY=<seed>`), then a combine job.
- 12c for rods runs as a 56-task array, one partition per task (`jacc_rods.sbatch`).
- 14 runs one job per normalization x distance x population, then `PART=c` combines them.
- The figure job (`FIGURES`) waits for every table.

Logs go to `logs/<job name>_<job id>.log`, and the job IDs to `logs/submitted.txt`.

| Variable | Default | Meaning |
|---|---|---|
| `SCCB_RAW_DATA` | parent of the repository | folder with the six samples |
| `SCCB_ENV` | `~/micromamba/envs/sccb` | R environment |
| `SCVI_PYTHON` | `~/micromamba/envs/scvi/bin/python` | Python with scvi-tools |
| `SCVI_GPU` | 0 | 1 trains scVI on a GPU node (`run_gpu.sbatch`) instead of CPU |
| `MAIL` | none | address for failure e-mails |
| `SBATCH_PARTITION` | cluster default | partition for all jobs |
| `SCVI_PARTITION` | `gpu` | partition for GPU scVI jobs |

## Running one step

```bash
mkdir -p logs
sbatch --cpus-per-task=16 --mem=96G --time=12:00:00 hpc/slurm/run.sbatch 04_phase2_rod_benchmark.R
sbatch --cpus-per-task=8 --mem=128G --time=24:00:00 hpc/slurm/run_scvi_cpu.sbatch 13s_scvi_family.R POP=Rods
```

Arguments after the script name are exported as environment variables.

## Notes

- CPU scVI: `run_scvi_cpu.sbatch` sets `OPENBLAS_NUM_THREADS=1`. A multi-threaded OpenBLAS in
  the same process as torch slowed training more than tenfold.
- Parallel loops inside a job use `SLURM_CPUS_PER_TASK` workers.
- Step 16 resumes from `output/tables/16_null_calibration_curves.csv`, which is tracked, so on a
  fresh clone its job only recomputes the summary tables. To recompute the simulations, add
  `FRESH=1` after `16_null_calibration.R` in its line of `submit_all.sh`.
- To keep jobs off particular nodes, add `--exclude=<nodes>` to the `sb()` call in `submit_all.sh`.
