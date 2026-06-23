# Simulation Code for More-Powerful K-Means Inference

This directory contains the research code for the more-powerful selective
k-means experiments. The package API in `R/` exports the known-variance tests;
the unknown-variance R-fiber code, large sweeps, real-data examples, and figure
generation live here.

Run scripts from the package root:

```bash
Rscript sims/<script>.R
```

Generated outputs go to `sims/results/`, which is intentionally git-ignored.
The tracked source of truth is the code in `sims/` plus the manifest in
`sims/RESULTS_MANIFEST.md`.

## Method Map

| test | variance | status | main code |
|---|---|---|---|
| union, more powerful | known sigma | valid, package API | `KmeansInference::kmeans_inference_union` |
| path, Chen-Witten | known sigma | valid, package API | `KmeansInference::kmeans_inference` |
| path, Yun-Barber style | unknown sigma | valid F-pivot post-processing | `unknown_var_union.R` |
| union/path on R-fiber | unknown sigma | valid studentized test | `kmeans_union_unknownvar.R` |
| exact R-fiber solver | unknown sigma | valid, any K | `kmeans_union_unknownvar_exact.R`, `rfiber_exact.R` |

Important: the valid unknown-variance union test must be computed on the
studentized F-pivot's own R-fiber. Remapping the known-sigma phi-ray truncation
set into an F statistic is anti-conservative and should not be used.

## Recommended Result Layout

Use these directories for current results:

| directory | role |
|---|---|
| `sims/results/sweep_cw_uv/` | headline CW-style sweep: known sigma, R-fiber unknown sigma, plug-ins, detection, conditional power |
| `sims/results/sweep_standard/` | extensions: K sweep, unbalanced clusters, AR covariance, sigma grid |
| `sims/results/*.rds` | focused experiments such as heavy tails, data thinning, real data, and validation |
| `sims/results/*.pdf`, `*.png` | figures generated from the above |

`sweep_cw_uv` is canonical for all paper figures because it includes the current
R-fiber columns (`p_union_rfib`, `p_path_rfib`) and plug-in columns. The older
`sweep_cw/` run has been archived to
`sims/results/_archive/sweep_cw_superseded_20260622.tar.gz` and removed from the
live tree.

## Core Scripts

Infrastructure:

- `checkpoint.R` - interrupt-safe checkpointing and resume.
- `sweep_harness.R` - one parameter cell to per-replicate rows and summaries.
- `sweep_run.R` - checkpointed driver and parameter grids.
- `house_style.R` - shared figure grammar and palette.

Method and validation scripts:

- `unknown_var_union.R` - valid unknown-variance path p-value from a known-sigma
  path set.
- `kmeans_union_unknownvar.R` - numerical R-fiber union/path test.
- `rfiber_exact.R` and `kmeans_union_unknownvar_exact.R` - exact R-fiber solver.
- `test_rfiber_exact.R`, `test_path_arc.R`, `test_f_region.R` - solver checks.
- `validate_Rfiber.R`, `validate_exact.R`, `validate_pivot_calibration.R` -
  calibration checks.

Main experiments:

- `power_type1_union_vs_path.R` - small known-variance Type I/power comparison.
- `sweep_run.R cw ...` - headline CW-style sweep.
- `sweep_run.R standard ...` - broader extension sweep.
- `exp_heavytail.R` - Gaussian/t5/t10 null robustness.
- `exp_datathin.R` - comparison to data thinning.
- `real_data_penguins.R`, `real_data_scrna.R` - real-data p-value tables.
- `union_multiplicity.R` - effective number of high-probability union regions.
- `repro_yunbarber.R` - matched known/unknown variance comparison.

Figure scripts:

- Headline sweeps: `plot_cw_typeI.R`, `plot_cw_power.R`, `plot_cw.R`,
  `plot_standard_sweep.R`, `plot_power_type1.R`.
- Unknown-variance geometry and comparisons: `plot_cw_path_known_vs_unknown.R`,
  `plot_repro_yunbarber.R`, `plot_repro_conditioning.R`,
  `plot_studentized_away.R`, `plot_decomposition.R`, `plot_circle_geometries.R`,
  `plot_theta_toy.R`, `plot_theta_matching.R`, `plot_theta_sphere.R`.
- Robustness/real-data/supporting figures: `plot_heavytail.R`,
  `plot_datathin.R`, `plot_gencov_typeI.R`, `plot_penguins.R`,
  `plot_scaling_overview.R`.

## Reproduce the Current Headline Set

One-time package setup:

```bash
R CMD INSTALL .
Rscript -e 'install.packages(c("intervals","ggplot2","gridExtra","datathin"))'
```

Fast correctness checks:

```bash
Rscript sims/test_rfiber_exact.R
Rscript sims/test_path_arc.R
Rscript -e 'library(testthat); library(KmeansInference); test_dir("tests/testthat")'
```

Headline CW-style sweep:

```bash
Rscript sims/sweep_run.R cw 7 50 sims/results/sweep_cw_uv
Rscript sims/sweep_run.R cw 7 50 sims/results/sweep_cw_uv finalize
Rscript sims/plot_cw_typeI.R  sims/results/sweep_cw_uv
Rscript sims/plot_cw_power.R  sims/results/sweep_cw_uv
Rscript sims/plot_cw.R        sims/results/sweep_cw_uv
Rscript sims/plot_cw_path_known_vs_unknown.R sims/results/sweep_cw_uv
```

Extension sweep:

```bash
Rscript sims/sweep_run.R standard 7 50 sims/results/sweep_standard
Rscript sims/sweep_run.R standard 7 50 sims/results/sweep_standard finalize
Rscript sims/plot_standard_sweep.R
Rscript sims/plot_gencov_typeI.R
```

Focused experiments:

```bash
Rscript sims/exp_heavytail.R
Rscript sims/plot_heavytail.R
Rscript sims/exp_datathin.R
Rscript sims/plot_datathin.R
Rscript sims/real_data_penguins.R
Rscript sims/plot_penguins.R
```

## Current Readiness Snapshot

The local ignored results currently include complete `sweep_cw_uv` and
`sweep_standard` summaries. The main story is coherent:

- Known-sigma null calibration is close to 0.05 across `q = 2, 10, 50, 100`.
- R-fiber unknown-variance null calibration is also close to 0.05.
- The known-sigma union test improves power over the path test in the main
  rising-power region, especially around moderate separations.
- The standard sweep covers K changes, unbalanced clusters, and AR covariance.
- Heavy-tail, data-thinning, penguins, and scRNA outputs exist as focused
  top-level results.

Remaining gaps before a polished paper run:

- Add or rerun the sigma `0.25` cells if a full Chen-Witten-style sigma panel is
  needed.
- Decide whether the exact R-fiber solver should replace the numerical R-fiber
  solver in the headline sweep, or remain a validation/accuracy result.
- Keep `sims/RESULTS_MANIFEST.md` updated after any long rerun so readers know
  which ignored outputs are current.

## Cluster Execution

The sweep is single-node shared-memory (`mclapply`) and checkpointed. Re-run the
same command to resume.

```bash
tmux new-session -d -s cw 'caffeinate -i Rscript sims/sweep_run.R cw 7 50 sims/results/sweep_cw_uv > sims/results/sweep_cw_uv.log 2>&1'
```

SLURM example:

```bash
#!/bin/bash
#SBATCH -c 16
#SBATCH --time=06:00:00
#SBATCH --requeue
module load R
cd "$SLURM_SUBMIT_DIR"
Rscript sims/sweep_run.R cw "$SLURM_CPUS_PER_TASK" 50 sims/results/sweep_cw_uv
```
