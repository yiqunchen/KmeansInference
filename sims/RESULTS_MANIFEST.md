# Results Manifest

This file records how to interpret the generated files in `sims/results/`.
The result files themselves are git-ignored and should be regenerated from the
scripts in `sims/`.

Snapshot audited locally: 2026-06-22.

## Headline Sweep: `sims/results/sweep_cw_uv/`

Status: complete local summary.

- Logical cells: 28.
- Materialized chunks: 176.
- Type I cells: `cw_typeI_q{2,10,50,100}`, each with 1,000 valid replicates.
- Power cells: `cw_power_q{2,10,50}_d{1,...,8}`, each with 200 valid replicates.
- Columns include known-sigma union/path, naive, valid R-fiber unknown-sigma
  union/path, MED/sample plug-ins, detection probability, and conditional power.

Null rejection rates at alpha = 0.05:

| label | known union | known path | R-fiber union | R-fiber path | naive |
|---|---:|---:|---:|---:|---:|
| `cw_typeI_q2` | 0.061 | 0.049 | 0.047 | 0.048 | 1.000 |
| `cw_typeI_q10` | 0.063 | 0.058 | 0.046 | 0.044 | 1.000 |
| `cw_typeI_q50` | 0.032 | 0.037 | 0.029 | 0.027 | 0.994 |
| `cw_typeI_q100` | 0.063 | 0.061 | 0.054 | 0.054 | 0.986 |

Largest known-sigma union gains in the current power cells:

| label | union | path | gain | detection | conditional gain |
|---|---:|---:|---:|---:|---:|
| `cw_power_q2_d4` | 0.795 | 0.400 | 0.395 | 0.990 | 0.394 |
| `cw_power_q10_d4` | 0.670 | 0.300 | 0.370 | 0.900 | 0.394 |
| `cw_power_q2_d3` | 0.495 | 0.215 | 0.280 | 0.905 | 0.293 |
| `cw_power_q50_d5` | 0.670 | 0.405 | 0.265 | 0.755 | 0.331 |
| `cw_power_q10_d5` | 0.910 | 0.655 | 0.255 | 0.930 | 0.269 |

Use this directory for:

- `plot_cw_typeI.R`
- `plot_cw_power.R`
- `plot_cw.R`
- `plot_cw_path_known_vs_unknown.R`

## Extension Sweep: `sims/results/sweep_standard/`

Status: complete local summary.

- Logical cells: 36.
- Materialized chunks: 168.
- Includes Type I by `q`, power by `q` and `sigma`, K sweep, AR covariance, and
  unbalanced clusters.

Selected extension cells:

| label | union | path | gain | note |
|---|---:|---:|---:|---|
| `typeI_q2` | 0.058 | 0.053 | 0.005 | null calibration |
| `typeI_q10` | 0.050 | 0.045 | 0.005 | null calibration |
| `typeI_q50` | 0.038 | 0.040 | -0.003 | null calibration |
| `k2` | 0.955 | 0.880 | 0.075 | K sweep |
| `k4` | 0.470 | 0.245 | 0.225 | K sweep |
| `k5` | 0.355 | 0.155 | 0.200 | K sweep |
| `unbal_d5` | 0.875 | 0.575 | 0.300 | unbalanced sizes |
| `gencov_ar0.5_d0` | 0.068 | 0.045 | 0.023 | AR covariance null |
| `gencov_ar0.9_d0` | 0.068 | 0.035 | 0.033 | AR covariance null |

Use this directory for:

- `plot_standard_sweep.R`
- `plot_gencov_typeI.R`

## Focused Top-Level Outputs

These are generated directly under `sims/results/`:

| output stem | generator | plotter | purpose |
|---|---|---|---|
| `power_type1_*` | `power_type1_union_vs_path.R` | `plot_power_type1.R` | fast known-sigma smoke comparison |
| `heavytail_typeI` | `exp_heavytail.R` | `plot_heavytail.R` | Gaussian/t5/t10 Type I robustness |
| `datathin_compare` | `exp_datathin.R` | `plot_datathin.R` | data-thinning comparison |
| `penguins_pvalues_k3` | `real_data_penguins.R` | `plot_penguins.R` | penguins real-data example |
| `scrna_pvalues_k*` | `real_data_scrna.R` | none separate beyond `scrna_k*` | scRNA real-data example |
| `repro_yunbarber_q*` | `repro_yunbarber.R` | `plot_repro_yunbarber.R`, `plot_repro_conditioning.R` | matched known/unknown comparison |
| `union_multiplicity` | `union_multiplicity.R` | none | effective region-count diagnostic |

## Readiness Assessment

Ready:

- Package-level union and path tests pass.
- Exact R-fiber primitive tests pass.
- Headline known-sigma Type I, R-fiber Type I, power, detection, and conditional
  power summaries are materialized.
- Extension summaries for K, balance, and AR covariance are materialized.
- Heavy-tail, data-thinning, penguins, and scRNA focused outputs exist locally.

Resolved 2026-06-22:

- EXACT R-fiber backfill done. `sweep_cw_uv/` is now the EXACT arc-sweep run
  (`kmeans_union_unknownvar_exact`, the harness default; set
  `KM_RFIBER_SOLVER=numerical` for the grid solver). The prior numerical run is
  archived at `_archive/sweep_cw_uv_numerical_20260622.tar.gz`. Determinism
  verified: known-variance + plug-in columns reproduce the numerical run exactly
  (max|diff| = 0); only the `*_rfib` columns moved (rejection rates by up to
  ~0.026). Validation: `sims/validate_exact_backfill.R`.
- sigma `0.25` cells: SKIP (power depends on delta/sigma, so it only relabels the
  x-axis; adds no conclusion).
- The superseded `sweep_cw/` run was archived to
  `_archive/sweep_cw_superseded_20260622.tar.gz` and removed; all figures use
  `sweep_cw_uv/`.

Paper figures (writeup/paper.tex): the Type-I and power TABLES were replaced by a
single combined boxplot `cw_box_combined.png` (delta=0 calibration | delta>0
power; `sims/plot_cw_box_combined.R`); the data-thinning table by the supplementary
`datathin_sigma_box.png` (`sims/plot_datathin_sigma_box.R`).
