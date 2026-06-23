# Figure and Experiment Audit

This audit reflects the current tracked scripts and the local ignored summaries
documented in `sims/RESULTS_MANIFEST.md`.

## Current Figure Set

There are 19 tracked `plot_*.R` scripts.

Headline figures:

- `plot_cw_typeI.R` - null QQ curves for known sigma and R-fiber unknown sigma.
- `plot_cw_power.R` - known-sigma union/path power, detection probability, and
  conditional power.
- `plot_cw.R` - compact CW-style panels from the current `sweep_cw_uv` summary.
- `plot_power_type1.R` - fast small-sweep Type I and power sanity figure.
- `plot_standard_sweep.R` - K, balance, and covariance extension summary.

Unknown-variance and geometry figures:

- `plot_cw_path_known_vs_unknown.R`
- `plot_repro_yunbarber.R`
- `plot_repro_conditioning.R`
- `plot_studentized_away.R`
- `plot_decomposition.R`
- `plot_circle_geometries.R`
- `plot_theta_toy.R`
- `plot_theta_matching.R`
- `plot_theta_sphere.R`

Robustness, real-data, and supporting figures:

- `plot_heavytail.R`
- `plot_datathin.R`
- `plot_gencov_typeI.R`
- `plot_penguins.R`
- `plot_scaling_overview.R`

## Visual Grammar Status

`house_style.R` now defines the shared grammar:

- Colour encodes treatment or variance handling:
  `naive`, `oracle`, `studentized`, `med`, `sample`.
- Linetype encodes conditioning:
  `path` is solid, `union` is dashed.
- Facets encode swept dimensions such as `q`, `sigma`, or geometry.
- Grey `#4D4D4D` is reserved for reference lines.
- Shapes are avoided except for point-only dot plots or real-data panels where
  linetype cannot carry the distinction.

The major old conflicts have been resolved in the current plotting scripts:

- Known/oracle sigma is black, not grey.
- `q` is no longer encoded with the same blue used for the proposed method in
  the headline CW plots.
- The current CW plot defaults point to `sims/results/sweep_cw_uv`.

## Experiment Coverage

Ready for paper figures:

- Type I across `q = 2, 10, 50, 100` for known sigma and R-fiber unknown sigma.
- Power across `q = 2, 10, 50` and `delta = 1,...,8` for known sigma.
- Detection probability and conditional power for the same CW-style power cells.
- K sweep (`K = 2, 3, 4, 5` via `sweep_standard`), unbalanced clusters, and AR
  covariance checks.
- Heavy-tailed Type I robustness via `exp_heavytail.R`.
- Data-thinning comparison via `exp_datathin.R`.
- Penguins and scRNA real-data examples.

Partial or needs a decision:

- Sigma `0.25` is not in the current CW headline sweep; only `sigma = 1` appears
  in `sweep_cw_uv`, and `sigma = 0.5, 1` appear in `sweep_standard`.
- The headline unknown-variance sweep uses the numerical R-fiber solver. The
  exact solver is validated, but exact-vs-numerical values are not yet
  backfilled into the headline summary.
- `sweep_cw/` (the older superseded run) has been archived to
  `sims/results/_archive/sweep_cw_superseded_20260622.tar.gz` and removed from
  the live tree. `sweep_cw_uv/` is canonical for all paper figures.

## Recommended Paper Story

1. Start with Type I calibration (`plot_cw_typeI.R`): naive over-rejects,
   selective tests calibrate.
2. Show known-sigma power (`plot_cw_power.R`): union gains over path in the
   moderate-signal region, with detection and conditional power separated.
3. Explain why unknown variance changes the geometry (`plot_decomposition.R`,
   `plot_circle_geometries.R`, `plot_studentized_away.R`).
4. Use `plot_heavytail.R`, `plot_standard_sweep.R`, and `plot_gencov_typeI.R`
   as robustness/extension figures.
5. Use `plot_datathin.R` plus `plot_penguins.R` or `real_data_scrna.R` for
   comparison and application.
