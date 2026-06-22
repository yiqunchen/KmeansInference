# Methods-paper draft: a more powerful selective test for k-means clustering

This directory holds a structured LaTeX skeleton for a statistics methods paper
introducing the *union* (condition-on-less) and *studentized R-fiber*
(unknown-variance) selective inference tests implemented in the
`KmeansInference` R package. The narrative and section structure mirror
Chen, Jewell & Witten (2023), *More powerful selective inference for the graph
fused lasso* (JCGS) — the "condition on less ⇒ more power" template.

- `paper.tex` — self-contained `article`-class draft (abstract, introduction,
  background, the union test, the unknown-variance R-fiber, computation,
  simulations, real data, discussion, figures). Every empirical number is pulled
  from the repository: `power_type1_summary.rds`, `repro_yunbarber_q{2,10}.rds`,
  `datathin_compare.rds`, `penguins_pvalues_k3.csv`, `scrna_pvalues_k3.csv`, and
  the two memory files. Figures are `\includegraphics` of the real PNGs in
  `../sims/results/` (via `\graphicspath`).
- `refs.bib` — all references, verified by web search (June 2026).

## Compiling

```sh
pdflatex paper
bibtex   paper
pdflatex paper
pdflatex paper
```

(LaTeX was **not** run here, and the figure PNGs are referenced by relative path
to `../sims/results/`; compile from inside `writeup/`.)

## Open TODOs (search the source for `% TODO:` / `\todo`)

- **Author list & affiliations** (`\author`).
- **Formal validity proposition** for the union test — state and prove that the
  conditioning event is φ-measurable so the truncated-χ pivot survives (mirror
  Prop. 1 of Chen, Jewell & Witten 2023). `paper.tex` §4.2.
- **H0′ vs H0 discussion** — spell out the pointwise-null requirement of the
  studentized test and the conservative/anti-conservative directions of a
  violation (source: `union-method-status.md`). §5.1.
- **Unit-test counts** — confirm 23/23 (F-region solver, `sims/test_f_region.R`),
  14/14 (7 primitive R-fiber modules, `sims/test_rfiber_exact.R`), 18/18
  (`path_arc` assembler, `sims/test_path_arc.R`) before stating them. §6.
- **K = 2..5 union-gain numbers** (0.075/0.115/0.225/0.200) come from the
  now-dropped `varying_k` figure folded into `standard_sweep`; re-confirm against
  `sims/results/standard_sweep` source data. §5.2.
- **Foreground (or not) the "do not remap the φ-union" caveat** — methodologically
  the most error-prone point; decide prominence. §5.
- **Future work** bullets (valid more-powerful union under unknown variance;
  multiple-pairs FWER per Yun & He 2024; randomized/data-fission variants).

Numbers that are **not** marked TODO are real and taken from the sources above;
do not treat them as placeholders.
