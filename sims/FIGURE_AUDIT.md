Confirmed: noise is strictly Gaussian (`rnorm`), so heavy-tails (t5/t10) have zero infrastructure. All findings are grounded. I'll now write the report.

# K-means Selective-Inference Figure Set — Audit Report

This report audits our figure set (`sims/`) against Yun & Barber (2023) and Chen & Witten (2023, JMLR), grounded in the actual repo state: `sims/sweep_harness.R`, `sims/sweep_run.R` (`build_grid`), the materialized `sims/results/sweep_cw_uv/summary.csv` and `sims/results/sweep_standard/summary.csv`, `sims/house_style.R`, the 14 `plot_*.R` scripts, and `vignettes/real_data_example.Rmd`.

Two facts up front that shape everything below:
- The harness **draws noise only via `rnorm`** (`sweep_harness.R:204` `Z <- matrix(rnorm(...))`, then `%*% SigChol` or `* sd_sim`). There is **no t-distribution / heavy-tail mechanism anywhere** — so t5/t10 robustness is not "partial," it is structurally absent.
- The harness **does** compute detection probability + conditional power (`detect_prob`, `cpow_*`) and the two plug-ins (`p_union_med/p_path_med`, `p_union_samp/p_path_samp`), and these are materialized in `sweep_cw_uv/summary.csv`. So several "power" items are further along than the figure set reveals.

---

## 1. GAP ANALYSIS

Legend: **Have** = run + figure exists; **Partial** = infra/data exist but no figure, or only a thin slice run; **Missing** = neither. Priority: **P1** = must-have for the paper; **P2** = nice; **P3** = optional.

| # | Analysis | Who has it | Why it matters | Our status | Priority |
|---|----------|-----------|----------------|-----------|----------|
| 1 | **Type I across q** (global-null QQ, q∈{2,10,50,100}) | C&W Fig 4; Y&B Fig 3 | Core validity claim; reviewers check calibration first | **Have** — `sweep_cw_uv` has `cw_typeI_q{2,10,50,100}` @1000 reps; `plot_cw_typeI.R` → `cw_typeI.{pdf,png}`, faceted `variance ~ q` (known + unknown). Strong. | P1 |
| 2 | **Power vs delta** (separated triangle) | C&W Fig 5; Y&B Fig 4 | Headline benefit; our union gain lives here | **Have** — `cw_power_q{2,10,50}_d{1..8}` @200; `plot_cw_power.R` → `cw_power.{pdf,png}` (union vs path, `variance ~ q`). | P1 |
| 3 | **Power vs delta swept over sigma** (σ∈{0.25,0.5,1}) | C&W Fig 5/8/9 | C&W's main multi-σ panel; shows σ̂_Sample bias grows with delta | **Partial** — `standard` grid has σ∈{0.5,1} only; `cw` grid is σ=1 only; **σ=0.25 never run**. `plot_cw.R` panel B draws σ but on partial data (tagged "partial-data"). Add σ=0.25 cells and a clean multi-σ figure. | P1 |
| 4 | **Detection probability vs delta** (eq. 24) | C&W Fig 5-left, 9 | Separates "did k-means find true clusters" from "did the test reject" — required to interpret conditional power honestly | **Partial** — `detect_prob`/`n_detected` computed and in `sweep_cw_uv`; rendered as `cw_detection.{pdf,png}` left panel (`plot_cw_power.R` p_det). Exists but only at σ=1, one geometry. | P1 |
| 5 | **Conditional power vs delta** (eq. 23) | C&W Fig 5-right, 9 | The estimand under which union/path/plug-ins are honestly comparable | **Partial** — `cpow_union/path/_unk` computed; `cw_detection` right panel (p_cp). Same σ=1 / single-geometry limitation as #4; no σ stratification. | P1 |
| 6 | **Heavy-tailed robustness — t5** (global-null Type I) | Y&B Fig 5 | Standard misspecification stress test; reviewers expect it | **Missing** — no `rt`/`rmvt` path in harness; noise is pure Gaussian. Needs a `noise_dist` knob. | P1 |
| 7 | **Heavy-tailed robustness — t10** | Y&B Fig 6 | Milder-tail companion to #6 | **Missing** — same as #6. | P2 |
| 8 | **Non-isotropic / unequal-variance Type I** (diag(1,2)) | Y&B Fig 7 | Tests calibration when isotropy fails — directly relevant since our unknown-σ test assumes isotropic scale | **Partial** — `Sigma` knob exists; `standard` grid runs `gencov_ar0.5/0.9_d0` (AR(1) global-null Type I), which is a *stronger* anisotropy test than Y&B's diagonal. But there is **no figure** rendering the genCov Type-I cells, and no diagonal-heteroscedastic cell to match Y&B exactly. Add a genCov Type-I QQ figure. | P1 |
| 9 | **Equidistant vs collinear/horizontal geometry** | Y&B Fig 4B vs 4C | Shows power conclusions aren't an artifact of one geometry | **Missing** — only `geometry="triangle"` (equilateral) is run. `single-blob`/`elongated` exist but **no collinear geometry** and no horizontal/triangle contrast figure. Add a `collinear` geometry + side-by-side. | P2 |
| 10 | **Spherical vs elongated/anisotropic geometry** (power) | local extension (neither paper) | Our `elongated` knob is novel infra; demonstrates union robustness to anisotropic separation | **Partial** — `geometry="elongated"` (`elong=3`) is wired in `sweep_harness.R:168` but **no grid cell uses it** and no figure. Cheap to add. | P2 |
| 11 | **Cluster balance (unbalanced sizes)** | local extension | Reviewers ask "what if clusters differ in size" | **Partial** — `unbal_d0`/`unbal_d5` (balance c(1,2,3)) run in `standard`; surfaced only inside `plot_standard_sweep.R` panel-2 dumbbell. No dedicated balance figure. | P2 |
| 12 | **Varying K** (K∈{2,4,5}) | local extension (C&W/Y&B fix K=3/K=2) | Our exact R-fiber solver works for any K — this is a *selling point* with no paper analog | **Partial** — `k2`/`k4`/`k5` run in `standard`, shown only in `plot_standard_sweep.R` dumbbell. Deserves a standalone "union gain vs K" figure. | P1 (novelty) |
| 13 | **Varying n** | C&W use n=150/100/30; we fix n=90 | Sensitivity to sample size; cheap robustness | **Missing** — single n=90 (cw) / n=60 (standard). No n-sweep cell. | P3 |
| 14 | **Sigma misspecification** (sig_test ≠ sig, or wrong SigInv) | local extension | Probes the unknown-σ test under a deliberately wrong σ | **Partial** — `sig_test`/`SigInv_test` knobs fully wired (`sweep_harness.R:51-60`), but **every run cell uses the correct value** (`sig_test==sig` throughout both summaries). Infra ready, never exercised, no figure. | P3 |
| 15 | **Variance-decomposition / motivation figure** | Y&B Fig 2 (X=P0X+P1X+P2X); C&W Fig 1-3 | The conceptual figure reviewers use to *understand* the construction; our R-fiber `x'(r)` perturbation is exactly Y&B Fig 3's spirit and is our novel coordinate | **Partial** — `plot_studentized_away.R` visualizes the truncation-set mass (good, novel), but there is **no perturbation/decomposition schematic** showing the φ-fiber vs R-fiber and the P0/P1/P2 split. This is our cleanest "why union+R-fiber" pedagogical figure and it's absent. | P1 |
| 16 | **Real-data example + p-value table** | C&W Figs 6/7/10+Table 1 (penguins/MNIST/scRNA); Y&B Fig 8+Table 1 (penguins) | **Every** competing paper has ≥1 real-data application with a p-value table; a methods paper without one is a guaranteed reviewer complaint | **Missing (effectively)** — `vignettes/real_data_example.Rmd` is the *inherited* C&W scRNA-seq vignette and calls only the old `kmeans_inference` (path) + naive σ̂_MED. It does **not** use our union test, has no `union` column, produces no figure/table into `sims/results/`. We have **zero** real-data output demonstrating the union/R-fiber contribution. | **P1** |
| 17 | **QQ that naive over-rejects + selective fixes it** (motivation) | C&W Fig 1; Y&B Fig 1 | The "why selective inference at all" opener | **Have** — `plot_power_type1.R` panel-2 ECDF (union/path) + `cw_typeI` include naive over-rejection visually. Adequate. | P2 |
| 18 | **Unconditional power vs effect size ‖μᵀν‖** (spline) | C&W Fig 8 | Smooth power-vs-signal summary | **Missing** — we plot vs delta only. Low value given #2/#5. | P3 |

### Where our novel work has NO analog in either paper (reviewers will want these standalone — feature them)

- **Union vs path power gain under KNOWN σ** (`power_gain`, `cpow_gain` columns; `plot_standard_sweep.R` panel-1, `plot_cw_power.R`). This is *the* novel result. Currently scattered across figures; deserves one crisp headline figure. **Have, but under-featured.**
- **"Studentized away"**: union gain shrinks under the studentized-F vs known-σ (`plot_studentized_away.R`, `cw_typeI_conditioning.R`, `repro_conditioning.R`). Genuinely novel; no paper analog. **Have.**
- **Fixed-conditioning-set ordering oracle ≥ studentized-F ≥ plug-ins** (`plot_repro_yunbarber.R`; `or_p/uf_p/md_p/sa_p`). This is our reframing of Y&B; **Have** but only q∈{2,10}.
- **Exact R-fiber studentized-F solver for K>2** (`rfiber_exact.R`, `kmeans_union_unknownvar_exact.R`) replacing Y&B's importance sampling. **Validated** (`validate_exact.R`) but has **no figure** quantifying its advantage (accuracy vs IS, or runtime). A "exact vs importance-sampling" calibration/cost figure would be a strong standalone — currently **Missing**.

---

## 2. AESTHETIC AUDIT

Every cross-figure channel collision, with figures and hex.

### COLOUR — blue `#1F5AA6` carries **five** meanings
1. **union method** — `plot_cw_typeI.R`, `plot_cw_power.R` (Fig 2 power), `plot_power_type1.R`, `plot_standard_sweep.R` panel 2, `plot_cw.R` panel C. *(canonical, via `km_pal[["union"]]`)*
2. **studentized / unknown-σ data series** — `plot_cw_path_known_vs_unknown.R`, `plot_cw_typeI_conditioning.R`, `plot_repro_yunbarber.R` (`studentized-F`), `plot_repro_conditioning.R` (`unknown (R-fiber F)`).
3. **q = 2 level** — `plot_cw_power.R` (p_det/p_cp), `plot_cw.R` panel A (`qcol`), `plot_standard_sweep.R` panel 1.
4. **σ = 0.25 level** — `plot_cw.R` panel B (`scol`).
5. **"union gain (below obs)" fill region** — `plot_studentized_away.R`.

Worst single-figure offenders: `plot_standard_sweep.R` (blue = q=2 in panel 1, union in panel 2 of the **same** figure) and `plot_cw.R` (blue = q=2 in A, σ=0.25 in B, union in C).

### COLOUR — grey `#4D4D4D` (`km_ref`) overloaded as reference line AND data
- Correctly the 45°/null/zero reference in `plot_cw_typeI.R`, `plot_power_type1.R`, `plot_standard_sweep.R`, `plot_cw.R`.
- But a **data series** ("known/oracle σ") in `plot_cw_path_known_vs_unknown.R`, `plot_cw_typeI_conditioning.R`, `plot_repro_yunbarber.R`, `plot_repro_conditioning.R`.
- **Worst case:** `plot_cw_typeI_conditioning.R` uses `#4D4D4D` *both* as the "known" data line *and* as the 45° `geom_abline` reference in the same panel → visually indistinguishable.

### COLOUR — rose `#C48A97` (`km_pal[["path"]]`)
- "path method" in canonical figures; but "**path set** region-fill" in `plot_studentized_away.R`.

### COLOUR — q=10 has two near-identical oranges
- `km_pal[["naive"]] = #E69F00` (orange) but every q-palette renders **q=10 = `#C98A2E`** (`plot_cw_power.R`, `plot_cw.R` `qcol`, `plot_standard_sweep.R`). Two oranges, different concepts.

### COLOUR — house-style-orphan hexes introduced ad hoc
- `plot_cw.R`: `#B55D4C` (red, σ=1), `#9A4D8E` (purple, q=100), `#2B7A78` (teal, q=50) — none in `house_style.R`.
- `plot_cw_power.R` p_det/p_cp: `#2B7A78` (teal, q=50).
- `plot_repro_yunbarber.R`: `#009E73` (green, σ_MED) — fine as a *new role* but undeclared in the palette.
- `plot_scaling_overview.R`: abandons manual palette entirely for `scale_colour_viridis_d(option="C")` to encode delta — a totally different colour system.
- `plot_studentized_away.R`: observed-stat vline hard-coded `#222222` (not `km_ref`); plus grey `#BBBBBB` for "union extra (above obs)".

### LINETYPE — union flips solid↔dashed
- union = **solid (1)** / path = dashed (2): `plot_cw_power.R` (Fig 2), `plot_repro_conditioning.R`.
- "weaker (union / R-fiber)" = **dashed (2)** / path = solid (1): `plot_cw_path_known_vs_unknown.R`, `plot_cw_typeI_conditioning.R`.
- `plot_cw.R` panel C: linetype encodes **variance** (known=1/unknown=2), not method at all.
- Reference-line linetype inconsistent: dashed (2) in most; `plot_cw.R` uses dotted (3) for the α hline in panel B while panel A uses dashed (2) — mixed within one figure. Reference alpha varies 0.55 / 0.6 / 0.65.

### SHAPE — values inconsistent + almost always redundant
- **path**: 17 (triangle) in `plot_cw_power.R`/`plot_cw.R`; **1 (open circle)** in `plot_repro_conditioning.R`; **16 (same as union → indistinguishable)** in `plot_power_type1.R`.
- **naive**: 4 (×) in `plot_power_type1.R`; never shaped elsewhere.
- plug-ins/oracle: 15/16/18/17 in `plot_repro_yunbarber.R`.
- Redundancy: shape duplicates an already-encoded channel in essentially every figure (redundant with colour=method in cw_power Fig 2 *and* with linetype; with colour=q in p_det / standard_sweep p1; with linetype=cond in repro_conditioning; with linetype in cw panels B/C).

### LABEL / NAMING — four names for the unknown-σ method
- "studentized (unknown)" (cw_path_known_vs_unknown, cw_typeI_conditioning) / "unknown (R-fiber F)" (repro_conditioning) / "studentized-F (unknown sigma)" (repro_yunbarber) / plain "unknown" (cw panel C). And "weaker (union / R-fiber)" vs "union" for the conditioning level. Canonical `km_lab` strings are used **only** in the `km_pal`-based figures.

### STYLE-RULE violations of `house_style.R` header
- Header says **NO subtitles** + **ONE shared bottom legend**, yet `plot_cw_path_known_vs_unknown.R`, `plot_cw_typeI_conditioning.R`, `plot_repro_yunbarber.R`, `plot_repro_conditioning.R` all add `subtitle=`.
- **Inset legends** instead of shared bottom: `plot_repro_yunbarber.R` `c(0.13,0.80)`, `plot_repro_conditioning.R` `c(0.12,0.78)`, `plot_studentized_away.R` `c(0.66,0.78)`, `plot_cw.R` panel A `c(0.02,0.98)`.

**Net:** only the canonical `km_pal` family (`plot_cw_typeI`, `plot_cw_power` Fig 2 power panel, `plot_power_type1`, `plot_standard_sweep` panel 2) is mutually consistent. The conditioning/variance family redefines blue+grey; the q-faceted family redefines blue/orange/teal/purple; `plot_studentized_away` reuses union-blue and path-rose as region fills.

---

## 3. UNIFIED VISUAL GRAMMAR

### The core tension and its resolution
Some figures compare **selection methods** (naive / path / union); others compare **variance handling** (known / studentized-F / σ_MED / σ_Sample). These are **two orthogonal axes** and must never share a channel. Decision:

> **COLOUR = the "treatment" being compared = variance handling.**
> **LINETYPE = conditioning level (path / union/weaker).**
> **FACET = the sweep dimension (q, σ, or geometry — whichever varies).**
> **SHAPE = dropped globally** (it is redundant with linetype or colour in every current figure; the user's rule is to drop redundant shape). Re-introduce shape *only* if a figure needs a third categorical that is not already on colour or linetype — none currently does.

Rationale: variance handling is the dimension with the most levels (4: known/studentized-F/MED/Sample) and the one whose *colour identity* readers must track across the whole paper, so it earns colour. Conditioning is binary-ish (path vs union, occasionally "weaker") and reads cleanly as solid-vs-dashed. q and σ are *experimental sweeps*, not treatments — they belong on facets/x-axis, never on colour. **naive** is a selection-method baseline, not a variance method; treat it as a fifth colour role used only in the motivation/Type-I figures where it appears.

### Canonical semantic palette (role → hex)

| Role | Hex | Channel | Notes |
|------|-----|---------|-------|
| **naive (invalid)** | `#E69F00` | colour | keep `km_pal[["naive"]]`; only in motivation/Type-I |
| **known / oracle σ** | `#000000` (black) — *or* keep `#4D4D4D` ONLY as a data colour and forbid grey for reference lines | colour | black makes "oracle" unmistakable and frees grey for references |
| **studentized-F (exact unknown σ)** | `#1F5AA6` (blue) | colour | the proposed unknown-σ method = blue, our headline colour |
| **σ_MED plug-in** | `#009E73` (green) | colour | already used in `repro_yunbarber`; promote to canon |
| **σ_all / sample plug-in** | `#D55E00` (vermillion) | colour | distinct from naive-orange; Okabe-Ito-safe |
| **path (single Lloyd path)** | — | **linetype solid (1)** | conditioning, never a colour |
| **union (proposed conditioning)** | — | **linetype dashed (2)** | union/weaker/R-fiber = dashed, fixed forever |
| **weaker (union / R-fiber)** | — | **linetype dashed (2)** | same as union (they are the same conditioning concept) |
| **reference lines (45°, null, zero, α)** | `#4D4D4D` (grey), `linetype=2`, `alpha=0.6` | annotation | grey is **reference-only** — never a data series |

Fixed global rules:
- **union/weaker = dashed (2); path = solid (1)** everywhere. (This flips `plot_cw_power.R` and `plot_repro_conditioning.R`, which currently make union solid. The conditioning family is already dashed-for-union, so adopt that convention — it also reads correctly as "the augmented set adds the dashed extension.")
- **grey `#4D4D4D` is forbidden as a data colour.** Oracle/known moves to **black `#000000`**. This single change resolves the worst collision (`cw_typeI_conditioning` grey-line-vs-grey-reference).
- **blue `#1F5AA6` means exactly one thing: the proposed unknown-σ (studentized-F) / and, in pure selection-method figures, the union method.** It must **never** encode q or σ. q and σ move to facets or x-axis.
- **q → facet columns; σ → facet rows or linetype-free small multiples; geometry → facet.** Never colour.
- **Shape scales deleted** from all scripts.
- **Labels** standardize to: `"Naive (invalid)"`, `"Known/oracle σ"`, `"Studentized-F (unknown σ)"`, `"σ_MED plug-in"`, `"σ_sample plug-in"`; conditioning: `"Path"`, `"Union"`. Retire "R-fiber F", "weaker", "studentized (unknown)" variants.
- **No subtitles; single shared bottom legend** via `km_panels()` (enforce the existing house rule; remove all inset `legend.position=c(...)` and `subtitle=`).

### How each existing figure re-maps

| Selection-method figures (colour = method, naive/path/union) | becomes |
|---|---|
| `plot_cw_typeI.R`, `plot_power_type1.R`, `plot_standard_sweep.R` p2, `plot_cw.R` C | colour=method (naive `#E69F00` / path: drop colour, use solid line / union: blue). Cleaner: keep these as **colour=method** since variance handling isn't varied here, but force union=dashed/path=solid so they're consistent with the variance figures. |

| Variance-handling figures (colour = variance, linetype = conditioning) | becomes |
|---|---|
| `plot_repro_yunbarber.R` | colour = {oracle black, studentized-F blue, σ_MED green, σ_sample vermillion}; **drop shape (15/16/18/17)**; facet = q. Already 90% there — just recolour σ_sample to `#D55E00`, oracle to black, delete shape, move legend to bottom, drop subtitle. |
| `plot_cw_path_known_vs_unknown.R` | colour = sig (known=black, studentized=blue); linetype = cond (path solid, union dashed). Already correct **except** known must go black, drop subtitle. |
| `plot_cw_typeI_conditioning.R` | same as above; **critical fix**: known → black so it's distinct from the grey 45° reference. |
| `plot_repro_conditioning.R` | colour = sig (known black / studentized blue); linetype = cond — **flip union to dashed**; delete shape scale (16/1). |

### Per-script REFACTOR CHECKLIST

- **`house_style.R`** — Add to `km_pal`/`km_lab` the canonical variance roles: `oracle="#000000"`, `studentized="#1F5AA6"`, `sig_med="#009E73"`, `sig_samp="#D55E00"`. Add a `km_cond` linetype map `c(path=1, union=2)`. Add `km_var_lab`/`km_cond_lab` strings. Document: grey `#4D4D4D` = reference-only; shape never used. Add a `scale_km_q_facet` helper so q is always a facet.
- **`plot_cw_typeI.R`** — keep (it's canonical). Drop any shape; confirm reference grey only.
- **`plot_cw_power.R`** (Fig 2 / p_det / p_cp) — **flip union to dashed(2), path solid(1)**; **delete `scale_shape_manual`**; move q off colour onto facet for p_det/p_cp (currently q=colour `#1F5AA6/#C98A2E/#2B7A78`) → facet by q, colour by method. Remove teal/orange q-hexes.
- **`plot_cw.R`** — biggest rebuild: panel A move q from colour (`qcol`) to **facet**; panel B move σ from colour (`scol` incl. `#B55D4C`) to **facet or linetype-free panels**; panel C set known=black/unknown=blue + union dashed. Remove `#9A4D8E`, `#2B7A78`, `#B55D4C`, `#C98A2E`. Unify the α-line linetype to 2 (not 3). Move panel-A inset legend to shared bottom.
- **`plot_standard_sweep.R`** — panel 1: move q off blue/`#C98A2E` onto facet or keep as a *q-only* figure with an explicit q-colour scale that is **never** reused for method; panel 2: keep method colours but force path solid / union dashed; drop shape; reconcile so blue isn't "q=2" and "union" in one figure (make panel-1 q a facet).
- **`plot_repro_yunbarber.R`** — recolour σ_sample → `#D55E00`, oracle → black; **delete `scale_shape_manual(15,16,18,17)`**; legend to bottom; remove subtitle.
- **`plot_repro_conditioning.R`** — **flip union to dashed(2)**; delete `scale_shape_manual(16,1)`; known→black; add error bars or note absence; legend bottom; remove subtitle.
- **`plot_cw_path_known_vs_unknown.R`** — known→black (was grey); confirm union=dashed; remove subtitle; legend bottom.
- **`plot_cw_typeI_conditioning.R`** — **known→black** (resolves grey-on-grey); union=dashed; remove subtitle; legend bottom.
- **`plot_power_type1.R`** — fix shape: currently union=16/path=16 (indistinguishable). **Delete shape**; rely on path-solid/union-dashed + colour. naive stays `#E69F00`.
- **`plot_studentized_away.R`** — fills currently reuse union-blue/path-rose as *regions*. Keep blue/rose for the path-set vs union-gain fills (this is the one figure where blue=union-gain is defensible **because there is no method/variance series to confuse it with**) but change the observed-stat vline from `#222222` to `km_ref`. Document the exception in a header comment.
- **`plot_scaling_overview.R`** — delta currently on viridis-C. This is a *cost* figure with no method/variance/conditioning series, so a sequential delta scale is fine; just switch to a single declared sequential palette in `house_style.R` (e.g. `scale_colour_viridis_d` is acceptable here) and note the exception. No grammar conflict because none of colour's reserved meanings appear.

### Where the unified grammar is impossible / needs an exception (flagged)
- **`plot_studentized_away.R`**: fill encodes *regions of a truncation set* (path set / union-gain / no-help), not methods. There is no method or variance series in this figure, so reusing union-blue for "union-gain region" cannot be confused. **Allowed exception**, documented in-script.
- **`plot_scaling_overview.R`**: colour = delta is a legitimate *sweep* encoding in a pure cost figure with no treatment/conditioning series. **Allowed exception** (sequential palette), documented; it is the only figure permitted a non-`km_pal` colour system.
- **`plot_cw.R` panel A (Type-I QQ across q)**: q genuinely needs to be visible and there are up to 4 q-levels and a method axis. Resolution = **facet by q**, colour by method — fully compatible, but it forces a layout change (4 facets) rather than overlaying q on colour. Not impossible, just a structural edit.

**Single highest-leverage change:** move **oracle/known from grey `#4D4D4D` to black `#000000`** and reserve grey strictly for reference lines. That one rule eliminates the worst collision (the indistinguishable grey-line-vs-grey-45° in `cw_typeI_conditioning`) and frees grey to mean "reference" unambiguously across all 14 scripts.

---

Relevant absolute paths:
`/Users/yiqun/Desktop/chen-lab/more-powerful-kmeans/sims/sweep_harness.R` (geometry/detection/cpow/plug-in/sig_test wiring; noise = `rnorm` only, line 204), `/Users/yiqun/Desktop/chen-lab/more-powerful-kmeans/sims/sweep_run.R` (`build_grid`: cw + standard grids; no t-tail/elong/sig_test≠sig/σ=0.25 cells), `/Users/yiqun/Desktop/chen-lab/more-powerful-kmeans/sims/results/sweep_cw_uv/summary.csv` (detection/cpow/MED/sample columns materialized), `/Users/yiqun/Desktop/chen-lab/more-powerful-kmeans/sims/results/sweep_standard/summary.csv` (k2/k4/k5, AR(1) genCov, unbalanced), `/Users/yiqun/Desktop/chen-lab/more-powerful-kmeans/sims/house_style.R` (palette to extend), `/Users/yiqun/Desktop/chen-lab/more-powerful-kmeans/vignettes/real_data_example.Rmd` (inherited scRNA vignette — path+naive only, no union, no figure/table output → the real-data gap).