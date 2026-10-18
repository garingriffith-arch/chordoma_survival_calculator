# Independent local rerun: 2024 NCDB chordoma analysis

**Status:** Completed locally on 18 October 2026. Manuscript release and repository-history review remain pending.

The ordered pipeline in `R/run_all.R` was executed against the locally licensed 2024 NCDB Bone/Joint Participant User File (R 4.6.0, Windows). This audit compares the locally exported, aggregate results with the previously locked 100-mm analysis. The NCDB source file, patient-level cohorts, model objects, and resampling checkpoints are not shared here.

## Execution

All six stages completed. The primary model used 1,243 patients and 310 deaths; the 946-patient common cohort had 242 deaths. Bootstrap optimism correction completed 300 of 300 iterations. The uncertainty analysis completed 120 of 120 outer iterations, with 30 requested inner iterations each. Multiple-imputation sensitivity and aggregate report generation completed.

The initial reporting run stopped on extraction of the reverse-Kaplan–Meier median. This was corrected in `R/06_tables_figures.R` (commit `db1b243`); reporting then completed without error.

## Comparison with the locked results

| Check | Result |
|---|---|
| Cohort sizes and event counts | Exact agreement |
| All 24 fitted Cox coefficients, three models | Exact agreement |
| Primary reference survival at 60 and 120 months | Exact agreement (0.880862 and 0.762530) |
| Primary-model equation vs fitted Cox prediction | Maximum absolute error: 4.44e-16 |
| All 19 pooled multiple-imputation coefficients | Exact agreement |
| Tumor-size sensitivity hazard ratios and nonlinear tests | Exact agreement |
| Apparent performance estimates | Exact agreement |
| Optimism-corrected performance | Largest absolute difference from earlier bootstrap run: 0.002643 |
| Nested-bootstrap percentile bounds | Largest absolute bound difference: 0.011587 |
| Reverse-KM median follow-up | Exact agreement: 80.92 and 82.60 months (primary/common) |

Independent resampling uses a different random-number stream than the earlier locked analysis. The corrected-performance and uncertainty differences above are consistent with Monte Carlo variation; they are not bitwise identical reproductions. All listed findings are based on the aggregate review package exported from the local run.

### Primary model: earlier and independent rerun

| Corrected measure | Earlier run | Independent rerun |
|---|---:|---:|
| Harrell C-index | 0.743661 | 0.742899 |
| 5-year AUC | 0.761793 | 0.760143 |
| 10-year AUC | 0.760265 | 0.759789 |
| 5-year Brier score | 0.129094 | 0.128481 |
| 10-year Brier score | 0.176958 | 0.176675 |

The manuscript should use one designated analysis run consistently rather than mix its point estimates with uncertainty intervals from another.

## Release boundaries

- This is an **internal reproducibility audit**, not independent validation on an external patient population.
- Stage 07 now **runs to completion locally** and reports Tables 1–4, S1–S9 and Figures 1–2, S1–S2 in CSV/DOCX and PDF/PNG. The user supplied a complete file inventory and run manifest dated 18 October 2026, but the contents of those generated exhibits have **not yet been independently audited**. Table S9 records unresolved histology/metastasis field adjudication rather than effect estimates.
- Aggregate files labeled `INTERNAL_DISCLOSURE_REVIEW`, and size-threshold tables permitting small-cell recovery by subtraction, must not be published without NCDB disclosure review.
- Earlier public Git history includes fitted `.rds` artifacts. Removing them from the branch working tree does not erase them from history. Inspect and remediate history before treating this repository as a fully clean public release.
- The linked calculator remains investigational. Tumor site C41.0 does not confirm intracranial location, and registry tumor size is not assured to be a pretreatment measurement.

## 18 October 2026: final exhibit code added

`R/07_publication_exhibits.R` passed syntax checks and was executed locally against the saved, locked 2024 PUF analysis on 18 October 2026 (R 4.6.0). The run finished and the output inventory verified the expected file names and formats. Content-level comparison with the manuscript and NCDB disclosure review are pending. It reads the completed private fitted models and local aggregate outputs without rerunning the 300-iteration or 120-outer bootstrap. The canonical performance series is the completed GitHub rerun shown above; older independently seeded estimates should not be intermixed. Publication exhibits remain ignored by Git under `outputs/private/manuscript_exhibits_DISCLOSURE_REVIEW` pending a local end-to-end and disclosure review.

## Numeric and figure audit of the generated exhibit set

The user provided a disclosure-limited review ZIP containing nine CSV tables and three PNG figures from Stage 07. The primary Cox coefficients and reference design agree with the rerun coefficient file within CSV rounding (maximum coefficient discrepancy below 5e-10). All 15 internal performance estimates agree with the canonical GitHub aggregate outputs at three-decimal publication precision. The reverse-KM medians are 80.92 and 82.60 months; the multiple-imputation coefficients and the earlier fixed-100-mm Riley-style sample-size scenario outputs agree.

A review of presentation identified three revisions to the *generator*: add baseline categorical characteristics with conservative whole-block suppression when nested cohorts may disclose small cells; collapse redundant five-/ten-year sample-size scenario rows when the estimated required N is identical; and add attrition arrows and clearer facet/variable labels to figures and supplementary coefficients. These improvements were committed after the audited output ZIP was generated. Their newly rendered output remains to be checked locally.

The three reviewed figures render and have appropriate internal-validation labels. The original cohort-flow figure depicts boxes without connecting arrows; a new connected generator was committed. The uploaded ZIP does not include Tables S1, S4, S5, S8, S9, or Figure S1; no claim is made that these excluded exhibits received a content-level disclosure audit. In particular, tumour-size threshold counts must not be publicized until small-cell and subtraction review is complete.

**Publication status remains conditional:** reconcile the final regenerated exhibits with the manuscript, complete NCDB disclosure review, resolve histology/metastasis reporting and the historical fitted-object exposure, then select one definitive release commit for the TRIPOD+AI code link.
