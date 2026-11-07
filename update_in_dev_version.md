# Changes in the development version of prCGGM

This file tracks every change made since version 0.1.0:

- **Version 0.1.0:** branch `main`, commit `c8393c8`, the CRAN submission.
- **Where the changes are:** branch `v.0.2` (committed 2026-11-06), version `0.1.0.9000`.
- **Release plan:** to be released as 0.2.0 once CRAN has accepted 0.1.0. The user-facing summary is in `NEWS.md`.
- **Build:** this file is excluded from the build (`.Rbuildignore`).
- **Companion file:** the paper repository has its own tracking file, `Complete_projects/Prior-informed-conditional-GGMs/update_in_dev_version.md`. The package and the template `ggReg_main_functions.R` implement the same estimator and are checked against each other.

Last updated: 2026-11-07. The changes of 2026-11-07 (3b) are committed on `v.0.2`.

---

## Summary

| # | Area | Files | Breaking? |
|---|---|---|---|
| 1 | Selection rules for λ and the hyperparameters | `R/fit.R`, `R/node_estimation.R`, `R/slurm.R` | default changes |
| 2 | OLS refit of the selected predictors | `R/fit.R`, `R/node_estimation.R` | default changes |
| 3 | Relaxed BIC / EBIC | `R/node_estimation.R`, `R/fit.R`, `R/slurm.R` | default changes |
| 3b | Final defaults `selection = "bic_relaxed"`, `symm_method = "AND"` | `R/fit.R`, `R/slurm.R`, tests, docs | default changes |
| 4 | Lasso without covariates fitted by `sparsegl` | `R/node_estimation.R` | yes, slightly |
| 5 | `print()` / `summary()` report the selection rule | `R/methods.R` | no |
| 6 | `simulate_prcggm_data()` follows the paper's design | `R/simulate.R` | yes |
| 6b | Default hyperparameter grids aligned with the paper | `R/fit.R`, `R/slurm.R` | default changes |
| 7 | Documentation | `man/`, `vignettes/prCGGM.Rmd`, `NEWS.md`, `DESCRIPTION`, `inst/WORDLIST` | — |
| 8 | Tests | `tests/testthat/` | — |

---

## 1. Selection rules (2026-11-04, extended 2026-11-06)

- **New argument `selection`** of `prcggm()` and `prcggm_slurm_prepare()`:

  | Rule | λ | α, ϑ | Extra cost |
  |---|---|---|---|
  | `"bic_relaxed"` (default, see 3b) | min BIC of the least-squares refit along the path | jointly with λ | one Gram matrix per node |
  | `"ebic_relaxed"` | min EBIC (γ = 0.5) of the least-squares refit | jointly | as above |
  | `"bic"`, `"ebic"` | min BIC / EBIC of the penalised fit | jointly | none |
  | `"cv_min"`, `"cv_1se"` | `nfolds`-fold CV, `lambda.min` / `lambda.1se` (rule of 0.1.0) | BIC at the selected λ | CV |
  | `"stability_bic"` | stability selection | those of `"bic_relaxed"` | 50 subsample paths |
  | `"stability_cv"` | stability selection | those of `"cv_min"` | CV + 50 subsample paths |

- **New argument `selection_control`:**
  - settings: `ebic_gamma` 0.5, `stab_B` 50, `stab_fraction` 0.5, `stab_pi` 0.75, `stab_pfer` 1;
  - validation is done by `check_selection_control()`;
  - unknown names are an error.
- **Stability selection** (Meinshausen & Bühlmann, 2010):
  - the λ path of the chosen (α, ϑ) is refitted on `stab_B` subsamples of size `stab_fraction`·n;
  - the path is restricted to the region where the average number of selected predictors is ≤ √(PFER(2π − 1)P);
  - predictors with selection frequency ≥ π somewhere in that region are kept.

  The subsamples are drawn in `prepare_fit()` before fitting, so they are shared by all nodes and reproducible with `set.seed()`, also in parallel and SLURM runs.
- **`lambda_prec_type` deprecated:**
  - it is now the last argument;
  - `"min"`/`"1se"` are mapped to `"cv_min"`/`"cv_1se"` with a warning (`resolve_selection()`);
  - the mapping applies only when `selection` is not given.
- **Fixed `lambda_prec`:** the combination is chosen by BIC at that λ.
- **`tuning`** gains `df`, the number of selected predictors. Its `lambda` is `NA` for stability selection.

## 2. OLS refit (2026-11-04)

- New argument `refit` (default `TRUE`): the selected predictors of each node are refitted by least squares (`refit_ols()`).
  - This removes the shrinkage of the edge weights; the support is unchanged.
  - σ²_j = RSS/(n − df).
  - The refit is always applied with stability selection.
- New element `coefficients_penalised` of the fit: the penalised coefficients of the selected model. It is `NULL` for stability selection.

## 3. Relaxed BIC / EBIC and new default (2026-11-06)

- **Problem found in the paper simulations.** With BIC/EBIC computed on the RSS of the penalised fit:
  - the lasso shrinkage inflates the RSS by an amount that grows with n;
  - so the selected λ falls like 1/n;
  - so false positives increase with n (e.g. baseline F1 0.95 → 0.87 from n = 1k to 20k, perfect prior, p = 50).
- **New internal function `relaxed_rss(S, gram)`:**
  - computes the RSS of the least-squares refit of each distinct support (column of S) of the path;
  - RSS = y'y − b'X_s'y with b = (X_s'X_s)⁻¹X_s'y, from the Gram matrix X'X computed once per node (`fit_node()`);
  - QR fallback when X_s'X_s is singular;
  - identical to `lm.fit` to 1e-14.
- **`fit_node()`** computes `bic_relaxed` and `ebic_relaxed` along each path when the rule needs them.
- **Default `selection = "ebic_relaxed"`** on 2026-11-06 (was `"bic"` in the development version, `cv_min` without refit in 0.1.0); replaced by `"bic_relaxed"` on 2026-11-07 (see 3b).
- **Evidence for the choice** (paper repository, `Simulation_analysis_v2/selection_test/` and `Add_simulation_analysis/`):
  - the relaxed criteria remove the decline with n;
  - EBIC relaxed is best or tied in 23 of 24 prior × OR/AND × n cells (overall F1);
  - it matches stability selection at about one fifth of the cost;
  - BIC relaxed over-selects when P is large relative to n.
- **`stability_bic`** now takes its hyperparameters from `"bic_relaxed"`, which adds under 2% to the cost.

## 3b. Final defaults: relaxed BIC with AND (2026-11-07)

- `selection` default `"bic_relaxed"` (was `"ebic_relaxed"`) and `symm_method` default `"AND"` (was `"OR"`, also in 0.1.0), in `prcggm()` and `prcggm_slurm_prepare()`; order of `rules` in `resolve_selection()` updated.
- **Reason** (paper simulations, Supplementary Parts 1 and 2, 100 replications): the relaxed BIC with AND is the most accurate combination for n ≥ 5000 and keeps improving with n, in particular for binary covariate effects (the setting of the paper's case study). The relaxed EBIC with OR is best for small samples (n ≈ 1000) and is recommended there. F1, n = 1k / 5k / 10k / 20k (p = 100, q = 3, noisy prior):

  | Component | BIC relaxed + AND | EBIC relaxed + OR |
  |---|---|---|
  | Baseline | 0.90 / 0.98 / 0.99 / 0.99 | 0.91 / 0.97 / 0.97 / 0.96 |
  | Binary effect | 0.63 / 0.97 / 0.99 / 0.99 | 0.64 / 0.95 / 0.96 / 0.95 |
  | Continuous effect | 0.08 / 0.53 / 0.83 / 0.97 | 0.06 / 0.58 / 0.86 / 0.97 |
  | Overall | 0.71 / 0.89 / 0.95 / 0.98 | 0.73 / 0.89 / 0.94 / 0.96 |

- Documentation: "Selection rules" section of `prcggm()` (BIC relaxed first, with the pairing guideline), `symm_method` argument, vignette ("Hyperparameters": table and paragraph on the default and the guideline), `NEWS.md` (new breaking-change item on `symm_method`; 0.1.0 reproduced with `selection = "cv_min", refit = FALSE, symm_method = "OR"`), `DESCRIPTION`.
- Test: the default fit has `selection = "bic_relaxed"` and `symm_method = "AND"`.

## 4. Lasso without covariates

- Without covariates, the network regressions are now fitted with `sparsegl` (one group per predictor, no intercept), as in all other cases.
- Before, they used `glmnet::cv.glmnet()` with intercept.

## 5. Methods

- `print.prcggm()` and `summary.prcggm()` report the selection rule and the refit. `summary()$fit` gains `selection`.
- **New internal `selection_label()`:** fits made with 0.1.0 (no `selection` setting) are labelled `cv_min`/`cv_1se`, penalised coefficients.

## 6. `simulate_prcggm_data()` (paper design)

- **Baseline network:**
  - power-law degree distribution, fitness model, new argument `exponent` (default 2.5);
  - `p` edges;
  - `power` is deprecated and ignored (it gave super-linear preferential attachment, i.e. a star).
- **Covariate networks:**
  - edge probability `edge_prob`, default 1/(p − 1);
  - weight range `effect_range`.
- **Common diagonal:**
  - set from the minimum eigenvalue over the corners of the covariate box (`eig_margin`), so the precision matrix is positive definite for every subject;
  - no subject-specific correction, so the model holds exactly;
  - handles 0 covariates;
  - returned as `diagonal`.
- **Covariates:** continuous ones are Uniform(0, 1), no longer standard normal; binary ones are factors.
- **Truth:**
  - the true `Baseline`, and the simulated prior, refer to the network at the mean of the continuous covariates and the reference level of the binary ones;
  - true effects of continuous covariates are per standard deviation, the scale of the estimates.

## 6b. Default hyperparameter grids (2026-11-06)

- `alpha_grid = c(0.75, 0.5, 0.9, 0.99)` and `theta_grid = c(1, 0.8, 1.1, 1.3)` in `prcggm()` and `prcggm_slurm_prepare()`.
  - They were `c(0.75, 0.5, 0.9, 0.95)` and `c(1, 0.8, 1.1, 1.3, 1.5)`.
  - These are the values evaluated in the paper's simulations: 16 combinations per node instead of 20.
- α = 0.99 was the value chosen most often in the simulations; ϑ = 1.5 was never evaluated.
- The first values (0.75, 1) are unchanged, so `tune = FALSE` fits are unaffected.
- `random_search` stays `FALSE` (deterministic full grid). The paper's simulations use 50% random search for cost reasons.

## 7. Documentation

- `prcggm()`:
  - new section "Selection rules";
  - new arguments documented;
  - references Chen & Chen (2008) and Meinshausen & Bühlmann (2010).
- `prcggm_slurm_prepare()`: new arguments (inherited).
- Vignette: "Hyperparameters" section with the rules table, the relaxed criteria and the default.
- `NEWS.md`: 0.2.0 (development) entry with new features and breaking changes.
- `DESCRIPTION`: version 0.1.0.9000; the description names the default rule.
- `inst/WORDLIST`: new terms.
- `man/` regenerated with roxygen2.

## 8. Tests

`tests/testthat/test-prcggm.R`:
- relaxed BIC, AND and refit are the defaults, and the refit keeps the support;
- the refit is least squares (residuals orthogonal to the selected predictors);
- all eight rules run;
- EBIC is never denser than BIC, on the penalised and on the relaxed RSS;
- `relaxed_rss()` equals `lm.fit` for each support, including the rank-deficient fallback;
- stability selection is reproducible;
- `lambda_prec_type` is deprecated but works;
- invalid `selection` / `selection_control` / `refit` are rejected.

`tests/testthat/test-predict-methods.R`: the new simulator:
- common diagonal;
- prior built from the baseline target;
- Uniform(0, 1) continuous covariates;
- deprecated `power`;
- `prior = "none"`.

The selection label in `print()`/`summary()` is exercised by the existing print/summary tests but not checked explicitly.

---

## Verification (2026-11-06)

| Check | Result |
|---|---|
| `devtools::test()` | all pass (1 skip: parallel test under `load_all`) |
| `R CMD check --as-cran --no-manual` (local, RStudio Pandoc) | 0 errors, 0 warnings, 2 NOTEs (new submission / maintainer; unable to verify current time) |
| Package vs template, BIC/EBIC with and without refit | identical (2026-11-04) |
| Package vs template, `ebic_relaxed`, `bic_relaxed`, `bic`, `stability_bic` (`prcggm_example`, `mean_estimation = FALSE`) | same λ, α, ϑ and network supports |

## Open points

- **Release:** branch `v.0.2`. Bump to 0.2.0 and merge into `main` after CRAN accepts 0.1.0.
- **Default grids:** aligned with the paper (6b, decided 2026-11-06).
- **Mean step:** still a cross-validated lasso (`lambda_mean_type`). Unchanged for now.
