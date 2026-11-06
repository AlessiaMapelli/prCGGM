# prCGGM 0.2.0 (development version)

## New features

* `prcggm()` and `prcggm_slurm_prepare()` gain a `selection` argument with
  eight rules for the penalty of the node-wise network regressions:
  `"ebic_relaxed"` (new default) and `"bic_relaxed"` (extended BIC / BIC of
  the least-squares refit of the predictors selected along the lambda path,
  jointly with `alpha` and `theta`, no cross-validation), `"bic"` and
  `"ebic"` (the same criteria on the penalised fit), `"cv_min"` and
  `"cv_1se"` (cross-validation, as in version 0.1.0) and `"stability_bic"`
  and `"stability_cv"` (stability selection at the hyperparameters chosen by
  the relaxed BIC or by cross-validation). Their settings are set with
  `selection_control`. The relaxed criteria are computed from the Gram
  matrix of each node, at a small extra cost.
* New argument `refit` (default `TRUE`): the predictors selected for each
  node are refitted by least squares, which removes the shrinkage of the
  edge weights without changing the selected edges. The penalised
  coefficients are kept in the new element `coefficients_penalised`.
* `tuning` gains the number of selected predictors `df`; `print()` and
  `summary()` report the selection rule.

## Breaking changes

* The default selection rule changes from cross-validation with
  `lambda.min` without refit to the relaxed extended BIC with refit, so
  default fits are sparser. In simulations, `lambda.min`, and BIC computed
  on the penalised fit, selected increasingly dense networks as the sample
  size grew. Use `selection = "cv_min", refit = FALSE` to reproduce
  version 0.1.0.
* The default hyperparameter grids are now `alpha_grid = c(0.75, 0.5, 0.9,
  0.99)` and `theta_grid = c(1, 0.8, 1.1, 1.3)` (were `c(0.75, 0.5, 0.9,
  0.95)` and `c(1, 0.8, 1.1, 1.3, 1.5)`), the grid evaluated in the
  simulations of the paper: 16 combinations per node instead of 20. With
  `tune = FALSE` the first values (0.75 and 1) are used, as before.
* `lambda_prec_type` is deprecated in favour of `selection`; `"min"` and
  `"1se"` are mapped to `"cv_min"` and `"cv_1se"` with a warning.
* Without covariates and without prior, the network regressions are now
  fitted with 'sparsegl' without intercept (as in all other cases) instead of
  `glmnet::cv.glmnet()` with intercept.
* `simulate_prcggm_data()` follows the simulation design of the paper:
  * the baseline network has a power-law degree distribution (fitness
    model, new argument `exponent`; `power` is deprecated) and `p` edges;
  * covariate networks have edge probability `1 / (p - 1)` by default and
    their own weight range `effect_range`;
  * a common diagonal makes the precision matrix positive definite for every
    subject (no subject-specific correction), so the model holds exactly;
  * continuous covariates are Uniform(0, 1) instead of standard normal;
  * the true baseline, and the prior, refer to the network at the mean of
    the continuous covariates, and the true effects of continuous covariates
    are given per standard deviation, the scale of the estimates;
  * the common diagonal is returned as `diagonal`.

# prCGGM 0.1.0

* First CRAN release. The package implements the method of Mapelli et al.
  (2026), arXiv:2606.31805, and replaces the `GGReg_*` scripts of the paper
  repository:
  * `prcggm()` replaces `GGReg_full_estimation()` and returns an S3 object
    with `print()`, `summary()`, `plot()` and `predict()` methods;
    `predict()` replaces `predict_personalized_network()`.
  * Sequential and parallel (`ncores`) fits keep all results in memory and
    write no files; the SLURM workflow is split into the non-interactive
    `prcggm_slurm_prepare()`, `prcggm_node()` and `prcggm_slurm_collect()`.
  * Results are reproducible with `set.seed()`; the package no longer sets
    the random seed itself.
