## Submission

This is an update of prCGGM 0.1.0 (on CRAN since 2026-11-17) to 0.2.0.

It is submitted soon after the first release because it brings
the package in line with the revised version of the paper it implements
(Mapelli et al. 2026, arXiv:2606.31805). The simulations of the revision
showed that the previous default selection rule (cross-validation with
`lambda.min`) selects increasingly dense networks as the sample size grows.
Version 0.2.0 therefore:

* adds the `selection` argument with relaxed BIC / extended BIC,
  BIC / extended BIC, cross-validation and stability-selection rules, and a
  least-squares refit of the selected edges (`refit`);
* changes the defaults to `selection = "bic_relaxed"`, `refit = TRUE` and
  `symm_method = "AND"`, the combination recommended in the paper;
* aligns `simulate_prcggm_data()` with the simulation design of the paper.

The changes are listed in `NEWS.md`. The settings of version 0.1.0 can be
reproduced with `selection = "cv_min", refit = FALSE, symm_method = "OR"`;
the deprecated argument `lambda_prec_type` still works, with a warning.

## Test environments

* local macOS (aarch64), R 4.3.0
* win-builder: Windows, R-devel
* win-builder: Windows, R-release
* R-hub v2 (GitHub Actions), R-devel (2026-11-17 to 2026-11-25,
  r90926–r90985): Ubuntu Linux (x86_64), Windows (x86_64), macOS (x86_64)
  and macOS (arm64); Status: OK on all four

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Days since last update: 8

  We apologise for the quick update. As explained above, the default
  selection rule of 0.1.0 gives overly dense networks for large samples,
  which is the main use case of the package (biobank-scale proteomics). We
  would like users to have the corrected defaults, which the revised paper
  describes, as soon as possible.

## Reverse dependencies

There are no reverse dependencies.

## Notes for the reviewers

* The examples of `prcggm_slurm_prepare()` that need a SLURM cluster are in
  `\dontrun{}`; the same workflow is run locally in the remaining examples,
  writing only to `tempdir()`.
* Parallel computation uses at most the number of cores requested through
  the `ncores` argument (default 1); examples and tests use at most 2.
