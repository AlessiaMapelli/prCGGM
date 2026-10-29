## Submission

This is a new submission.

## Test environments

* local macOS (aarch64), R 4.3.0
* win-builder: Windows Server 2022, R-devel (2026-10-28 r90783 ucrt)
* win-builder: Windows Server 2022, R-release (4.6.2)
* macOS builder (R-release)            <!-- TODO -->
* R-hub (Linux, R-devel)               <!-- TODO -->

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Possibly misspelled words in DESCRIPTION:
    Mapelli, al, et, proteomics

  These are the author's surname in the reference to the method
  ("Mapelli et al. (2026)") and a standard technical term.

## Notes for the reviewers

* The examples of `prcggm_slurm_prepare()` that need a SLURM cluster are in
  `\dontrun{}`; the same workflow is run locally in the remaining examples,
  writing only to `tempdir()`.
* Parallel computation uses at most the number of cores requested through
  the `ncores` argument (default 1); examples and tests use at most 2.
