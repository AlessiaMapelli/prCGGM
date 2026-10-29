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
