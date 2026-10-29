test_that("the SLURM workflow reproduces prcggm()", {
  d <- small_data()
  dir <- tempfile("prcggm_slurm")
  set.seed(5)
  expect_message(
    job <- prcggm_slurm_prepare(d$x, covariates = d$covariates,
                                prior = d$prior, tune = FALSE, dir = dir),
    "sbatch --array=1-6")
  expect_true(file.exists(job$input))
  expect_true(dir.exists(file.path(dir, "logs")))
  expect_true(all(file.exists(file.path(dir, c("prcggm_node.sbatch",
                                               "prcggm_node_job.R")))))
  expect_error(prcggm_slurm_collect(dir), "missing for 6 node")
  for (i in seq_len(job$p)) {
    prcggm_node(job$input, i)
  }
  fit_slurm <- prcggm_slurm_collect(dir)
  set.seed(5)
  fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                tune = FALSE)
  expect_equal(fit_slurm$networks, fit$networks)
  expect_error(prcggm_node(job$input, 7))
  expect_error(prcggm_slurm_prepare(d$x), "dir")
  unlink(dir, recursive = TRUE)
})
