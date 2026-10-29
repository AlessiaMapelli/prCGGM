small_data <- function() {
  d <- prcggm_example
  idx <- seq_len(300)
  list(x = d$x[idx, 1:6],
       covariates = d$covariates[idx, , drop = FALSE],
       prior = d$prior[1:6, 1:6])
}

expect_valid_network <- function(net, p) {
  expect_true(is.matrix(net))
  expect_equal(dim(net), c(p, p))
  expect_equal(net, t(net))
  expect_true(all(diag(net) == 0))
}
