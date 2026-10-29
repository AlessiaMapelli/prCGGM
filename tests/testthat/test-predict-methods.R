test_that("predict combines baseline and covariate components", {
  d <- small_data()
  set.seed(1)
  fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                tune = FALSE)
  nets <- predict(fit, data.frame(x1 = factor(c("0", "1"))))
  expect_named(nets, c("Subject_1", "Subject_2"))
  expect_equal(nets$Subject_1, fit$networks$Baseline)
  expect_equal(nets$Subject_2, fit$networks$Baseline + fit$networks$x1_1)
  expect_equal(predict(fit), list(Baseline = fit$networks$Baseline))
  prec <- predict(fit, data.frame(x1 = "1"), type = "precision")
  expect_equal(prec$Subject_1,
               fit$precision$Baseline + fit$precision$x1_1)
  expect_error(predict(fit, data.frame(x1 = "2")), "not seen")
  expect_error(predict(fit, data.frame(z = 1)), "Missing covariate")
})

test_that("numeric covariates are scaled with the training parameters", {
  set.seed(2)
  sim <- simulate_prcggm_data(n = 300, p = 5, n_binary = 0,
                              n_continuous = 1)
  sim$covariates$x1 <- 50 + 10 * sim$covariates$x1
  fit <- prcggm(sim$x, covariates = sim$covariates, tune = FALSE)
  m <- mean(sim$covariates$x1)
  s <- sd(sim$covariates$x1)
  nets <- predict(fit, data.frame(x1 = c(m, m + s)))
  expect_equal(nets$Subject_1, fit$networks$Baseline)
  expect_equal(nets$Subject_2, fit$networks$Baseline + fit$networks$x1)
})

test_that("print, summary and plot methods work", {
  d <- small_data()
  set.seed(1)
  fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                tune = FALSE)
  expect_output(print(fit), "Features: 6")
  s <- summary(fit)
  expect_s3_class(s, "summary.prcggm")
  expect_equal(s$components$component, c("Baseline", "x1_1"))
  expect_output(print(s), "Network components")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, component = "x1_1", threshold = 0.05))
  expect_invisible(plot(fit, newdata = data.frame(x1 = "1")))
  expect_error(plot(fit, component = "nope"), "must be one of")
})

test_that("network utilities", {
  a <- matrix(c(0, 0.4, 0.4, 0), 2, 2)
  b <- matrix(c(0, 0.1, 0.1, 0), 2, 2)
  expect_equal(compare_networks(a, b), mean(abs(a - b)))
  expect_equal(compare_networks(a, b, "frobenius"), norm(a - b, "F"))
  expect_error(compare_networks(a, diag(3)), "same dimensions")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot_network(a, threshold = 0.5)
  expect_true(all(plotted == 0))
  expect_equal(unname(plot_network_difference(a, b)), a - b)
})

test_that("simulate_prcggm_data returns consistent objects", {
  set.seed(1)
  sim <- simulate_prcggm_data(n = 50, p = 6, n_binary = 1, n_continuous = 1,
                              prior = "noisy")
  expect_equal(dim(sim$x), c(50, 6))
  expect_named(sim$truth_precision, c("Baseline", "x1_1", "x2"))
  expect_named(sim$truth_networks, c("Baseline", "x1_1", "x2"))
  for (h in names(sim$truth_networks)) {
    net <- sim$truth_networks[[h]]
    prec <- sim$truth_precision[[h]]
    expect_equal(net, t(net))
    expect_true(all(diag(net) == 0))
    off <- row(net) != col(net)
    expect_true(all(sign(net[off]) == -sign(prec[off])))
  }
  expect_true(is.factor(sim$covariates$x1))
  expect_true(all(sim$prior >= 0 & sim$prior <= 1))
  expect_equal(sim$prior, t(sim$prior))
})
