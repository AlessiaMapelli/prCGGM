test_that("all model variants return valid networks", {
  d <- small_data()
  set.seed(1)
  variants <- list(
    ggm = prcggm(d$x),
    weighted = prcggm(d$x, prior = d$prior),
    conditional = prcggm(d$x, covariates = d$covariates, tune = FALSE),
    full = prcggm(d$x, covariates = d$covariates, prior = d$prior,
                  tune = FALSE)
  )
  for (fit in variants) {
    expect_s3_class(fit, "prcggm")
    for (net in fit$networks) {
      expect_valid_network(net, 6)
      expect_equal(rownames(net), colnames(d$x))
    }
    expect_equal(nrow(fit$tuning), 6)
    expect_length(fit$sigma2, 6)
  }
  expect_named(variants$ggm$networks, "Baseline")
  expect_named(variants$full$networks, c("Baseline", "x1_1"))
  expect_equal(dim(variants$full$mean_effects), c(6, 2))
  expect_equal(diag(variants$full$precision$Baseline),
               unname(1 / variants$full$sigma2), ignore_attr = TRUE)
})

test_that("fits are reproducible with set.seed()", {
  d <- small_data()
  set.seed(42)
  a <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
              alpha_grid = c(0.5, 0.9), theta_grid = c(1, 1.3),
              random_search = TRUE)
  set.seed(42)
  b <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
              alpha_grid = c(0.5, 0.9), theta_grid = c(1, 1.3),
              random_search = TRUE)
  expect_identical(a$networks, b$networks)
  expect_true(all(a$tuning$alpha %in% c(0.5, 0.9)))
  expect_true(all(a$tuning$theta %in% c(1, 1.3)))
})

test_that("fixed hyperparameters use the first grid values", {
  d <- small_data()
  set.seed(1)
  fit <- prcggm(d$x, covariates = d$covariates, tune = FALSE,
                alpha_grid = c(0.9, 0.5), theta_grid = c(1.3, 1))
  expect_true(all(fit$tuning$alpha == 0.9))
  expect_true(all(fit$tuning$theta == 1.3))
})

test_that("the AND rule is not denser than the OR rule", {
  d <- small_data()
  set.seed(1)
  or <- prcggm(d$x, prior = d$prior, symm_method = "OR")
  set.seed(1)
  and <- prcggm(d$x, prior = d$prior, symm_method = "AND")
  expect_lte(sum(and$networks$Baseline != 0), sum(or$networks$Baseline != 0))
})

test_that("screening removes predictors but keeps prior-supported pairs", {
  d <- small_data()
  none <- matrix(0, 6, 6)
  set.seed(1)
  fit <- prcggm(d$x, screening_matrix = none)
  expect_true(all(fit$networks$Baseline == 0))
  set.seed(1)
  fit_prior <- prcggm(d$x, prior = d$prior, screening_matrix = none)
  outside <- d$prior == 0
  expect_true(all(fit_prior$networks$Baseline[outside] == 0))
  set.seed(1)
  fit_cor <- prcggm(d$x, screening = TRUE)
  expect_true(fit_cor$settings$screening)
  expect_true(is.numeric(fit_cor$settings$gamma))
})

test_that("a fixed lambda can be supplied", {
  d <- small_data()
  set.seed(1)
  fit <- prcggm(d$x, covariates = d$covariates, tune = FALSE,
                lambda_prec = 0.01, lambda_prec_type = "1se")
  expect_true(all(abs(fit$tuning$lambda - 0.01) <= 0.0001 + 1e-12))
})

test_that("invalid inputs are rejected", {
  d <- small_data()
  expect_error(prcggm(d$x[, 1:2]), "at least 3")
  x_na <- d$x
  x_na[1, 1] <- NA
  expect_error(prcggm(x_na), "missing values")
  expect_error(prcggm(d$x, prior = d$prior[1:5, 1:5]), "6 x 6")
  bad_names <- d$prior
  rownames(bad_names) <- paste0("G", 1:6)
  expect_error(prcggm(d$x, prior = bad_names), "Row names")
  expect_error(prcggm(d$x, covariates = d$covariates[1:10, , drop = FALSE]),
               "one row per row")
  expect_error(prcggm(d$x, alpha_grid = 2), "alpha_grid")
  expect_error(prcggm(d$x, symm_method = "XOR"))
  set.seed(1)
  expect_warning(prcggm(d$x, prior = -d$prior), "negative")
  set.seed(1)
  expect_warning(prcggm(d$x, prior = 2 * d$prior), "greater than 1")
})

test_that("parallel and sequential fits agree", {
  skip_on_cran()
  # workers load the installed package, so skip under devtools::load_all()
  skip_if(exists(".__DEVTOOLS__", envir = asNamespace("prCGGM")))
  d <- small_data()
  set.seed(3)
  seq_fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                    tune = FALSE)
  set.seed(3)
  par_fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                    tune = FALSE, ncores = 2)
  expect_equal(seq_fit$networks, par_fit$networks)
})
