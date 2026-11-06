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
                lambda_prec = 0.01)
  expect_true(all(abs(fit$tuning$lambda - 0.01) <= 0.0001 + 1e-12))
})

test_that("relaxed EBIC with refit is the default and the refit keeps the support", {
  d <- small_data()
  set.seed(1)
  fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                tune = FALSE)
  expect_equal(fit$settings$selection, "ebic_relaxed")
  expect_true(fit$settings$refit)
  expect_equal(fit$coefficients != 0, fit$coefficients_penalised != 0)
  set.seed(1)
  pen <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                tune = FALSE, refit = FALSE)
  expect_equal(pen$coefficients, fit$coefficients_penalised)
  expect_equal(pen$networks$Baseline != 0, fit$networks$Baseline != 0)
  # least squares on the selected predictors: residuals orthogonal to them
  z <- fit$residuals
  i <- which.max(rowSums(fit$coefficients != 0))
  b <- fit$coefficients[i, -1]
  X <- cbind(z, as.numeric(as.character(d$covariates$x1)) * z)
  keep <- b != 0
  res <- z[, i] - drop(X %*% b)
  expect_lt(max(abs(crossprod(X[, keep, drop = FALSE], res))), 1e-6)
})

test_that("all selection rules run", {
  d <- small_data()
  rules <- c("ebic_relaxed", "bic_relaxed", "bic", "ebic", "cv_min", "cv_1se",
             "stability_bic", "stability_cv")
  for (rule in rules) {
    set.seed(1)
    fit <- prcggm(d$x, covariates = d$covariates, prior = d$prior,
                  alpha_grid = c(0.5, 0.9), theta_grid = 1,
                  selection = rule,
                  selection_control = list(stab_B = 10))
    expect_equal(fit$settings$selection, rule)
    for (net in fit$networks) expect_valid_network(net, 6)
    if (grepl("^stability", rule)) {
      expect_null(fit$coefficients_penalised)
      expect_true(all(is.na(fit$tuning$lambda)))
    } else {
      expect_false(anyNA(fit$tuning$lambda))
    }
  }
  # EBIC is never less sparse than BIC on the same path
  for (type in c("", "_relaxed")) {
    set.seed(1)
    bic <- prcggm(d$x, covariates = d$covariates, tune = FALSE,
                  selection = paste0("bic", type))
    set.seed(1)
    ebic <- prcggm(d$x, covariates = d$covariates, tune = FALSE,
                   selection = paste0("ebic", type))
    expect_true(all(ebic$tuning$df <= bic$tuning$df))
  }
})

test_that("the relaxed RSS equals the least-squares RSS of each support", {
  set.seed(2)
  X <- matrix(rnorm(200 * 8), 200, 8)
  y <- drop(X[, 1:3] %*% c(1, -0.5, 0.3)) + rnorm(200)
  S <- cbind(rep(FALSE, 8), c(TRUE, rep(FALSE, 7)), rep(c(TRUE, FALSE), 4),
             rep(TRUE, 8), c(TRUE, rep(FALSE, 7)))
  gram <- list(X = X, y = y, XtX = crossprod(X), Xty = drop(crossprod(X, y)),
               yty = sum(y^2))
  ols <- apply(S, 2, function(s) {
    if (!any(s)) sum(y^2) else sum(lm.fit(X[, s, drop = FALSE], y)$residuals^2)
  })
  expect_equal(relaxed_rss(S, gram), ols, tolerance = 1e-10)
  # rank-deficient support: least-squares fallback
  Xd <- cbind(X, X[, 1])
  gd <- list(X = Xd, y = y, XtX = crossprod(Xd), Xty = drop(crossprod(Xd, y)),
             yty = sum(y^2))
  s <- c(TRUE, rep(FALSE, 7), TRUE)
  expect_equal(relaxed_rss(matrix(s), gd),
               sum(lm.fit(Xd[, s], y)$residuals^2), tolerance = 1e-10)
})

test_that("stability selection is reproducible", {
  d <- small_data()
  set.seed(5)
  a <- prcggm(d$x, prior = d$prior, selection = "stability_bic",
              selection_control = list(stab_B = 10))
  set.seed(5)
  b <- prcggm(d$x, prior = d$prior, selection = "stability_bic",
              selection_control = list(stab_B = 10))
  expect_identical(a$networks, b$networks)
})

test_that("lambda_prec_type is deprecated but still works", {
  d <- small_data()
  set.seed(1)
  expect_warning(fit <- prcggm(d$x, covariates = d$covariates, tune = FALSE,
                               lambda_prec_type = "1se"), "deprecated")
  expect_equal(fit$settings$selection, "cv_1se")
})

test_that("invalid selection settings are rejected", {
  d <- small_data()
  expect_error(prcggm(d$x, selection = "aic"))
  expect_error(prcggm(d$x, selection_control = list(stab_b = 10)),
               "Unknown")
  expect_error(prcggm(d$x, selection_control = list(stab_pi = 0.3)),
               "stab_pi")
  expect_error(prcggm(d$x, refit = NA), "refit")
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
