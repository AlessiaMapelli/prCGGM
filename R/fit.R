#' Fit a prior-informed conditional Gaussian graphical model
#'
#' Estimates a Gaussian graphical model in which both the mean and the
#' precision matrix depend linearly on subject-level covariates, while prior
#' knowledge on the network (for example 'STRING' interaction confidence
#' scores) reduces the penalty on supported edges of the population-level
#' network.
#'
#' Throughout, \eqn{n} is the number of subjects, \eqn{p} the number of
#' features (network nodes) and \eqn{q} the number of covariate columns after
#' encoding (see *Covariate encoding*).
#'
#' Estimation follows two steps. First, each feature is regressed on the
#' covariates with a lasso penalty and the residuals are kept. Second, for each node
#' \eqn{j} the residuals are regressed on the residuals of the other nodes and
#' on their interactions with the covariates, using the sparse-group lasso
#' penalty
#' \deqn{\alpha\lambda \|(1 - W_{j,-j}) \odot d_j^0\|_1 +
#'   \vartheta\alpha\lambda \|d_j^{-0}\|_1 +
#'   (1 - \alpha)\lambda \sum_{h=1}^q \|d_j^h\|_2,}{alpha * lambda *
#'   ||(1 - W[j, -j]) * d0_j||_1 + theta * alpha * lambda * ||d-0_j||_1 +
#'   (1 - alpha) * lambda * sum_h ||dh_j||_2,}
#' where \eqn{d_j^0} are the population-level (baseline) coefficients,
#' \eqn{d_j^h} the coefficients modulated by covariate \eqn{h} and \eqn{W} the
#' prior weight matrix. The prior therefore only affects the baseline network.
#' The penalty \eqn{\lambda} is chosen by cross-validation, and the structural
#' parameters \eqn{\alpha}{alpha} and \eqn{\vartheta}{theta} are chosen separately for each
#' node by minimising the node-wise BIC
#' \eqn{n \log(RSS_j / n) + \log(n) df_j}{n * log(RSS_j / n) + log(n) * df_j}.
#' Node-wise estimates are finally combined into symmetric matrices with the
#' "OR" or "AND" rule.
#'
#' Depending on the inputs, the same function fits a standard
#' neighbourhood-selection GGM (no covariates, no prior), a prior-weighted GGM
#' (prior only), a conditional GGM (covariates only) or the full
#' prior-informed conditional GGM.
#'
#' **Covariate encoding.** Numeric covariates are standardised and give one
#' column each. Factors (and character or logical columns) with \eqn{C} levels
#' are dummy coded into \eqn{C - 1} columns with the first level as the
#' reference, so the baseline network refers to the reference group; use
#' [factor()] or [relevel()] to choose the reference level. The total number
#' of resulting columns is \eqn{q}.
#'
#' **Reproducibility.** Cross-validation folds (shared by all nodes and
#' hyperparameter values) and the random-search subsets are drawn from the
#' current random number generator state before any model is fitted, so
#' calling [set.seed()] beforehand makes the fit reproducible, also when
#' `ncores > 1`.
#'
#' @param x Numeric \eqn{n \times p}{n x p} matrix or data frame of expression data,
#'   with the \eqn{n} subjects in rows and the \eqn{p} features (network
#'   nodes) in columns. Missing values are not allowed.
#' @param covariates Optional data frame of subject-level covariates with
#'   \eqn{n} rows, in the same order as the rows of `x`. Columns can be
#'   numeric, factor, character or logical; they are encoded into \eqn{q}
#'   columns as described in *Covariate encoding*. If `NULL` (default), a
#'   single population network is estimated (\eqn{q = 0}).
#' @param prior Optional symmetric \eqn{p \times p}{p x p} matrix of prior edge
#'   confidence in \[0, 1\], with rows and columns in the same order as the
#'   columns of `x`. Edges with entries close to 1 are penalised less in the
#'   baseline network. If `NULL` (default), no prior knowledge is used.
#' @param screening Logical; if `TRUE`, predictors weakly correlated with the
#'   response node (absolute correlation of the residuals not above `gamma`)
#'   are excluded from each node-wise regression, unless they have prior
#'   support. Default `FALSE`.
#' @param gamma Correlation threshold in \[0, 1\] for screening. If `NULL`
#'   (default), the 20th percentile of the \eqn{p(p - 1)/2}{p(p - 1)/2} absolute pairwise
#'   correlations is used.
#' @param screening_matrix Optional \eqn{p \times p}{p x p} matrix with entries in
#'   \{0, 1\}; entry \eqn{(j, k) = 0} excludes feature \eqn{k} from the
#'   regression of node \eqn{j}. When supplied, it is used for screening
#'   instead of the correlation rule. Prior-supported pairs are never
#'   excluded.
#' @param mean_estimation Logical; if `TRUE` (default) the covariate effects
#'   on the mean are estimated and removed first. Set to `FALSE` only if `x`
#'   has already been centred conditionally on the covariates.
#' @param lambda_mean Optional fixed lasso penalty for the mean regressions.
#'   If `NULL` (default), it is selected by cross-validation with
#'   [glmnet::cv.glmnet()].
#' @param lambda_mean_type Cross-validation rule used to select
#'   `lambda_mean`, with the same meaning as in [glmnet::cv.glmnet()]:
#'   `"min"` (default) selects `lambda.min`, the value that minimises the
#'   cross-validated error; `"1se"` selects `lambda.1se`, the largest value
#'   whose error is within one standard error of the minimum, which gives
#'   sparser models. Ignored if `lambda_mean` is supplied or if \eqn{q = 1},
#'   in which case the mean is estimated by least squares.
#' @param lambda_prec Optional fixed penalty for the node-wise network
#'   regressions. If `NULL` (default), it is selected by cross-validation
#'   with [sparsegl::cv.sparsegl()] (or [glmnet::cv.glmnet()] when there are
#'   neither covariates nor prior).
#' @param lambda_prec_type Cross-validation rule used to select
#'   `lambda_prec`: `"min"` (default, `lambda.min`) or `"1se"`
#'   (`lambda.1se`), with the same meaning as for `lambda_mean_type`.
#'   Ignored if `lambda_prec` is supplied.
#' @param tune Logical; if `TRUE` (default) \eqn{\alpha}{alpha} and \eqn{\vartheta}{theta}
#'   are tuned by BIC over `alpha_grid` and `theta_grid`. If `FALSE`, the
#'   first value of each grid is used.
#' @param alpha_grid Candidate values in \[0, 1\] for \eqn{\alpha}{alpha}, the balance
#'   between element-wise sparsity (\eqn{\alpha \to 1}{alpha -> 1}) and group sparsity
#'   (\eqn{\alpha \to 0}{alpha -> 0}).
#' @param theta_grid Positive candidate values for \eqn{\vartheta}{theta}, the
#'   penalty on covariate-specific effects relative to the baseline penalty.
#' @param random_search Logical; if `TRUE`, each node evaluates a random
#'   subset of the \eqn{\alpha}{alpha} and \eqn{\vartheta}{theta} hyperparameter grid instead of the full grid. Default
#'   `FALSE`.
#' @param p_random Proportion of the grid evaluated when `random_search` is
#'   `TRUE`. Default 0.5.
#' @param nfolds Number of cross-validation folds. Default 5.
#' @param symm_method Symmetrisation rule: `"OR"` (default) keeps an edge if
#'   it is selected in either node-wise regression, `"AND"` only if it is
#'   selected in both.
#' @param ncores Number of cores used to fit the node-wise regressions in
#'   parallel with the \pkg{parallel} package. Default 1 (sequential). For
#'   very large networks see [prcggm_slurm_prepare()].
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return An object of class `"prcggm"`, a list with elements:
#' \describe{
#'   \item{`networks`}{Named list of symmetric \eqn{p \times p}{p x p} adjacency
#'     matrices of partial-correlation-type edge weights: `"Baseline"` (the
#'     population network) followed by one matrix per covariate effect,
#'     named after the encoded covariate (for example `"sex_Male"`).}
#'   \item{`precision`}{Named list with the corresponding precision matrix
#'     components \eqn{\hat\Delta_h}. The diagonal of `"Baseline"` holds the
#'     inverse residual variances.}
#'   \item{`mean_effects`}{\eqn{p \times (q + 1)}{p x (q + 1)} matrix of covariate effects
#'     on the mean (intercept first), or `NULL` if `mean_estimation = FALSE`.}
#'   \item{`residuals`}{Matrix of residuals from the mean step, used as input
#'     to the network step.}
#'   \item{`coefficients`}{Matrix of node-wise regression coefficients
#'     (one row per node) before symmetrisation.}
#'   \item{`sigma2`}{Residual variance of each node-wise regression.}
#'   \item{`tuning`}{Data frame with, for each node, the selected `lambda`,
#'     `alpha`, `theta`, the `bic` and the computing time in minutes.}
#'   \item{`covariates`}{Encoding information used by [predict.prcggm()].}
#'   \item{`settings`}{List of the main estimation settings.}
#'   \item{`n`, `p`, `q`}{Number of subjects, features and encoded
#'     covariates.}
#'   \item{`call`}{The matched call.}
#' }
#'
#' @references
#' Mapelli, A., Massi, M. C., Cuccuru, G., Di Angelantonio, E. and Ieva, F.
#' (2026). Prior-informed conditional Gaussian graphical models: an
#' application to protein interaction network reconstruction.
#' arXiv:2606.31805. \doi{10.48550/arXiv.2606.31805}
#'
#' Zhang, J. and Li, Y. (2023). High-dimensional Gaussian graphical
#' regression models with covariates. *Journal of the American Statistical
#' Association*, 118(543), 2088--2100.
#'
#' @seealso [predict.prcggm()], [summary.prcggm()], [plot.prcggm()],
#'   [prcggm_slurm_prepare()].
#'
#' @examples
#' data(prcggm_example)
#' set.seed(1)
#' fit <- prcggm(prcggm_example$x,
#'               covariates = prcggm_example$covariates,
#'               prior = prcggm_example$prior,
#'               tune = FALSE)
#' fit
#' summary(fit)
#'
#' \donttest{
#' # Tune alpha and theta by BIC on a small grid
#' fit_tuned <- prcggm(prcggm_example$x,
#'                     covariates = prcggm_example$covariates,
#'                     prior = prcggm_example$prior,
#'                     alpha_grid = c(0.5, 0.9), theta_grid = c(1, 1.3))
#' fit_tuned$tuning
#' }
#' @export
prcggm <- function(x,
                   covariates = NULL,
                   prior = NULL,
                   screening = FALSE,
                   gamma = NULL,
                   screening_matrix = NULL,
                   mean_estimation = TRUE,
                   lambda_mean = NULL,
                   lambda_mean_type = c("min", "1se"),
                   lambda_prec = NULL,
                   lambda_prec_type = c("min", "1se"),
                   tune = TRUE,
                   alpha_grid = c(0.75, 0.5, 0.9, 0.95),
                   theta_grid = c(1, 0.8, 1.1, 1.3, 1.5),
                   random_search = FALSE,
                   p_random = 0.5,
                   nfolds = 5,
                   symm_method = c("OR", "AND"),
                   ncores = 1L,
                   verbose = FALSE) {
  call <- match.call()
  ncores <- check_number(ncores, "ncores", lower = 1, null_ok = FALSE)
  prep <- prepare_fit(
    x = x, covariates = covariates, prior = prior, screening = screening,
    gamma = gamma, screening_matrix = screening_matrix,
    mean_estimation = mean_estimation, lambda_mean = lambda_mean,
    lambda_mean_type = match.arg(lambda_mean_type),
    lambda_prec = lambda_prec,
    lambda_prec_type = match.arg(lambda_prec_type),
    tune = tune, alpha_grid = alpha_grid, theta_grid = theta_grid,
    random_search = random_search, p_random = p_random, nfolds = nfolds,
    symm_method = match.arg(symm_method), verbose = verbose)
  results <- run_nodes(prep, ncores = as.integer(ncores), verbose = verbose)
  fit <- assemble_fit(prep, results)
  fit$call <- call
  fit
}

# Validates the inputs, runs the mean step and stores everything the
# node-wise regressions need.
prepare_fit <- function(x, covariates, prior, screening, gamma,
                        screening_matrix, mean_estimation, lambda_mean,
                        lambda_mean_type, lambda_prec, lambda_prec_type, tune,
                        alpha_grid, theta_grid, random_search, p_random,
                        nfolds, symm_method, verbose) {
  x <- check_expression(x)
  n <- nrow(x)
  p <- ncol(x)
  nodes <- colnames(x)
  prior <- check_prior(prior, p, nodes)
  screening <- check_flag(screening, "screening")
  mean_estimation <- check_flag(mean_estimation, "mean_estimation")
  tune <- check_flag(tune, "tune")
  random_search <- check_flag(random_search, "random_search")
  verbose <- check_flag(verbose, "verbose")
  gamma <- check_number(gamma, "gamma", lower = 0, upper = 1)
  lambda_mean <- check_number(lambda_mean, "lambda_mean", lower = 0)
  lambda_prec <- check_number(lambda_prec, "lambda_prec", lower = 0)
  p_random <- check_number(p_random, "p_random", lower = 0, upper = 1,
                           null_ok = FALSE)
  nfolds <- check_number(nfolds, "nfolds", lower = 3, upper = n,
                         null_ok = FALSE)
  if (!is.numeric(alpha_grid) || length(alpha_grid) == 0L ||
      anyNA(alpha_grid) || any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("`alpha_grid` must contain values in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(theta_grid) || length(theta_grid) == 0L ||
      anyNA(theta_grid) || any(theta_grid <= 0)) {
    stop("`theta_grid` must contain positive values.", call. = FALSE)
  }

  enc <- encode_covariates(covariates, n)
  U <- enc$U
  q <- if (is.null(U)) 0L else ncol(U)

  # Random draws happen here, before any fit (see Reproducibility).
  foldid <- sample(rep(seq_len(nfolds), length.out = n))
  grid <- if (tune) {
    expand.grid(alpha = alpha_grid, theta = theta_grid)
  } else {
    data.frame(alpha = alpha_grid[1], theta = theta_grid[1])
  }
  combos <- lapply(seq_len(p), function(i) {
    if (tune && random_search && nrow(grid) > 1L) {
      size <- max(1L, floor(p_random * nrow(grid)))
      grid[sort(sample(nrow(grid), size)), , drop = FALSE]
    } else {
      grid
    }
  })

  if (mean_estimation) {
    mean_fit <- fit_mean(x, U, lambda_mean, lambda_mean_type, foldid,
                         verbose)
    z <- mean_fit$z
    mean_effects <- mean_fit$effects
  } else {
    z <- x
    mean_effects <- NULL
  }

  keep <- NULL
  if (!is.null(screening_matrix)) {
    screening_matrix <- check_square(screening_matrix, p, nodes,
                                     "screening_matrix")
    if (!all(screening_matrix %in% c(0, 1))) {
      stop("`screening_matrix` must contain only 0 and 1.", call. = FALSE)
    }
    keep <- screening_matrix == 1
    screening <- TRUE
  } else if (screening) {
    r <- abs(cor(z))
    if (is.null(gamma)) {
      gamma <- unname(quantile(r[upper.tri(r)], probs = 0.2))
    }
    keep <- r > gamma
    vmessage(verbose, "Screening threshold (gamma): ", signif(gamma, 3))
  }
  if (!is.null(keep)) {
    if (!is.null(prior)) {
      keep <- keep | prior > 0
    }
    diag(keep) <- FALSE
    dimnames(keep) <- list(nodes, nodes)
  }

  list(
    n = n, p = p, q = q, nodes = nodes,
    z = z, U = U, covariate_info = enc$info,
    mean_effects = mean_effects,
    penalty = if (is.null(prior)) NULL else 1 - prior,
    keep = keep, foldid = foldid, combos = combos,
    lambda_prec = lambda_prec, lambda_prec_type = lambda_prec_type,
    settings = list(
      prior = !is.null(prior), screening = screening, gamma = gamma,
      mean_estimation = mean_estimation, lambda_mean = lambda_mean,
      lambda_mean_type = lambda_mean_type, lambda_prec = lambda_prec,
      lambda_prec_type = lambda_prec_type, tune = tune,
      alpha_grid = alpha_grid, theta_grid = theta_grid,
      random_search = random_search, p_random = p_random, nfolds = nfolds,
      symm_method = symm_method),
    verbose = verbose
  )
}

run_nodes <- function(prep, ncores, verbose) {
  p <- prep$p
  if (ncores > 1L) {
    vmessage(verbose, "Fitting ", p, " node-wise regressions on ", ncores,
             " cores.")
    cl <- parallel::makeCluster(min(ncores, p))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    return(parallel::parLapply(cl, seq_len(p), fit_node, prep = prep))
  }
  lapply(seq_len(p), function(i) {
    vmessage(verbose, "Fitting node ", i, " of ", p, ".")
    fit_node(i, prep)
  })
}

# Combines the node-wise results into a "prcggm" object.
assemble_fit <- function(prep, results) {
  p <- prep$p
  q <- prep$q
  nodes <- prep$nodes
  results <- results[order(vapply(results, `[[`, numeric(1), "node"))]
  components <- c("Baseline", if (q > 0) colnames(prep$U))
  coefficients <- do.call(rbind, lapply(results, `[[`, "coef_row"))
  colnames(coefficients) <- c(
    "(Intercept)", nodes,
    unlist(lapply(components[-1], function(h) paste0(h, ":", nodes))))
  rownames(coefficients) <- nodes
  sigma2 <- stats::setNames(vapply(results, `[[`, numeric(1), "sigma2"),
                            nodes)
  networks <- list()
  precision <- list()
  for (h in seq_along(components)) {
    delta <- coefficients[, 1L + (h - 1L) * p + seq_len(p), drop = FALSE]
    dimnames(delta) <- list(nodes, nodes)
    method <- prep$settings$symm_method
    networks[[components[h]]] <- symmetrize(delta, method)
    # Omega_jk = -delta_jk / Var(eps_j), divided row-wise
    precision[[components[h]]] <- symmetrize(-delta / sigma2, method,
                                             ref = delta)
  }
  diag(precision$Baseline) <- 1 / sigma2
  tuning <- do.call(rbind, lapply(results, `[[`, "tuning"))
  tuning <- data.frame(node = nodes, tuning,
                       time = vapply(results, `[[`, numeric(1), "time"),
                       row.names = NULL)
  structure(
    list(
      networks = networks,
      precision = precision,
      mean_effects = prep$mean_effects,
      residuals = prep$z,
      coefficients = coefficients,
      sigma2 = sigma2,
      tuning = tuning,
      covariates = prep$covariate_info,
      settings = prep$settings,
      n = prep$n, p = p, q = q,
      call = NULL
    ),
    class = "prcggm"
  )
}
