#' Predict personalised networks
#'
#' Computes the personalised network of new individuals as the linear
#' combination of the estimated baseline and covariate-specific components,
#' \eqn{\hat A(x) = \hat A_0 + \sum_h x_h \hat A_h}.
#'
#' Covariates in `newdata` are encoded exactly as in the estimation data:
#' numeric covariates are standardised with the training mean and standard
#' deviation, and factors are dummy coded against the training reference
#' level.
#'
#' @param object A `"prcggm"` object returned by [prcggm()].
#' @param newdata Data frame with the covariates of the new individuals,
#'   containing every covariate used in the fit. If `NULL` (default), or if
#'   the model has no covariates, the baseline network is returned.
#' @param type `"network"` (default) for the adjacency matrices of edge
#'   weights, or `"precision"` for the precision matrices.
#' @param ... Not used.
#'
#' @return A list with one matrix per row of `newdata`, named
#'   `Subject_1`, `Subject_2`, and so on. Without `newdata`, a list holding
#'   the baseline matrix only.
#' @seealso [prcggm()], [compare_networks()].
#'
#' @examples
#' data(prcggm_example)
#' set.seed(1)
#' fit <- prcggm(prcggm_example$x,
#'               covariates = prcggm_example$covariates,
#'               prior = prcggm_example$prior,
#'               tune = FALSE)
#' new_subjects <- data.frame(x1 = factor(c("0", "1")))
#' nets <- predict(fit, new_subjects)
#' compare_networks(nets$Subject_1, nets$Subject_2, method = "frobenius")
#' @export
predict.prcggm <- function(object, newdata = NULL,
                           type = c("network", "precision"), ...) {
  type <- match.arg(type)
  components <- if (type == "network") object$networks else object$precision
  if (is.null(newdata) || is.null(object$covariates)) {
    return(list(Baseline = components$Baseline))
  }
  if (is.matrix(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  U <- covariate_design(newdata, object$covariates)
  out <- lapply(seq_len(nrow(U)), function(i) {
    net <- components$Baseline
    for (h in colnames(U)) {
      net <- net + U[i, h] * components[[h]]
    }
    net
  })
  names(out) <- paste0("Subject_", seq_len(nrow(U)))
  out
}

#' @export
print.prcggm <- function(x, ...) {
  cat("Prior-informed conditional Gaussian graphical model\n")
  cat("  Subjects:", x$n, " Features:", x$p, " Covariates:", x$q,
      "\n")
  cat("  Prior knowledge:", if (x$settings$prior) "yes" else "no",
      " Symmetrisation:", x$settings$symm_method, "\n")
  cat("  Selection:", selection_label(x$settings), "\n")
  cat("  Network components:", paste(names(x$networks), collapse = ", "),
      "\n")
  invisible(x)
}

#' Summarise a fitted prCGGM model
#'
#' Gives an overview of a fitted `"prcggm"` object: the size of the problem,
#' the model that was fitted, the structure of each estimated network
#' component and the hyperparameters selected across nodes.
#'
#' The summary reports:
#' * the number of subjects \eqn{n}, features \eqn{p} and encoded covariate
#'   columns \eqn{q}, whether prior knowledge was used and the symmetrisation
#'   rule;
#' * for each network component (the baseline network and one component per
#'   covariate effect), the number of edges, the density (edges divided by
#'   the \eqn{p(p - 1)/2}{p(p - 1)/2} possible pairs) and the mean absolute
#'   weight of the non-zero edges;
#' * how often each combination of \eqn{\alpha}{alpha} and
#'   \eqn{\vartheta}{theta} was selected by BIC across the \eqn{p} node-wise
#'   regressions.
#'
#' Printing the returned object shows these results; the full estimates
#' remain available in the fitted object (see [prcggm()]).
#'
#' @param object A `"prcggm"` object returned by [prcggm()] or
#'   [prcggm_slurm_collect()].
#' @param ... Not used.
#'
#' @return An object of class `"summary.prcggm"`, a list with elements:
#' \describe{
#'   \item{`components`}{Data frame with one row per network component and
#'     columns `component`, `edges`, `density` and `mean_abs_weight`.}
#'   \item{`tuning`}{Contingency table of the selected `alpha` (rows) and
#'     `theta` (columns) across nodes.}
#'   \item{`fit`}{List with `n`, `p`, `q`, `prior` (logical),
#'     `symm_method` and `selection` (the selection rule, with `refit`).}
#' }
#'
#' @seealso [prcggm()], [plot.prcggm()], [predict.prcggm()].
#'
#' @examples
#' data(prcggm_example)
#' set.seed(1)
#' fit_tuned <- prcggm(prcggm_example$x,
#'               covariates = prcggm_example$covariates,
#'               prior = prcggm_example$prior,
#'               alpha_grid = c(0.5, 0.9), theta_grid = c(1, 1.3))
#' summary(fit_tuned)
#' @export
summary.prcggm <- function(object, ...) {
  pairs <- object$p * (object$p - 1) / 2
  components <- do.call(rbind, lapply(names(object$networks), function(h) {
    w <- object$networks[[h]][upper.tri(object$networks[[h]])]
    nz <- w != 0
    data.frame(component = h,
               edges = sum(nz),
               density = sum(nz) / pairs,
               mean_abs_weight = if (any(nz)) mean(abs(w[nz])) else 0)
  }))
  out <- list(
    components = components,
    tuning = table(alpha = object$tuning$alpha, theta = object$tuning$theta),
    fit = list(n = object$n, p = object$p, q = object$q,
               prior = object$settings$prior,
               symm_method = object$settings$symm_method,
               selection = selection_label(object$settings))
  )
  class(out) <- "summary.prcggm"
  out
}

#' @export
print.summary.prcggm <- function(x, ...) {
  cat("Prior-informed conditional Gaussian graphical model\n")
  cat("  Subjects:", x$fit$n, " Features:", x$fit$p,
      " Covariates:", x$fit$q, "\n")
  cat("  Prior knowledge:", if (x$fit$prior) "yes" else "no",
      " Symmetrisation:", x$fit$symm_method, "\n")
  if (!is.null(x$fit$selection)) {
    cat("  Selection:", x$fit$selection, "\n")
  }
  cat("\n")
  cat("Network components:\n")
  print(x$components, row.names = FALSE, digits = 3)
  cat("\nSelected (alpha, theta) across nodes:\n")
  print(x$tuning)
  invisible(x)
}

#' Plot a network component of a fitted model
#'
#' Draws one network component of a fitted model, or the personalised
#' network of one individual, as a heatmap with [plot_network()].
#'
#' @param x A `"prcggm"` object returned by [prcggm()].
#' @param component Name of the component to draw: `"Baseline"` (default) or
#'   one of the covariate effects in `names(x$networks)`. Ignored if
#'   `newdata` is supplied.
#' @param newdata Optional one-row data frame of covariates; if supplied, the
#'   personalised network of that individual is drawn.
#' @param ... Further arguments passed to [plot_network()], such as
#'   `threshold` or `main`.
#'
#' @return Called for its side effect; returns `x` invisibly.
#'
#' @examples
#' data(prcggm_example)
#' set.seed(1)
#' fit <- prcggm(prcggm_example$x,
#'               covariates = prcggm_example$covariates,
#'               prior = prcggm_example$prior,
#'               tune = FALSE)
#' plot(fit)
#' plot(fit, component = "x1_1")
#' @export
plot.prcggm <- function(x, component = "Baseline", newdata = NULL, ...) {
  if (!is.null(newdata)) {
    if (nrow(newdata) != 1L) {
      stop("`newdata` must have exactly one row.", call. = FALSE)
    }
    network <- predict(x, newdata)[[1]]
    main <- "Personalised network"
  } else {
    if (!component %in% names(x$networks)) {
      stop("`component` must be one of: ",
           paste(names(x$networks), collapse = ", "), ".", call. = FALSE)
    }
    network <- x$networks[[component]]
    main <- component
  }
  args <- list(...)
  if (is.null(args$main)) {
    args$main <- main
  }
  do.call(plot_network, c(list(network = network), args))
  invisible(x)
}

# Description of the selection rule of a fit. Fits made with version 0.1.0
# have no `selection` setting: they used cross-validation without refit.
selection_label <- function(settings) {
  if (!is.null(settings$lambda_prec)) {
    return(paste0("fixed lambda (", signif(settings$lambda_prec, 3), ")",
                  if (isTRUE(settings$refit)) ", OLS refit" else ""))
  }
  selection <- settings$selection
  refit <- settings$refit
  if (is.null(selection)) {
    selection <- paste0("cv_", if (is.null(settings$lambda_prec_type)) "min"
                        else settings$lambda_prec_type)
    refit <- FALSE
  }
  paste0(selection, if (isTRUE(refit) || grepl("^stability", selection))
    ", OLS refit" else ", penalised coefficients")
}
