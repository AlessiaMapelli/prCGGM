#' Simulate data from a conditional Gaussian graphical model
#'
#' Generates expression data whose precision matrix depends linearly on the
#' covariates, \eqn{\Omega(x) = \Delta_0 + \sum_h x_h \Delta_h}, together with a
#' prior weight matrix of the chosen quality, following the simulation design
#' of Mapelli et al. (2026).
#'
#' * **Baseline network** \eqn{\Delta_0}: `p` edges drawn from a fitness
#'   model (Goh et al., 2001), whose degree distribution follows a power law
#'   with exponent `exponent`; non-zero entries
#'   \eqn{\pm}{+-}U(`coef_range`).
#' * **Covariate networks** \eqn{\Delta_h}, one per covariate: Erdos-Renyi
#'   graphs with edge probability `edge_prob` (by default \eqn{1/(p - 1)},
#'   about \eqn{p/2} edges); non-zero entries \eqn{\pm}{+-}U(`effect_range`) and
#'   zero diagonal.
#' * **Diagonal**: one constant \eqn{d} for all nodes and all covariate values,
#'   the smallest that keeps the minimum eigenvalue of \eqn{\Omega(x)} at
#'   least `eig_margin` over the whole covariate range. Because
#'   \eqn{\Omega(x)} is linear in \eqn{x}, checking the \eqn{2^q} corners of
#'   the covariate range is enough, so \eqn{\Omega(x)} is positive definite for
#'   every subject. The conditional variances \eqn{1/d} do not depend on the
#'   covariates, so the node-wise regressions are exactly linear in the
#'   covariates, as assumed by [prcggm()].
#' * **Covariates**: binary ones are Bernoulli(0.5), continuous ones
#'   Uniform(0, 1). Every covariate has an effect.
#'
#' **Target of the baseline.** [prcggm()] standardises numeric covariates and
#' codes factors against their first level, so its baseline network is the
#' network of a subject with the continuous covariates at their mean and the
#' binary ones at the reference level 0:
#' \eqn{T_0 = \Delta_0 + \sum_{h \in continuous} \bar{x}_h \Delta_h}{T0 =
#' Delta0 + sum over continuous h of mean(x_h) * Delta_h}. The true baseline
#' returned in `truth_precision` and `truth_networks` is \eqn{T_0}, and the prior
#' is generated from the edges of \eqn{T_0}. Likewise, the true effects of
#' continuous covariates are returned per standard deviation of the covariate,
#' \eqn{sd(x_h) \Delta_h}, the scale of the estimates of [prcggm()].
#'
#' @param n Number of subjects.
#' @param p Number of features (nodes).
#' @param n_binary,n_continuous Number of binary (returned as factors with
#'   levels `"0"` and `"1"`) and continuous covariates.
#' @param edge_prob Edge probability of the covariate networks. If `NULL`
#'   (default), \eqn{1/(p - 1)}.
#' @param coef_range Range of the absolute values of the non-zero
#'   off-diagonal entries of the baseline network; signs are random.
#' @param effect_range Range of the absolute values of the non-zero entries
#'   of the covariate networks; signs are random.
#' @param exponent Exponent of the power-law degree distribution of the
#'   baseline network.
#' @param prior Quality of the simulated prior: `"perfect"` (weight 1 on the
#'   edges of the true baseline, 0 elsewhere), `"noisy"` (true edges kept with
#'   probability `sensitivity`, false edges added with probability
#'   `1 - specificity`, with random confidence scores) or `"none"`.
#' @param sensitivity,specificity Parameters of the noisy prior.
#' @param confidence_range Range of the confidence scores of true edges in
#'   the noisy prior; false edges get scores between 0.1 and
#'   `confidence_range[1] - 0.1`.
#' @param eig_margin Minimum eigenvalue of the precision matrix over the
#'   covariate range.
#' @param power Deprecated and ignored (preferential-attachment power of
#'   version 0.1.0); see `exponent`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`x`}{The \eqn{n \times p}{n x p} matrix of simulated expression
#'     data.}
#'   \item{`covariates`}{Data frame of the simulated covariates, named `x1`,
#'     `x2`, ...: first the binary ones (factors), then the continuous ones.}
#'   \item{`prior`}{The \eqn{p \times p}{p x p} prior weight matrix (all zero
#'     for `prior = "none"`).}
#'   \item{`truth_precision`}{Named list of the true precision components on
#'     the scale of the `precision` element of a fitted [prcggm()] model:
#'     `"Baseline"` (\eqn{T_0}, with diagonal \eqn{d}) and one component per
#'     covariate.}
#'   \item{`truth_networks`}{Named list of the true networks on the scale of
#'     the `networks` element of a fitted [prcggm()] model, the node-wise
#'     regression coefficients \eqn{-(\Delta_h)_{jk} / d}{-(Delta_h)_jk / d},
#'     with zero diagonal. Their sign is opposite to that of the precision
#'     components.}
#'   \item{`diagonal`}{The common diagonal \eqn{d}.}
#' }
#' Both lists use the same names as the components of a fitted [prcggm()]
#' model (`"Baseline"`, `"x1_1"`, `"x2"`, ...), so that, for example,
#' `fit$networks$x1_1` estimates `sim$truth_networks$x1_1`.
#'
#' @references
#' Goh, K.-I., Kahng, B. and Kim, D. (2001). Universal behavior of load
#' distribution in scale-free networks. *Physical Review Letters*, 87(27),
#' 278701.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_prcggm_data(n = 200, p = 8)
#' str(sim, max.level = 1)
#'
#' # Compare the estimated baseline network with the truth
#' fit <- prcggm(sim$x, covariates = sim$covariates, prior = sim$prior,
#'               mean_estimation = FALSE, tune = FALSE)
#' compare_networks(fit$networks$Baseline, sim$truth_networks$Baseline,
#'                  method = "correlation")
#' @export
simulate_prcggm_data <- function(n = 500,
                                 p = 10,
                                 n_binary = 1,
                                 n_continuous = 0,
                                 edge_prob = NULL,
                                 coef_range = c(0.35, 0.5),
                                 effect_range = c(0.175, 0.25),
                                 exponent = 2.5,
                                 prior = c("perfect", "noisy", "none"),
                                 sensitivity = 0.8,
                                 specificity = 0.9,
                                 confidence_range = c(0.7, 1),
                                 eig_margin = 0.1,
                                 power = NULL) {
  prior <- match.arg(prior)
  check_number(n, "n", lower = 2, null_ok = FALSE)
  check_number(p, "p", lower = 3, null_ok = FALSE)
  check_number(n_binary, "n_binary", lower = 0, null_ok = FALSE)
  check_number(n_continuous, "n_continuous", lower = 0, null_ok = FALSE)
  if (is.null(edge_prob)) edge_prob <- 1 / (p - 1)
  check_number(edge_prob, "edge_prob", lower = 0, upper = 1, null_ok = FALSE)
  check_number(exponent, "exponent", lower = 2, null_ok = FALSE)
  check_number(eig_margin, "eig_margin", lower = 0, null_ok = FALSE)
  if (!is.null(power)) {
    warning("`power` is deprecated and ignored; the baseline network now ",
            "has a power-law degree distribution with exponent `exponent`.",
            call. = FALSE)
  }
  nodes <- paste0("V", seq_len(p))
  ut <- upper.tri(diag(p))

  random_weights <- function(adj, range) {
    w <- matrix(0, p, p)
    up <- which(ut & adj)
    w[up] <- stats::runif(length(up), range[1], range[2]) *
      sample(c(-1, 1), length(up), replace = TRUE)
    w + t(w)
  }

  # Baseline: fitness model, P(edge jk) proportional to f_j f_k with
  # f_j = j^(-1 / (exponent - 1)), p edges without multiple edges
  fitness <- seq_len(p)^(-1 / (exponent - 1))
  pairs <- which(ut)
  n_edges <- min(p, length(pairs))
  prob <- outer(fitness, fitness)[pairs]
  adj0 <- matrix(FALSE, p, p)
  adj0[sample(pairs, n_edges, prob = prob)] <- TRUE
  adj0 <- adj0 | t(adj0)
  B0 <- random_weights(adj0, coef_range)

  # Covariates and covariate networks
  n_cov <- n_binary + n_continuous
  binary <- seq_len(n_cov) <= n_binary
  values <- matrix(0, n, n_cov)
  B <- vector("list", n_cov)
  for (h in seq_len(n_cov)) {
    values[, h] <- if (binary[h]) stats::rbinom(n, 1, 0.5) else stats::runif(n)
    adj <- matrix(stats::runif(p * p) < edge_prob, p, p)
    adj[!ut] <- FALSE
    B[[h]] <- random_weights(adj | t(adj), effect_range)
  }

  # Common diagonal: positive definite at all corners of the covariate range
  corner_eig <- function(u) {
    m <- B0
    for (h in seq_len(n_cov)) m <- m + u[h] * B[[h]]
    min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
  }
  min_eig <- if (n_cov == 0) corner_eig(numeric(0)) else
    min(apply(as.matrix(expand.grid(rep(list(0:1), n_cov))), 1, corner_eig))
  d <- eig_margin - min_eig
  delta0 <- B0 + diag(d, p)

  x <- matrix(0, n, p, dimnames = list(NULL, nodes))
  for (i in seq_len(n)) {
    omega <- delta0
    for (h in seq_len(n_cov)) {
      omega <- omega + values[i, h] * B[[h]]
    }
    x[i, ] <- backsolve(chol(omega), stats::rnorm(p))
  }

  # Truth on the scale of the estimates of prcggm()
  target <- delta0
  for (h in which(!binary)) target <- target + mean(values[, h]) * B[[h]]
  truth_precision <- list(Baseline = target)
  covariates <- data.frame(row.names = seq_len(n))
  for (h in seq_len(n_cov)) {
    name <- paste0("x", h)
    if (binary[h]) {
      covariates[[name]] <- factor(values[, h], levels = c(0, 1))
      truth_precision[[paste0(name, "_1")]] <- B[[h]]
    } else {
      covariates[[name]] <- values[, h]
      truth_precision[[name]] <- stats::sd(values[, h]) * B[[h]]
    }
  }
  truth_precision <- lapply(truth_precision, function(m) {
    dimnames(m) <- list(nodes, nodes)
    m
  })
  truth_networks <- lapply(truth_precision, function(m) {
    r <- -m / d
    diag(r) <- 0
    r
  })

  # Prior from the edges of the true baseline
  weights <- matrix(0, p, p, dimnames = list(nodes, nodes))
  true_edges <- abs(target[pairs]) > 1e-8
  if (prior == "perfect") {
    weights[pairs[true_edges]] <- 1
  } else if (prior == "noisy") {
    w <- numeric(length(pairs))
    keep <- true_edges & stats::runif(length(pairs)) < sensitivity
    add <- !true_edges & stats::runif(length(pairs)) > specificity
    w[keep] <- stats::runif(sum(keep), confidence_range[1],
                            confidence_range[2])
    w[add] <- stats::runif(sum(add), 0.1, confidence_range[1] - 0.1)
    weights[pairs] <- w
  }
  weights[lower.tri(weights)] <- t(weights)[lower.tri(weights)]

  list(x = x, covariates = covariates, prior = weights,
       truth_precision = truth_precision, truth_networks = truth_networks,
       diagonal = d)
}
