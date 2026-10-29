#' Simulate data from a conditional Gaussian graphical model
#'
#' Generates expression data whose precision matrix depends linearly on the
#' covariates, \eqn{\Omega(x) = \Delta_0 + \sum_h x_h \Delta_h}, together with a
#' prior weight matrix of the chosen quality. The baseline network
#' \eqn{\Delta_0} is scale-free (preferential attachment) and each covariate
#' effect \eqn{\Delta_h} is an Erdos-Renyi random graph, following the
#' simulation design of Mapelli et al. (2026).
#'
#' @param n Number of subjects.
#' @param p Number of features (nodes).
#' @param n_binary,n_continuous Number of binary (Bernoulli(0.5), returned as
#'   factors with levels `"0"` and `"1"`) and standard normal covariates.
#' @param edge_prob Edge probability of the covariate-effect networks.
#' @param coef_range Range of the absolute values of the non-zero
#'   off-diagonal entries; signs are random.
#' @param power Preferential attachment power of the baseline network.
#' @param prior Quality of the simulated prior: `"perfect"` (weight 1 on true
#'   baseline edges, 0 elsewhere), `"noisy"` (true edges kept with
#'   probability `sensitivity`, false edges added with probability
#'   `1 - specificity`, with random confidence scores) or `"none"`.
#' @param sensitivity,specificity Parameters of the noisy prior.
#' @param confidence_range Range of the confidence scores of true edges in
#'   the noisy prior; false edges get scores between 0.1 and
#'   `confidence_range[1] - 0.1`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`x`}{The \eqn{n \times p}{n x p} matrix of simulated expression
#'     data.}
#'   \item{`covariates`}{Data frame of the simulated covariates, named `x1`,
#'     `x2`, ...: first the binary ones (factors), then the continuous ones.}
#'   \item{`prior`}{The \eqn{p \times p}{p x p} prior weight matrix.}
#'   \item{`truth_precision`}{Named list of the true precision components
#'     \eqn{\Delta_0, \Delta_1, \ldots}{Delta_0, Delta_1, ...}. Compare with
#'     the `precision` element of a fitted [prcggm()] model.}
#'   \item{`truth_networks`}{Named list of the true networks on the scale of
#'     the `networks` element of a fitted [prcggm()] model, that is the
#'     partial correlations
#'     \eqn{-(\Delta_h)_{jk} / \sqrt{(\Delta_0)_{jj} (\Delta_0)_{kk}}}{-(Delta_h)_jk / sqrt((Delta_0)_jj (Delta_0)_kk)}
#'     with zero diagonal. Their sign is opposite to that of the precision
#'     components.}
#' }
#' Both lists use the same names as the components of a fitted [prcggm()]
#' model (`"Baseline"`, `"x1_1"`, `"x2"`, ...), so that, for example,
#' `fit$networks$x1_1` estimates `sim$truth_networks$x1_1`.
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
                                 edge_prob = 0.1,
                                 coef_range = c(0.35, 0.5),
                                 power = 2.5,
                                 prior = c("perfect", "noisy", "none"),
                                 sensitivity = 0.8,
                                 specificity = 0.9,
                                 confidence_range = c(0.7, 1)) {
  prior <- match.arg(prior)
  check_number(n, "n", lower = 2, null_ok = FALSE)
  check_number(p, "p", lower = 3, null_ok = FALSE)
  nodes <- paste0("V", seq_len(p))

  random_weights <- function(adj) {
    w <- matrix(0, p, p)
    up <- which(upper.tri(adj) & adj)
    w[up] <- stats::runif(length(up), coef_range[1], coef_range[2]) *
      sample(c(-1, 1), length(up), replace = TRUE)
    w + t(w)
  }
  make_pd <- function(m) {
    min_eig <- min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig <= 0) {
      diag(m) <- diag(m) + abs(min_eig) + 0.1
    }
    m
  }

  # Baseline: preferential attachment tree, probability ~ degree^power + 1
  adj0 <- matrix(FALSE, p, p)
  degree <- numeric(p)
  for (t in 2:p) {
    target <- sample.int(t - 1, 1, prob = degree[seq_len(t - 1)]^power + 1)
    adj0[t, target] <- adj0[target, t] <- TRUE
    degree[c(t, target)] <- degree[c(t, target)] + 1
  }
  delta0 <- random_weights(adj0)
  diag(delta0) <- 0.1
  delta0 <- make_pd(delta0)

  n_cov <- n_binary + n_continuous
  covariates <- data.frame(row.names = seq_len(n))
  truth_precision <- list(Baseline = delta0)
  for (h in seq_len(n_cov)) {
    name <- paste0("x", h)
    adj <- matrix(stats::runif(p * p) < edge_prob, p, p)
    adj[lower.tri(adj, diag = TRUE)] <- FALSE
    adj <- adj | t(adj)
    if (h <= n_binary) {
      covariates[[name]] <- stats::rbinom(n, 1, 0.5)
      truth_precision[[paste0(name, "_1")]] <- random_weights(adj)
    } else {
      covariates[[name]] <- stats::rnorm(n)
      truth_precision[[name]] <- random_weights(adj)
    }
  }
  truth_precision <- lapply(truth_precision, function(m) {
    dimnames(m) <- list(nodes, nodes)
    m
  })
  # Partial-correlation scale, comparable with prcggm()$networks
  scale0 <- sqrt(outer(diag(delta0), diag(delta0)))
  truth_networks <- lapply(truth_precision, function(m) {
    r <- -m / scale0
    diag(r) <- 0
    r
  })

  x <- matrix(0, n, p, dimnames = list(NULL, nodes))
  values <- as.matrix(covariates)
  for (i in seq_len(n)) {
    omega <- delta0
    for (h in seq_len(n_cov)) {
      omega <- omega + values[i, h] * truth_precision[[h + 1]]
    }
    omega <- make_pd(omega)
    x[i, ] <- backsolve(chol(omega), stats::rnorm(p))
  }
  for (h in seq_len(n_binary)) {
    covariates[[h]] <- factor(covariates[[h]], levels = c(0, 1))
  }

  true_edges <- adj0
  weights <- matrix(0, p, p, dimnames = list(nodes, nodes))
  if (prior == "perfect") {
    weights[true_edges] <- 1
  } else if (prior == "noisy") {
    for (i in seq_len(p - 1)) {
      for (j in (i + 1):p) {
        if (true_edges[i, j]) {
          if (stats::runif(1) < sensitivity) {
            weights[i, j] <- weights[j, i] <-
              stats::runif(1, confidence_range[1], confidence_range[2])
          }
        } else if (stats::runif(1) > specificity) {
          weights[i, j] <- weights[j, i] <-
            stats::runif(1, 0.1, confidence_range[1] - 0.1)
        }
      }
    }
  }
  list(x = x, covariates = covariates, prior = weights,
       truth_precision = truth_precision, truth_networks = truth_networks)
}
