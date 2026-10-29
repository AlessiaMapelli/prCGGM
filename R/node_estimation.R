# Step 2: node-wise prior-informed sparse-group lasso regressions.

# Predictors kept for node i: all other nodes, minus the screened-out ones.
# Prior-supported pairs are never screened out.
node_neighbours <- function(i, prep) {
  others <- seq_len(prep$p)[-i]
  if (!is.null(prep$keep)) {
    others <- others[prep$keep[i, others]]
  }
  others
}

# Fits the node-wise regression for node i. Returns the coefficient row
# (intercept followed by one block of p coefficients per network component),
# the residual variance and the selected tuning parameters.
fit_node <- function(i, prep) {
  time_start <- Sys.time()
  n <- prep$n
  p <- prep$p
  q <- prep$q
  y <- prep$z[, i]
  others <- node_neighbours(i, prep)
  m <- length(others)
  coef_row <- numeric(1 + (q + 1) * p)
  if (m == 0L) {
    rss <- sum(y^2)
    return(list(
      node = i, coef_row = coef_row, sigma2 = rss / n,
      tuning = data.frame(lambda = NA_real_, alpha = NA_real_,
                          theta = NA_real_, bic = n * log(rss / n)),
      time = as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    ))
  }
  zo <- prep$z[, others, drop = FALSE]
  X <- zo
  if (q > 0) {
    X <- cbind(zo, do.call(cbind, lapply(seq_len(q), function(h) {
      prep$U[, h] * zo
    })))
  }
  group <- rep(seq_len(q + 1), each = m)
  pf_baseline <- if (is.null(prep$penalty)) rep(1, m) else prep$penalty[i, others]
  s <- if (prep$lambda_prec_type == "min") "lambda.min" else "lambda.1se"
  lambda <- NULL
  if (!is.null(prep$lambda_prec)) {
    lambda <- prep$lambda_prec * c(1.01, 1, 0.99)
  }

  fit_one <- function(alpha, theta) {
    if (q == 0 && all(pf_baseline == 0)) {
      # every candidate edge has full prior support: no penalty at all
      beta <- c(0, stats::lm.fit(X, y)$coefficients)
      beta[is.na(beta)] <- 0
      rss <- sum((y - drop(cbind(1, X) %*% beta))^2)
      df <- sum(beta != 0)
      return(list(beta = beta, lambda = 0, sigma2 = rss / max(n - df, 1),
                  bic = n * log(rss / n) + log(n) * df))
    }
    if (q == 0 && is.null(prep$penalty) && m >= 2L) {
      mod <- glmnet::cv.glmnet(x = X, y = y, lambda = lambda,
                               foldid = prep$foldid)
    } else if (q == 0) {
      mod <- sparsegl::cv.sparsegl(x = X, y = y, pf_sparse = pf_baseline,
                                   asparse = 1, lambda = lambda,
                                   intercept = FALSE, foldid = prep$foldid)
    } else {
      mod <- sparsegl::cv.sparsegl(
        x = X, y = y, group = group,
        pf_group = c(0, rep(1, q)),
        pf_sparse = c(pf_baseline, rep(theta, q * m)),
        asparse = alpha, lambda = lambda,
        intercept = FALSE, foldid = prep$foldid)
    }
    beta <- as.numeric(coef(mod, s = s)[, 1])
    rss <- sum((y - drop(cbind(1, X) %*% beta))^2)
    df <- sum(beta != 0)
    list(beta = beta,
         lambda = if (s == "lambda.min") mod$lambda.min else mod$lambda.1se,
         sigma2 = rss / max(n - df, 1),
         bic = n * log(rss / n) + log(n) * df)
  }

  combos <- if (q == 0) data.frame(alpha = 1, theta = 1) else prep$combos[[i]]
  best <- NULL
  best_combo <- 1L
  for (k in seq_len(nrow(combos))) {
    fit <- fit_one(combos$alpha[k], combos$theta[k])
    if (is.null(best) || fit$bic < best$bic) {
      best <- fit
      best_combo <- k
    }
  }
  # Columns of the coefficient row: intercept, then for each component h the
  # block (1 + h * p + 1):(1 + (h + 1) * p) indexed by node.
  idx <- c(1L, unlist(lapply(0:q, function(h) 1L + h * p + others)))
  coef_row[idx] <- best$beta
  list(
    node = i,
    coef_row = coef_row,
    sigma2 = best$sigma2,
    tuning = data.frame(lambda = best$lambda,
                        alpha = combos$alpha[best_combo],
                        theta = combos$theta[best_combo],
                        bic = best$bic),
    time = as.numeric(difftime(Sys.time(), time_start, units = "mins"))
  )
}
