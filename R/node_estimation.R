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

# Default settings of the selection rules (see prcggm(), `selection_control`).
selection_defaults <- function() {
  list(ebic_gamma = 0.5, stab_B = 50, stab_fraction = 0.5, stab_pi = 0.75,
       stab_pfer = 1)
}

# OLS refit of y on the columns of X in `support`; returns the coefficients
# (zero outside the support).
refit_ols <- function(X, y, support) {
  coef <- numeric(ncol(X))
  if (any(support)) {
    b <- stats::lm.fit(X[, support, drop = FALSE], y)$coefficients
    b[is.na(b)] <- 0
    coef[support] <- b
  }
  coef
}

# Residual sum of squares of the OLS refit on each distinct support of a
# lambda path (columns of the logical matrix S), from the Gram matrix:
# RSS = y'y - b'X_s'y with b = (X_s'X_s)^-1 X_s'y. Falls back to a QR fit if
# X_s'X_s is singular. Returns one value per column of S.
relaxed_rss <- function(S, gram) {
  one <- function(s) {
    if (!any(s)) {
      return(gram$yty)
    }
    R <- tryCatch(chol(gram$XtX[s, s, drop = FALSE]), error = function(e) NULL)
    if (is.null(R)) {
      return(sum(stats::lm.fit(gram$X[, s, drop = FALSE], gram$y)$residuals^2))
    }
    b <- backsolve(R, forwardsolve(t(R), gram$Xty[s]))
    max(gram$yty - sum(gram$Xty[s] * b), 0)
  }
  key <- apply(S, 2, function(s) paste(which(s), collapse = ","))
  ukey <- unique(key)
  rss <- vapply(ukey, function(u) one(S[, match(u, key)]), numeric(1))
  unname(rss[match(key, ukey)])
}

# Stability selection (Meinshausen and Buhlmann, 2010): refits the lambda
# path on each subsample, restricts the path to the region where the average
# number of selected predictors is at most sqrt(pfer * (2 * pi - 1) * P) and
# returns the predictors with selection frequency >= pi in that region.
stability_support <- function(args, lambda, X, y, subsamples, ctrl) {
  P <- ncol(X)
  counts <- matrix(0, P, length(lambda))
  for (b in seq_len(nrow(subsamples))) {
    s <- subsamples[b, ]
    a <- args
    a$x <- X[s, , drop = FALSE]
    a$y <- y[s]
    a$lambda <- lambda
    fb <- do.call(sparsegl::sparsegl, a)
    L <- ncol(fb$beta)
    counts[, seq_len(L)] <- counts[, seq_len(L)] + (as.matrix(fb$beta) != 0)
  }
  freq <- counts / nrow(subsamples)
  q_lambda <- sqrt(ctrl$stab_pfer * (2 * ctrl$stab_pi - 1) * P)
  first_over <- which(colSums(freq) > q_lambda)[1]
  region <- if (is.na(first_over)) seq_along(lambda) else
    seq_len(first_over - 1L)
  if (length(region) == 0L) {
    return(rep(FALSE, P))
  }
  apply(freq[, region, drop = FALSE], 1, max) >= ctrl$stab_pi
}

# Fits the node-wise regression for node i. Returns the coefficient row
# (intercept followed by one block of p coefficients per network component),
# the residual variance and the selected tuning parameters.
fit_node <- function(i, prep) {
  time_start <- Sys.time()
  n <- prep$n
  p <- prep$p
  q <- prep$q
  selection <- prep$selection
  ctrl <- prep$selection_control
  y <- prep$z[, i]
  others <- node_neighbours(i, prep)
  m <- length(others)
  coef_row <- numeric(1 + (q + 1) * p)
  if (m == 0L) {
    rss <- sum(y^2)
    return(list(
      node = i, coef_row = coef_row, coef_row_penalised = coef_row,
      sigma2 = rss / n,
      tuning = data.frame(lambda = NA_real_, alpha = NA_real_,
                          theta = NA_real_, bic = n * log(rss / n), df = 0L),
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
  P <- ncol(X)
  pf_baseline <- if (is.null(prep$penalty)) rep(1, m) else prep$penalty[i, others]
  # Without covariates: (prior-weighted) lasso, one group per predictor
  group <- if (q == 0) seq_len(m) else rep(seq_len(q + 1), each = m)
  pf_group <- if (q == 0) rep(1, m) else c(0, rep(1, q))
  fixed_lambda <- !is.null(prep$lambda_prec)
  use_cv <- !fixed_lambda &&
    selection %in% c("cv_min", "cv_1se", "stability_cv")
  use_relaxed <- !fixed_lambda &&
    selection %in% c("bic_relaxed", "ebic_relaxed", "stability_bic")
  combos <- if (q == 0) data.frame(alpha = 1, theta = 1) else prep$combos[[i]]

  # Every candidate edge has full prior support and there are no covariates:
  # no penalty at all, ordinary least squares
  if (q == 0 && all(pf_baseline == 0)) {
    coef <- refit_ols(X, y, rep(TRUE, P))
    rss <- sum((y - drop(X %*% coef))^2)
    df <- sum(coef != 0)
    idx <- 1L + others
    coef_row[idx] <- coef
    return(list(
      node = i, coef_row = coef_row, coef_row_penalised = coef_row,
      sigma2 = rss / max(n - df, 1),
      tuning = data.frame(lambda = 0, alpha = 1, theta = 1,
                          bic = n * log(rss / n) + log(n) * df, df = df),
      time = as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    ))
  }

  # Gram matrix for the relaxed criteria (once per node)
  gram <- if (use_relaxed) {
    list(X = X, y = y, XtX = crossprod(X), Xty = drop(crossprod(X, y)),
         yty = sum(y^2))
  }

  # Lambda path (and cross-validation if needed) for each combination
  fits <- vector("list", nrow(combos))
  for (k in seq_len(nrow(combos))) {
    args <- list(x = X, y = y, group = group, pf_group = pf_group,
                 pf_sparse = if (q == 0) pf_baseline else
                   c(pf_baseline, rep(combos$theta[k], q * m)),
                 asparse = combos$alpha[k], intercept = FALSE)
    if (fixed_lambda) {
      args$lambda <- prep$lambda_prec * c(1.01, 1, 0.99)
    }
    if (use_cv) {
      cv <- do.call(sparsegl::cv.sparsegl, c(args, list(foldid = prep$foldid)))
      path <- cv$sparsegl.fit
    } else {
      cv <- NULL
      path <- do.call(sparsegl::sparsegl, args)
    }
    B <- as.matrix(path$beta)
    rss <- colSums((y - X %*% B)^2)
    df <- colSums(B != 0)
    bic <- n * log(rss / n) + log(n) * df
    fits[[k]] <- list(
      args = args, lambda = path$lambda, B = B, bic = bic,
      ebic = bic + 2 * ctrl$ebic_gamma * lchoose(P, df),
      idx_min = if (use_cv) which(path$lambda == cv$lambda.min)[1] else NA,
      idx_1se = if (use_cv) which(path$lambda == cv$lambda.1se)[1] else NA,
      idx_fixed = if (fixed_lambda) 2L else NA)
    if (use_relaxed) {
      # BIC / EBIC on the RSS of the OLS refit of each support
      bic_rel <- n * log(relaxed_rss(B != 0, gram) / n) + log(n) * df
      fits[[k]]$bic_relaxed <- bic_rel
      fits[[k]]$ebic_relaxed <- bic_rel + 2 * ctrl$ebic_gamma * lchoose(P, df)
    }
  }

  # Selected combination k and lambda index
  by_criterion <- function(field) {
    k <- unname(which.min(vapply(fits, function(f) min(f[[field]]),
                                 numeric(1))))
    c(k = k, idx = unname(which.min(fits[[k]][[field]])))
  }
  by_index <- function(field) {
    k <- unname(which.min(vapply(fits, function(f) f$bic[f[[field]]],
                                 numeric(1))))
    c(k = k, idx = unname(fits[[k]][[field]]))
  }
  pick <- if (fixed_lambda) {
    by_index("idx_fixed")
  } else {
    switch(selection,
           ebic_relaxed = by_criterion("ebic_relaxed"),
           bic_relaxed = , stability_bic = by_criterion("bic_relaxed"),
           bic = by_criterion("bic"),
           ebic = by_criterion("ebic"),
           cv_min = , stability_cv = by_index("idx_min"),
           cv_1se = by_index("idx_1se"))
  }
  f <- fits[[pick[["k"]]]]
  coef_pen <- f$B[, pick[["idx"]]]
  lambda <- f$lambda[pick[["idx"]]]
  support <- coef_pen != 0
  refit <- prep$refit
  if (!fixed_lambda && selection %in% c("stability_bic", "stability_cv")) {
    support <- stability_support(f$args, f$lambda, X, y, prep$subsamples, ctrl)
    lambda <- NA_real_
    refit <- TRUE
  }
  coef <- if (refit) refit_ols(X, y, support) else coef_pen
  rss <- sum((y - drop(X %*% coef))^2)
  df <- sum(coef != 0)

  # Columns of the coefficient row: intercept, then for each component h the
  # block (1 + h * p + 1):(1 + (h + 1) * p) indexed by node.
  idx <- unlist(lapply(0:q, function(h) 1L + h * p + others))
  coef_row[idx] <- coef
  # Penalised coefficients of the selected model (not defined for stability
  # selection)
  coef_row_penalised <- NULL
  if (!(!fixed_lambda && selection %in% c("stability_bic", "stability_cv"))) {
    coef_row_penalised <- numeric(length(coef_row))
    coef_row_penalised[idx] <- coef_pen
  }
  list(
    node = i,
    coef_row = coef_row,
    coef_row_penalised = coef_row_penalised,
    sigma2 = rss / max(n - df, 1),
    tuning = data.frame(lambda = lambda,
                        alpha = combos$alpha[pick[["k"]]],
                        theta = combos$theta[pick[["k"]]],
                        bic = n * log(rss / n) + log(n) * df,
                        df = df),
    time = as.numeric(difftime(Sys.time(), time_start, units = "mins"))
  )
}
