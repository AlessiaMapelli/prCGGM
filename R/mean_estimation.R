# Step 1: node-wise regression of each feature on the covariates.
# Returns the residuals Z and the matrix of covariate effects on the mean.
fit_mean <- function(x, U, lambda_mean, lambda_mean_type, foldid, verbose) {
  n <- nrow(x)
  p <- ncol(x)
  nodes <- colnames(x)
  if (is.null(U)) {
    vmessage(verbose, "No covariates: centring the features.")
    effects <- matrix(colMeans(x), p, 1,
                      dimnames = list(nodes, "(Intercept)"))
    z <- sweep(x, 2, effects[, 1])
    return(list(z = z, effects = effects))
  }
  q <- ncol(U)
  effects <- matrix(0, p, q + 1,
                    dimnames = list(nodes, c("(Intercept)", colnames(U))))
  z <- x
  s <- if (lambda_mean_type == "min") "lambda.min" else "lambda.1se"
  for (j in seq_len(p)) {
    y <- x[, j]
    if (q == 1L) {
      # glmnet needs at least two predictors: use least squares instead
      beta <- coef(lm(y ~ U[, 1]))
    } else if (is.null(lambda_mean)) {
      mod <- glmnet::cv.glmnet(x = U, y = y, foldid = foldid)
      beta <- coef(mod, s = s)[, 1]
    } else {
      mod <- glmnet::glmnet(x = U, y = y, lambda = lambda_mean)
      beta <- coef(mod, s = lambda_mean)[, 1]
    }
    effects[j, ] <- beta
    z[, j] <- y - drop(cbind(1, U) %*% beta)
  }
  vmessage(verbose, "Mean estimation completed.")
  list(z = z, effects = effects)
}
