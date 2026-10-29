# Internal helpers: input validation and covariate encoding.

check_expression <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop("All columns of `x` must be numeric.", call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains missing values; impute or remove them first.",
         call. = FALSE)
  }
  if (ncol(x) < 3) {
    stop("`x` must contain at least 3 features (columns).", call. = FALSE)
  }
  if (nrow(x) < 10) {
    warning("The sample size is very small; the network estimate may be ",
            "unreliable.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  if (anyDuplicated(colnames(x))) {
    stop("Column names of `x` must be unique.", call. = FALSE)
  }
  x
}

check_square <- function(m, p, nodes, name) {
  m <- as.matrix(m)
  if (!is.numeric(m) || !identical(dim(m), c(p, p))) {
    stop("`", name, "` must be a numeric ", p, " x ", p, " matrix.",
         call. = FALSE)
  }
  if (anyNA(m)) {
    stop("`", name, "` contains missing values.", call. = FALSE)
  }
  if (!is.null(rownames(m)) && !identical(rownames(m), nodes)) {
    stop("Row names of `", name, "` must match the column names of `x`.",
         call. = FALSE)
  }
  if (!is.null(colnames(m)) && !identical(colnames(m), nodes)) {
    stop("Column names of `", name, "` must match the column names of `x`.",
         call. = FALSE)
  }
  dimnames(m) <- list(nodes, nodes)
  m
}

check_prior <- function(prior, p, nodes) {
  if (is.null(prior)) {
    return(NULL)
  }
  prior <- check_square(prior, p, nodes, "prior")
  if (any(prior < 0)) {
    warning("`prior` contains negative values: setting them to 0.",
            call. = FALSE)
    prior[prior < 0] <- 0
  }
  if (any(prior > 1)) {
    warning("`prior` contains values greater than 1: dividing by the maximum.",
            call. = FALSE)
    prior <- prior / max(prior)
  }
  diag(prior) <- 0
  prior
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

check_number <- function(x, name, lower = -Inf, upper = Inf, null_ok = TRUE) {
  if (is.null(x) && null_ok) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower ||
      x > upper) {
    stop("`", name, "` must be a single number in [", lower, ", ", upper,
         "].", call. = FALSE)
  }
  x
}

# Encodes the covariates as the design matrix used in both estimation steps:
# numeric covariates are standardised, factors are dummy coded with treatment
# contrasts (the first level is the reference). Returns the design matrix and
# the information needed to encode new subjects in the same way.
encode_covariates <- function(covariates, n) {
  if (is.null(covariates)) {
    return(NULL)
  }
  if (is.matrix(covariates)) {
    covariates <- as.data.frame(covariates)
  }
  if (!is.data.frame(covariates)) {
    stop("`covariates` must be a data frame.", call. = FALSE)
  }
  if (nrow(covariates) != n) {
    stop("`covariates` must have one row per row of `x`.", call. = FALSE)
  }
  if (ncol(covariates) == 0L) {
    return(NULL)
  }
  if (anyNA(covariates)) {
    stop("`covariates` contains missing values.", call. = FALSE)
  }
  if (is.null(names(covariates)) || anyDuplicated(names(covariates))) {
    stop("`covariates` must have unique column names.", call. = FALSE)
  }
  numeric_vars <- character(0)
  scaling <- list()
  levels_list <- list()
  for (v in names(covariates)) {
    col <- covariates[[v]]
    if (is.numeric(col)) {
      s <- sd(col)
      if (s == 0) {
        stop("Covariate `", v, "` is constant.", call. = FALSE)
      }
      numeric_vars <- c(numeric_vars, v)
      scaling[[v]] <- list(mean = mean(col), sd = s)
    } else if (is.factor(col) || is.character(col) || is.logical(col)) {
      col <- factor(col)
      if (nlevels(col) < 2L) {
        stop("Covariate `", v, "` has a single level.", call. = FALSE)
      }
      levels_list[[v]] <- levels(col)
    } else {
      stop("Covariate `", v, "` must be numeric, factor, character or ",
           "logical.", call. = FALSE)
    }
  }
  info <- list(
    variables = names(covariates),
    scaling = scaling,
    levels = levels_list
  )
  U <- covariate_design(covariates, info)
  info$design_names <- colnames(U)
  list(U = U, info = info)
}

# Applies the stored encoding to a data frame of covariates.
covariate_design <- function(covariates, info) {
  missing_vars <- setdiff(info$variables, names(covariates))
  if (length(missing_vars) > 0L) {
    stop("Missing covariate(s): ", paste(missing_vars, collapse = ", "), ".",
         call. = FALSE)
  }
  df <- covariates[, info$variables, drop = FALSE]
  if (anyNA(df)) {
    stop("Covariates contain missing values.", call. = FALSE)
  }
  design_names <- character(0)
  for (v in info$variables) {
    if (v %in% names(info$scaling)) {
      if (!is.numeric(df[[v]])) {
        stop("Covariate `", v, "` must be numeric.", call. = FALSE)
      }
      df[[v]] <- (df[[v]] - info$scaling[[v]]$mean) / info$scaling[[v]]$sd
      design_names <- c(design_names, v)
    } else {
      levs <- info$levels[[v]]
      values <- as.character(df[[v]])
      unknown <- setdiff(unique(values), levs)
      if (length(unknown) > 0L) {
        stop("Covariate `", v, "` has level(s) not seen during estimation: ",
             paste(unknown, collapse = ", "), ".", call. = FALSE)
      }
      df[[v]] <- factor(values, levels = levs)
      design_names <- c(design_names, paste0(v, "_", levs[-1]))
    }
  }
  factor_vars <- names(info$levels)
  contrasts <- NULL
  if (length(factor_vars) > 0L) {
    contrasts <- stats::setNames(
      rep(list("contr.treatment"), length(factor_vars)), factor_vars)
  }
  U <- model.matrix(~ ., data = df, contrasts.arg = contrasts)
  U <- U[, -1, drop = FALSE]
  colnames(U) <- design_names
  rownames(U) <- NULL
  U
}

# Symmetrises a matrix of node-wise estimates. The choice between the (j, k)
# and (k, j) estimates is based on the magnitudes in `ref` (the regression
# coefficients): "OR" keeps the larger one, "AND" the smaller one.
symmetrize <- function(m, method, ref = m) {
  if (method == "OR") {
    out <- ifelse(abs(ref) >= abs(t(ref)), m, t(m))
  } else {
    out <- ifelse(abs(ref) < abs(t(ref)), m, t(m))
  }
  dimnames(out) <- dimnames(m)
  out
}

vmessage <- function(verbose, ...) {
  if (verbose) {
    message(...)
  }
}
