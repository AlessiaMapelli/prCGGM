#' Example data for prCGGM
#'
#' A small simulated data set with 1000 subjects, 10 proteins and one binary
#' covariate, generated from a conditional Gaussian graphical model with a
#' scale-free baseline network and an Erdos-Renyi covariate-specific network,
#' as in the simulation study of Mapelli et al. (2026). The prior weight
#' matrix is "perfect": it equals 1 on the true baseline edges and 0
#' elsewhere.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{`x`}{Numeric matrix with 1000 rows (subjects) and 10 columns
#'     (`Prot1` to `Prot10`) of expression values.}
#'   \item{`covariates`}{Data frame with one factor column, `x1`, with
#'     levels `"0"` (reference) and `"1"`.}
#'   \item{`prior`}{10 x 10 prior weight matrix with entries in \{0, 1\}.}
#' }
#' @source Simulated; see `data-raw/prcggm_example.R` in the package sources.
#' @examples
#' data(prcggm_example)
#' str(prcggm_example)
"prcggm_example"
