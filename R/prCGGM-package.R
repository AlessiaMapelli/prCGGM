#' prCGGM: Prior-Informed Conditional Gaussian Graphical Models
#'
#' Estimation of conditional Gaussian graphical models in which both the mean
#' and the precision matrix depend linearly on subject-level covariates, with
#' prior knowledge on edges entering the population-level network through a
#' weighted penalty.
#'
#' The main entry points are:
#' * [prcggm()] to fit the model;
#' * [predict.prcggm()] to obtain personalised networks for new individuals;
#' * [summary.prcggm()] and [plot.prcggm()] to inspect a fitted model;
#' * [prcggm_slurm_prepare()], [prcggm_node()] and [prcggm_slurm_collect()]
#'   to run the node-wise estimation as 'SLURM' array jobs.
#'
#' @references
#' Mapelli, A., Massi, M. C., Cuccuru, G., Di Angelantonio, E. and Ieva, F.
#' (2026). Prior-informed conditional Gaussian graphical models: an
#' application to protein interaction network reconstruction.
#' arXiv:2606.31805. \doi{10.48550/arXiv.2606.31805}
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats coef cor lm model.matrix predict quantile sd
NULL
