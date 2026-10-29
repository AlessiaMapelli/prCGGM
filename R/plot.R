#' Plot a network as a heatmap
#'
#' Draws an adjacency matrix (for example a component of a fitted model or a
#' personalised network from [predict.prcggm()]) as a heatmap. The
#' \pkg{corrplot} package is used if installed, otherwise
#' [stats::heatmap()].
#'
#' @param network Square numeric matrix.
#' @param labels Optional node labels. By default the row or column names of
#'   `network`, or `Node_1`, `Node_2`, and so on.
#' @param threshold Entries with absolute value below `threshold` are set to
#'   zero before plotting. Default 0.
#' @param main Plot title.
#'
#' @return Called for its side effect; returns the plotted matrix invisibly.
#'
#' @examples
#' net <- matrix(c(0, 0.4, 0, 0.4, 0, -0.3, 0, -0.3, 0), 3, 3)
#' plot_network(net, labels = c("A", "B", "C"))
#' @export
plot_network <- function(network, labels = NULL, threshold = 0,
                         main = "Network") {
  network <- as.matrix(network)
  if (!is.numeric(network) || nrow(network) != ncol(network)) {
    stop("`network` must be a square numeric matrix.", call. = FALSE)
  }
  if (is.null(labels)) {
    labels <- rownames(network)
    if (is.null(labels)) labels <- colnames(network)
    if (is.null(labels)) labels <- paste0("Node_", seq_len(nrow(network)))
  }
  if (length(labels) != nrow(network)) {
    stop("`labels` must have one entry per node.", call. = FALSE)
  }
  network[abs(network) < threshold] <- 0
  dimnames(network) <- list(labels, labels)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  if (requireNamespace("corrplot", quietly = TRUE)) {
    lim <- max(abs(network))
    if (lim == 0) lim <- 1
    # corrplot warns when the colour range is wider than the data range,
    # which is the intended behaviour of a symmetric scale
    withCallingHandlers(
      corrplot::corrplot(network, method = "color", is.corr = FALSE,
                         col.lim = c(-lim, lim),
                         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
                         title = main, mar = c(0, 0, 2, 0), cl.pos = "r"),
      warning = function(w) {
        if (grepl("col.lim interval too wide", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      })
  } else {
    stats::heatmap(network, main = main, Rowv = NA, Colv = NA,
                   scale = "none", cexRow = 0.8, cexCol = 0.8,
                   margins = c(10, 10),
                   col = grDevices::colorRampPalette(
                     c("blue", "white", "red"))(100))
  }
  invisible(network)
}

#' Plot the difference between two networks
#'
#' Draws `network1 - network2` with [plot_network()], for example to compare
#' the personalised networks of two individuals.
#'
#' @param network1,network2 Square numeric matrices of the same size.
#' @param labels Optional node labels, as in [plot_network()].
#' @param main Plot title.
#'
#' @return Called for its side effect; returns the difference matrix
#'   invisibly.
#'
#' @examples
#' a <- matrix(c(0, 0.4, 0.4, 0), 2, 2)
#' b <- matrix(c(0, 0.1, 0.1, 0), 2, 2)
#' plot_network_difference(a, b)
#' @export
plot_network_difference <- function(network1, network2, labels = NULL,
                                    main = "Network difference") {
  check_same_dim(network1, network2)
  diff_network <- as.matrix(network1) - as.matrix(network2)
  if (is.null(labels)) {
    labels <- rownames(network1)
    if (is.null(labels)) labels <- colnames(network1)
  }
  plot_network(diff_network, labels = labels, main = main)
}

#' Compare two networks
#'
#' Computes a summary measure of the difference between two adjacency
#' matrices.
#'
#' @param network1,network2 Square numeric matrices of the same size.
#' @param method `"difference"` (default) for the mean absolute difference,
#'   `"correlation"` for the Pearson correlation of the upper-triangular
#'   entries, or `"frobenius"` for the Frobenius norm of the difference.
#'
#' @return A single number.
#'
#' @examples
#' a <- matrix(c(0, 0.4, 0, 0.4, 0, -0.3, 0, -0.3, 0), 3, 3)
#' b <- matrix(c(0, 0.1, 0, 0.1, 0, -0.1, 0, -0.1, 0), 3, 3)
#' compare_networks(a, b, method = "frobenius")
#' @export
compare_networks <- function(network1, network2,
                             method = c("difference", "correlation",
                                        "frobenius")) {
  method <- match.arg(method)
  check_same_dim(network1, network2)
  network1 <- as.matrix(network1)
  network2 <- as.matrix(network2)
  switch(method,
    difference = mean(abs(network1 - network2)),
    correlation = {
      upper <- upper.tri(network1)
      cor(network1[upper], network2[upper])
    },
    frobenius = norm(network1 - network2, "F")
  )
}

check_same_dim <- function(network1, network2) {
  if (!identical(dim(as.matrix(network1)), dim(as.matrix(network2)))) {
    stop("Networks must have the same dimensions.", call. = FALSE)
  }
}
