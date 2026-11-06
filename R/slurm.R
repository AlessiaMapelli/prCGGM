#' Run the estimation as SLURM array jobs
#'
#' For large networks the node-wise regressions can be distributed over a
#' computing cluster, one 'SLURM' array task per node, in three steps:
#' 1. `prcggm_slurm_prepare()` validates the inputs, runs the mean step and
#'    saves everything the node-wise regressions need in `dir`. It also copies
#'    template job scripts to `dir` and returns the `sbatch` command to
#'    submit.
#' 2. Each array task calls `prcggm_node()` (through the copied
#'    `prcggm_node_job.R` script) and saves the result of one node in `dir`.
#' 3. Once all tasks have finished, `prcggm_slurm_collect()` reads the
#'    node results and returns the same `"prcggm"` object as [prcggm()].
#'
#' Edit the copied `prcggm_node.sbatch` file to match your cluster (partition,
#' memory, time, and how R is loaded) before submitting. The template files
#' are also available in `system.file("slurm", package = "prCGGM")`.
#'
#' @inheritParams prcggm
#' @param dir Directory where inputs, job scripts and node results are
#'   written. It is created if needed. There is no default: choose a location
#'   on a file system shared by the cluster nodes.
#' @param copy_templates Logical; copy the template job scripts to `dir`.
#'   Default `TRUE`.
#'
#' @return `prcggm_slurm_prepare()` returns, invisibly, a list with the path
#'   of the saved input file (`input`), the number of nodes (`p`) and the
#'   `sbatch` command to run (`command`), which is also printed as a message.
#'
#' @examples
#' \dontrun{
#' # On the cluster login node
#' job <- prcggm_slurm_prepare(expr, covariates = covs, prior = W,
#'                             dir = "/scratch/my_project/prcggm")
#' # Submit job$command from a shell, e.g.
#' #   sbatch --array=1-366 /scratch/my_project/prcggm/prcggm_node.sbatch ...
#' # When all array tasks have completed:
#' fit <- prcggm_slurm_collect("/scratch/my_project/prcggm")
#' }
#'
#' # The same three steps run locally, without SLURM:
#' data(prcggm_example)
#' dir <- file.path(tempdir(), "prcggm_slurm")
#' set.seed(1)
#' job <- prcggm_slurm_prepare(prcggm_example$x,
#'                             covariates = prcggm_example$covariates,
#'                             prior = prcggm_example$prior,
#'                             tune = FALSE, dir = dir,
#'                             copy_templates = FALSE)
#' for (i in seq_len(job$p)) prcggm_node(job$input, i)
#' fit <- prcggm_slurm_collect(dir)
#' unlink(dir, recursive = TRUE)
#' @export
prcggm_slurm_prepare <- function(x,
                                 covariates = NULL,
                                 prior = NULL,
                                 dir,
                                 screening = FALSE,
                                 gamma = NULL,
                                 screening_matrix = NULL,
                                 mean_estimation = TRUE,
                                 lambda_mean = NULL,
                                 lambda_mean_type = c("min", "1se"),
                                 lambda_prec = NULL,
                                 selection = c("ebic_relaxed", "bic_relaxed",
                                               "bic", "ebic", "cv_min",
                                               "cv_1se", "stability_bic",
                                               "stability_cv"),
                                 refit = TRUE,
                                 selection_control = list(),
                                 tune = TRUE,
                                 alpha_grid = c(0.75, 0.5, 0.9, 0.99),
                                 theta_grid = c(1, 0.8, 1.1, 1.3),
                                 random_search = FALSE,
                                 p_random = 0.5,
                                 nfolds = 5,
                                 symm_method = c("OR", "AND"),
                                 copy_templates = TRUE,
                                 verbose = FALSE,
                                 lambda_prec_type = NULL) {
  if (missing(dir) || !is.character(dir) || length(dir) != 1L) {
    stop("`dir` must be the path of a directory.", call. = FALSE)
  }
  call <- match.call()
  selection <- resolve_selection(selection, lambda_prec_type,
                                 missing(selection))
  prep <- prepare_fit(
    x = x, covariates = covariates, prior = prior, screening = screening,
    gamma = gamma, screening_matrix = screening_matrix,
    mean_estimation = mean_estimation, lambda_mean = lambda_mean,
    lambda_mean_type = match.arg(lambda_mean_type),
    lambda_prec = lambda_prec, selection = selection, refit = refit,
    selection_control = selection_control,
    tune = tune, alpha_grid = alpha_grid, theta_grid = theta_grid,
    random_search = random_search, p_random = p_random, nfolds = nfolds,
    symm_method = match.arg(symm_method), verbose = verbose)
  prep$call <- call
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir)
  input <- file.path(dir, "prcggm_input.rds")
  saveRDS(prep, input)
  script <- file.path(dir, "prcggm_node.sbatch")
  if (check_flag(copy_templates, "copy_templates")) {
    templates <- system.file("slurm", package = "prCGGM")
    file.copy(file.path(templates, c("prcggm_node.sbatch",
                                     "prcggm_node_job.R")),
              dir, overwrite = TRUE)
    # SLURM opens the log files of the template before the job starts
    dir.create(file.path(dir, "logs"), showWarnings = FALSE)
  }
  command <- sprintf("cd %s && sbatch --array=1-%d %s %s %s",
                     shQuote(dir), prep$p, shQuote(script), shQuote(input),
                     shQuote(dir))
  message("Inputs saved to ", input, "\nSubmit the array job with:\n  ",
          command, "\nthen run prcggm_slurm_collect(", deparse(dir),
          ") once all tasks have completed.")
  invisible(list(input = input, p = prep$p, command = command))
}

#' @rdname prcggm_slurm_prepare
#' @param input Path of the input file written by `prcggm_slurm_prepare()`.
#' @param node Index of the node (feature) to estimate, from 1 to `p`.
#' @param output_dir Directory where the node result is saved. Defaults to
#'   the directory of `input`.
#' @return `prcggm_node()` returns, invisibly, the path of the saved node
#'   result.
#' @export
prcggm_node <- function(input, node, output_dir = dirname(input)) {
  prep <- readRDS(input)
  node <- check_number(node, "node", lower = 1, upper = prep$p,
                       null_ok = FALSE)
  result <- fit_node(as.integer(node), prep)
  file <- file.path(output_dir, sprintf("prcggm_node_%d.rds", node))
  saveRDS(result, file)
  invisible(file)
}

#' @rdname prcggm_slurm_prepare
#' @return `prcggm_slurm_collect()` returns a `"prcggm"` object, as
#'   [prcggm()].
#' @export
prcggm_slurm_collect <- function(dir, output_dir = dir) {
  input <- file.path(dir, "prcggm_input.rds")
  if (!file.exists(input)) {
    stop("No prcggm_input.rds file found in `dir`.", call. = FALSE)
  }
  prep <- readRDS(input)
  files <- file.path(output_dir,
                     sprintf("prcggm_node_%d.rds", seq_len(prep$p)))
  missing_nodes <- which(!file.exists(files))
  if (length(missing_nodes) > 0L) {
    stop("Results missing for ", length(missing_nodes), " node(s): ",
         paste(utils::head(missing_nodes, 20), collapse = ", "),
         if (length(missing_nodes) > 20) ", ...", ".", call. = FALSE)
  }
  fit <- assemble_fit(prep, lapply(files, readRDS))
  fit$call <- prep$call
  fit
}
