#!/usr/bin/env Rscript
# Fits one node of a prCGGM model; called by each task of the SLURM array.
# Usage: Rscript prcggm_node_job.R <input_file> <node_index> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript prcggm_node_job.R <input_file> <node_index> <output_dir>")
}
prCGGM::prcggm_node(input = args[1], node = as.integer(args[2]),
                    output_dir = args[3])
