#!/bin/bash
#SBATCH --job-name=prcggm
#SBATCH --output=logs/prcggm_node_%A_%a.out
#SBATCH --error=logs/prcggm_node_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=5GB
##SBATCH --partition=<your_partition>
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=<your_email>

# SLURM array job template for prCGGM: each array task fits one node.
# Usage (as printed by prCGGM::prcggm_slurm_prepare()):
#   cd <dir> && sbatch --array=1-<p> prcggm_node.sbatch <input_file> <output_dir>
# Adjust the #SBATCH options above to your cluster before submitting.

# The logs/ directory must exist before submission (prcggm_slurm_prepare()
# creates it), because SLURM opens the output files before this script runs.

INPUT_FILE=$1
OUTPUT_DIR=$2
NODE_INDEX=$SLURM_ARRAY_TASK_ID

# Load R as required by your cluster, for example:
# module load R

# prCGGM must be installed in the R used by the jobs. If it is installed in
# a personal library, point R to it here (Rscript --vanilla ignores
# ~/.Renviron), for example:
# export R_LIBS=$HOME/R/prcggm_lib

# Avoids "Setting LC_CTYPE failed" warnings on nodes without a locale:
# export LC_ALL=C.UTF-8

echo "prCGGM: fitting node $NODE_INDEX"
Rscript --vanilla prcggm_node_job.R "$INPUT_FILE" "$NODE_INDEX" "$OUTPUT_DIR"
echo "prCGGM: node $NODE_INDEX completed"
