#!/bin/bash
# Run richness GAM on HPC.

#SBATCH --job-name=bam_ktest
#SBATCH --mem=512G
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --output=bam_ktest_%j.log

module load r/4.4.1

# Run from project root
Rscript script/08_hpc_install.R
Rscript script/08_richness_gam_hpc.R
