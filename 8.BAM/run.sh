#!/bin/bash
#SBATCH --job-name=bam_ktest
#SBATCH --mem=512G
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --output=bam_ktest_%j.log

module load r/4.4.1

Rscript install.R
Rscript 1500test.R
