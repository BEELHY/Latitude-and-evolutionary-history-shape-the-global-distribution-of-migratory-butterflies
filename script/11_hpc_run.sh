#!/bin/bash
# Run richness GAM on HPC.

module load r/4.4.1

Rscript 11_hpc_install.R
Rscript 11_richness_gam_hpc.R
