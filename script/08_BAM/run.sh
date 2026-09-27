#!/bin/bash
# Run spatial model on HPC.

module load r/4.4.1

Rscript install.R
Rscript 1500test.R
