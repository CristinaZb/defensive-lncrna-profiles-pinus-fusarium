#!/bin/bash

# We activate the conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate env

# Run
PLEK.py -fasta candidate_${spp}-lncRNA.fa -out ${spp}_plek_predicted -thread $SLURM_NTASKS
