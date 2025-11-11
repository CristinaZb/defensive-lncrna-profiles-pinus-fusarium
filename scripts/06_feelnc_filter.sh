#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate feelnc_env

FEELnc_filter.pl -p $SLURM_NTASKS -i unknown_${spp}.annotated.gtf -a Pita.2_01.gtf --monoex=-1 --biex=1 > candidate_${spp}-lncRNA.gtf

# deactivate the environment
conda deactivate
