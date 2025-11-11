#!/bin/bash

# We activate the conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate env

# Run
cpat.py -g candidate_${spp}-lncRNA.fa -d ${spp}.logit.RData -x ${spp}_Hexamer.tsv -o ${spp}_cpat

# We deactivate our environment
conda deactivate
