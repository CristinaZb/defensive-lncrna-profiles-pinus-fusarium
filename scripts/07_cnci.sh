#!/bin/bash

module load CNCI_06Mar2018

# Run
python CNCI.py -f candidate_${spp}-lncRNA.fa -o cnci_${spp} -m pl -p $SLURM_NTASKS
