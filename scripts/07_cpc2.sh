#!/bin/bash 

module load haswell/CPC2_1.0.1 

# Run 
CPC2.py -i $CANDIDATE/candidate_${spp}-lncRNA.fa -o CPC2_${spp}-lncRNA
