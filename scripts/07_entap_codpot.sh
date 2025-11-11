#!/bin/bash

# Define paths and environment variables
module load EnTAP_1.9.2_gcc8.2.0
module load jdk_11.0.2

# Run EnTAP
EnTAP --runP \
-i candidate_${spp}-lncRNA.fa \
-d EnTAP/1.9.2/databases/diamond/nr.dmnd \
--level 0 \
--ontology 0 \
--threads $SLURM_NTASKS \
-c fungi \
--taxon pinus
