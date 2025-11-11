#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate feelnc_env

srun FEELnc_codpot.pl -p $SLURM_NTASKS 
-i candidate_${spp}-lncRNA.gtf 
-a $KNOWN/Pita.gtf 
-g $GENOME/Pita.fa 
--kmer="1,2,3,6,9,12" 
--spethres=0.95,0.95 
--mode=shuffle 
--outname="${spp}"

# deactivate conda environment
conda deactivate
