#!/bin/bash

ls -1 *.gtf >> sample_assembly_gtf_list.txt

module load StringTie_2.1.4

stringtie --merge -p $SLURM_NTASKS \
-G Pita.2_01.gtf \
-F 0.5 sample_assembly_gtf_list.txt \
-o ${spp}_transcripts.gtf \
-l ${spp}
