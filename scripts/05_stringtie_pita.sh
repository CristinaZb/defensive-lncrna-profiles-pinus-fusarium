#!/bin/bash

file=$(ls *.bam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load haswell/StringTie_2.1.4

stringtie ${file}.dedup.bam -G $GTF/Pita.2_01.gtf -p $SLURM_NTASKS -l ${file}_Pinus --rf -o $OUT_PATH/${file}_transcripts.gtf
