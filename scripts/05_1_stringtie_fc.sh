#!/bin/bash

file=$(ls *.bam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load StringTie_2.1.4

stringtie ${file}.dedup.bam -p $SLURM_NTASKS -l ${file}_FCRIC --rf -o $OUT_PATH/${file}_transcripts.gtf
