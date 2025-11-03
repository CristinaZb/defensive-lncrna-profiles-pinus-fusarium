#!/bin/bash

file=$(ls *.bam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load haswell/samtools_1.7

samtools fixmate -m -@ $SLURM_NTASKS ${file}.name.bam ${file}.fix.bam
