#!/bin/bash

file=$(ls | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load haswell/samtools_1.7

samtools sort -n -@ $SLURM_NTASKS -o ${file}.name.bam ${file}.bam
