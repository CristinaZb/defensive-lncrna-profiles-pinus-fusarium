#!/bin/bash

file=$(ls *.bam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load samtools_1.7

samtools sort -@ $SLURM_NTASKS -o ${file}.sorted.bam ${file}.fix.bam
