#!/bin/bash

file=$(ls *.sam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load samtools_1.7

samtools view -@ $SLURM_NTASKS -o ${file}.bam ${file}.sam
