#!/bin/bash

file=$(ls *.bam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load haswell/samtools_1.7

samtools markdup -r -s -@ $SLURM_NTASKS ${file}.sorted.bam ${file}.dedup.bam
