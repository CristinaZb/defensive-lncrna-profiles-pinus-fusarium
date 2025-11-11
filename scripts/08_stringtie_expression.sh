#!/bin/bash

file=$(ls *.bam | cut -f1 -d. | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

echo ${file}
echo $SLURM_ARRAY_TASK_ID

module load StringTie_2.1.4

# Run StringTie
stringtie ${file}.dedup.bam -e -p $SLURM_NTASKS -G ${spp}_transcripts-mod.gtf -o ${file}.gtf -A ${file}_gene_abundance.out
