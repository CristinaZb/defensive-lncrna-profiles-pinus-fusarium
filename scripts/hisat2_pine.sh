#!/bin/bash

# Define the files
files_1=(*_1p.fastq.gz)
files_2=(*_2p.fastq.gz)
files_3=(*_1u.fastq.gz)

name1=${files_1[$SLURM_ARRAY_TASK_ID-1]}
name2=${files_2[$SLURM_ARRAY_TASK_ID-1]}
name3=${files_3[$SLURM_ARRAY_TASK_ID-1]}

echo $SLURM_ARRAY_TASK_ID
file=$(ls | cut -d_ -f1 | uniq | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load haswell/HISAT2_2.0.0

# Run Hisat2

hisat2 -p $SLURM_NTASKS --rg-id=${file} --rg SM:${file} --rg PL:ILLUMINA \
-x $PITA_PATH/genome_pita_index --dta --rna-strandness RF \
-1 ${name1} -2 ${name2} -U ${name3} \
-S $MAPPED_PATH/${file}.sam
