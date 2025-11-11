#!/bin/bash

# Define files
files_1=(*_1p.fastq.gz)
files_2=(*_2p.fastq.gz)
files_3=(*_1u.fastq.gz)

name1=${files_1[$SLURM_ARRAY_TASK_ID-1]}
name2=${files_2[$SLURM_ARRAY_TASK_ID-1]}
name3=${files_3[$SLURM_ARRAY_TASK_ID-1]}

echo $SLURM_ARRAY_TASK_ID
file=$(ls | cut -f1 -d_ | uniq | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

module load HISAT2_2.0.0

# Run Hisat2

srun -n 1 hisat2 -p $SLURM_NTASKS --rg-id=${file} --rg SM:${file} --rg PL:ILLUMINA \
-x $FCIRC_PATH/genome_fc072v_index --dta --rna-strandness RF \
-1 ${name1} -2 ${name2} -U ${name3} \
-S ${file}-fc072v.sam
