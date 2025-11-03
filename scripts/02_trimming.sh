#!/bin/bash

files_1=(*_1.fastq.gz)
files_2=(*_2.fastq.gz)

name1=${files_1[$SLURM_ARRAY_TASK_ID-1]}
name2=${files_2[$SLURM_ARRAY_TASK_ID-1]}

echo $SLURM_ARRAY_TASK_ID
file=$(ls | cut -f1 -d_ | uniq | tail -n +$SLURM_ARRAY_TASK_ID | head -n 1)

# Run Trimmomatic

java -jar $TRIMMOMATIC_PATH/trimmomatic-0.38.jar PE -threads $SLURM_NTASKS ${name1} ${name2} \
$TRIMMOMATIC_OUTPUT/${file}_1p.fastq.gz $TRIMMOMATIC_OUTPUT/${file}_1u.fastq.gz \
$TRIMMOMATIC_OUTPUT/${file}_2p.fastq.gz $TRIMMOMATIC_OUTPUT/${file}_2u.fastq.gz \
ILLUMINACLIP:$TRIMMOMATIC_PATH/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:10
