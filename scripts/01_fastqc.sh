#!/bin/bash

module load FastQC_0.11.9

# Run FastQC 
for i in $(ls *.fastq.gz) ; do 
echo $i 
srun fastqc -o ${i} -t 2 
done
