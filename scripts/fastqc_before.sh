PATH_RAW="/dual_lncrna/raw_data" 
PATH_OUT="/dual_lncrna/fastqc_before" 
module load haswell/FastQC_0.11.9 
# Run FastQC 
for i in $(ls $PATH_RAW/*.fastq.gz) ; do 
echo $i 
srun fastqc -o $PATH_OUT ${i} -t 2 
done
