# Pinus–Fusarium lncRNA (dual RNA-seq)

Reproducible dual RNA-seq analysis of defensive lncRNAs in *Pinus* challenged by *Fusarium circinatum*.  
**Scope:** code-only repository (no raw data).

## Overview
This repo contains the scripts and minimal configuration to reproduce:

## How to cite
If you use this code, please cite this repository and the related manuscript (when available).  
*(A `CITATION.cff` file will be added for automated citation metadata.)*

# Bioinformatics pipeline

## **Quality control**

We assessed read quality with FastQC on all raw FASTQ files (see [`scripts/fastqc.sh`](scripts/fastqc.sh)). 

## **Trimming**  

Reads were trimmed with Trimmomatic 0.38 using Illumina adapter removal and light head cropping (`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:10`). Post-trim QC was re-run to confirm improvement ([`scripts/fastqc.sh`](scripts/fastqc.sh)). All auxiliary *_2u.fastq.gz (unpaired 2) files were found empty and were removed. See [`scripts/trimming.sh`](scripts/trimming.sh)

Trimmed pairs were then split by species into pr_ (_P. radiata_) and pp_ (_P. pinea_) sets.

## **Alignment (host / pathogen)** → [`scripts/hisat2.sh`](scripts/hisat2.sh)




6) **SAM processing**
   - SAM to BAM → [`scripts/26_counts_featurecounts.sh`](scripts/26_counts_featurecounts.sh) 
   - Sorting BAM by name
   - Fixing BAM by removing duplicates
   - Sorting by coordinates
   - Removing duplicates
7) **Pine assembly**
8) **Pathogen assembly**
9) **Pine comparation** 
10) **Identification of long non-coding RNAs**  
   - Step 1:  → [`scripts/26_counts_featurecounts.sh`](scripts/26_counts_featurecounts.sh)  
   - Step 2:
  
11) 
12) **Differential expression (DESeq2)** → [`scripts/30_deseq2.R`](scripts/30_deseq2.R)  
13) **Co-expression (WGCNA, optional)** → [`scripts/40_wgcna.R`](scripts/40_wgcna.R)  
14) **lncRNA identification (notes + commands)** → [`scripts/50_lncrna_notes.md`](scripts/50_lncrna_notes.md)  
15) **Functional annotation (optional)** → [`scripts/60_annotation_notes.md`](scripts/60_annotation_notes.md)




## License
MIT — see `LICENSE`.

## Contact
Maintainer: Cristina Zamora Ballesteros (cristinazamoraballesteros@gmail.com) — issues and PRs welcome.
