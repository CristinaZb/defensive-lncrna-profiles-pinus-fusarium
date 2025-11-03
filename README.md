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

## **Alignment (host / pathogen)**

Trimmed reads were aligned with HISAT2 (`--dta`, `stranded RF`) against the host genome (_Pinus taeda_ Ptaeda2.0) and, separately, against _Fusarium circinatum_ FC072V (see [`scripts/hisat2_pine.sh`](scripts/hisat2_pine.sh) and [`scripts/hisat2_fc.sh`](scripts/hisat2_fc.sh); reference details in `data/DATA.md`). 

Typical host mapping rates were ~79–81% for _P. radiata_ and ~44–49% for _P. pinea_ to Ptaeda2.0, whereas pathogen alignments were expectedly low (<~4%) given plant-dominant libraries.

## **SAM processing**

Alignment outputs (SAM) were converted to BAM, coordinate-sorted, and indexed with samtools. Basic alignment QC was recorded with samtools flagstat. Intermediate SAM files were removed to save space. Host and pathogen BAMs are kept separately (one per sample).

   - SAM to BAM → [`scripts/samtools_stgtie.sh`](scripts/samtools_stgtie.sh) 
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
