# Pinus–Fusarium lncRNA (dual RNA-seq)

Reproducible dual RNA-seq analysis of defensive lncRNAs in *Pinus* challenged by *Fusarium circinatum*.  
**Scope:** code-only repository (no raw data).

## Overview
This repo contains the scripts and minimal configuration to reproduce:


## What’s here
- `scripts/` — Bash/R scripts used in the paper (cleaned and documented).
  
> **Note:** Raw FASTQs, genome indices, and large intermediates are intentionally **not** tracked.

## Quick start
1. Clone or download this repo.
2. Check `scripts/` headers for dependencies (R ≥ 4.3; packages listed at the top of each script).
3. Edit `config/` files (paths/contrasts) if you want to rerun parts of the analysis.
4. Run scripts in the numbered order (see filenames).

## Reproducibility
- Versions are tagged; releases will be archived in Zenodo to obtain a DOI.
- If you find any non-determinism, open an issue.

## How to cite
If you use this code, please cite this repository and the related manuscript (when available).  
*(A `CITATION.cff` file will be added for automated citation metadata.)*

## Pipeline (at a glance)

## **Quality control** → [`scripts/fastqc.sh`](scripts/fastqc.sh)

We assessed read quality with FastQC on all raw FASTQ files (see scripts/fastqc.sh). 

## **Trimming** → [`scripts/trimming.sh`](scripts/trimming.sh)  

Reads were trimmed with Trimmomatic 0.38 using Illumina adapter removal and light head cropping (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:10). Post-trim QC was re-run to confirm improvement. All auxiliary *_2u.fastq.gz (unpaired R2) files were found empty and were removed. See scripts/trimming.sh



5) **Alignment (host / pathogen)** → [`scripts/hisat2.sh`](scripts/hisat2.sh)
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
