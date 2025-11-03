# Pinus–Fusarium lncRNA (dual RNA-seq)

Reproducible dual RNA-seq analysis of defensive lncRNAs in *Pinus* challenged by *Fusarium circinatum*.  
**Scope:** code-only repository (no raw data).

## Data overview

This repository is **code-only**. Raw sequencing files and large intermediates are **not** distributed here due to size and licensing constraints. The goal of this document is to record provenance and make the analysis reproducible.

## Raw reads (not included)

- Type: paired-end RNA-Seq (dual RNA-seq: *Pinus* host + *Fusarium circinatum* pathogen).
- Number of original files: `<e.g., 32 paired libraries>`; total size: `<e.g., ~135 GB>`.
- Example filenames: `SampleX_1.fastq.gz`, `SampleX_2.fastq.gz`.
- Source: Raw reads have been deposited in the NCBI SRA Database under accession numbers SRR13737940-53 (BioProject PRJNA702546).

## References (not included)

- **Host genome/annotation**
  - Organism: *Pinus taeda* (loblolly pine)
  - Assembly name/version: **Pita v2.01**
  - NCBI Assembly record: **GCA_000404065.3**
  - Assembly level: Scaffold; submitter: **TreeGenes Database**; TreeGenes FTP: `treegenesdb.org/FTP/Genomes/Pita/v2.01/`
  - Genome FASTA: Pita.2_01.fa.gz
  - Gene annotation (GTF): Pita.2_01.gtf.gz 
 
- **Pathogen genome/annotation**
  - Name: Fusarium circinatum strain FC072V
  - Organism: *Fusarium circinatum* (TaxID: 48490)
  - BioProject: **PRJNA716280** — Genome sequencing and assembly; scope: monoisolate; registered 2021-04-27
  - Assembly accession: **GCA_018163655.1** (assembly level: Contig)
  - WGS master: **JAGGEA000000000**
  - BioSample: **SAMN18415966**
  - Source: NCBI BioProject / Assembly records

## How to cite
If you use this code, please cite this repository and the related manuscript (when available).  
*(A `CITATION.cff` file will be added for automated citation metadata.)*

# Bioinformatics pipeline

## **Quality control**

We assessed read quality with FastQC on all raw FASTQ files (see [`scripts/01_fastqc.sh`](scripts/01_fastqc.sh)). 

## **Trimming**  

Reads were trimmed with Trimmomatic 0.38 using Illumina adapter removal and light head cropping (`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:10`). Post-trim QC was re-run to confirm improvement ([`scripts/01_fastqc.sh`](scripts/01_fastqc.sh)). All auxiliary *_2u.fastq.gz (unpaired 2) files were found empty and were removed. See [`scripts/02_trimming.sh`](scripts/02_trimming.sh)

Trimmed pairs were then split by species into pr_ (_P. radiata_) and pp_ (_P. pinea_) sets.

## **Alignment (host / pathogen)**

Trimmed reads were aligned with HISAT2 (`--dta`, `stranded RF`) against the host genome (_Pinus taeda_ Ptaeda2.0) and, separately, against _Fusarium circinatum_ FC072V (see [`scripts/03_hisat2_pine.sh`](scripts/03_hisat2_pine.sh) and [`scripts/03_hisat2_fc.sh`](scripts/03_hisat2_fc.sh); reference details in `data/DATA.md`). 

Typical host mapping rates were ~79–81% for _P. radiata_ and ~44–49% for _P. pinea_ to Ptaeda2.0, whereas pathogen alignments were expectedly low (<~4%) given plant-dominant libraries.

## **SAM processing**

Alignment outputs (SAM) were converted to BAM, coordinate-sorted, and indexed with samtools. Intermediate SAM files were removed to save space. Host and pathogen BAMs are kept separately (one per sample).

   - SAM to BAM → [`scripts/04_1_samtools_stgtie.sh`](scripts/04_1_samtools_stgtie.sh) 
   - Sorting BAM by name → [`scripts/04_2_samtools_namesort.sh`](scripts/04_2_samtools_namesort.sh)
   - Fixing BAM by removing duplicates → [`scripts/04_3_samtools_fixmate.sh`](scripts/04_3_samtools_fixmate.sh)
   - Sorting by coordinates → [`scripts/04_4_samtools_sort.sh`](scripts/04_4_samtools_sort.sh)
   - Removing duplicates → [`scripts/04_5_samtools_rmdup.sh`](scripts/04_5_samtools_rmdup.sh)

##  **Pine assembly**

Per-sample BAMs were assembled with StringTie in reference-guided mode for the pine species (`Pita.2_01.gtf.gz`; [`scripts/05_stringtie_pita.sh`](scripts/05_stringtie_pita.sh)). Per-sample GTFs were then merged into a non-redundant consensus using `stringtie --merge` [`scripts/05_stringtie_merge.sh`](scripts/05_stringtie_merge.sh).

Number of assembled transcripts in each GTF:

```bash
awk '$3=="transcript"' ${file}_transcripts.gtf | wc -l
```

In order to get a fasta file with the sequences of the assembled transcripts of each organism, extract transcripts from the non-redundant GTF and convert to FASTA file using `gffread` and the reference genome: [`scripts/05_gffread_pine.sh`](scripts/05_gffread_pine.sh) and [`scripts/05_gffread_fc.sh`](scripts/05_gffread_fc.sh)



and compared to the reference with gffcompare to obtain class codes (e.g., =, u, x, i) for downstream lncRNA filtering. See scripts/28_stringtie_merge.sh




##  **Pathogen assembly**

Per-sample BAMs were assembled with StringTie in reference-free mode for the pathogen ([`scripts/05_stringtie_fc.sh`](scripts/05_stringtie_fc.sh). Per-sample GTFs were then merged into a non-redundant consensus using `stringtie --merge` [`scripts/05_stringtie_merge.sh`](scripts/05_stringtie_merge.sh)


10) **Pine comparation** 
11) **Identification of long non-coding RNAs**  
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
