# Pinus spp.–Fusarium long non-coding RNAs

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

Per-sample BAMs were assembled with StringTie in reference-guided mode for the pine species (`Pita.2_01.gtf.gz`; [`scripts/05_1_stringtie_pita.sh`](scripts/05_1_stringtie_pita.sh)). Per-sample GTFs were then merged into a non-redundant consensus using `stringtie --merge` [`scripts/05_2_stringtie_merge.sh`](scripts/05_2_stringtie_merge.sh).

Number of assembled transcripts in each GTF:

```bash
awk '$3=="transcript"' ${file}_transcripts.gtf | wc -l
```

In order to get a fasta file with the sequences of the assembled transcripts of each organism, extract transcripts from the non-redundant GTF and convert to FASTA file using `gffread` and the reference genome: [`scripts/05_3_gffread_pine.sh`](scripts/05_3_gffread_pine.sh)

Check the number of transcripts in the FASTA files:

```bash
grep '^>' transcripts_${spp}_consensus.fa | wc -l
```

Compare the _de novo_ transcripts to the reference GTF using ``gffcompare`` to obtain class codes (e.g., =, u, x, i) for downstream lncRNA filtering. See [`scripts/05_4_gffcompare_pine.sh`](scripts/05_4_gffcompare_pine.sh)

Output:

${spp}.missed_introns.gff
${spp}.R_missed.gff 
${spp}.tracking 
${spp}.${spp}_transcripts.gtf.tmap 
${spp}.${spp}_transcripts.gtf.refmap 
${spp}.loci 
${spp}.annotated.gtf 
${spp}.stats 

•	How many loci do we have?
```bash
awk '{print $1}' ${spp}.loci | sort | uniq | wc -l
```
•	Pipi.tracking --> Thirth column indicates the information of the most close reference annotation transcript.
```bash
awk '{print $4}' ${spp}.tracking | sort | uniq -c
```
•	How many transcripts do we have?
```bash
wc -l ${spp}.tracking
```
•	How many transcripts do we have with a full match "="?
```bash
awk '$4=="="' ${spp}.tracking | wc -l
```

Generate a file with the IDs of all known transcripts and then delete them from the gtf:
```bash
grep 'class_code \"=\"' ${spp}.annotated.gtf | wc -l
```
Extract the fist field (until ;) from the column 9:
```bash
grep 'class_code \"=\"' ${spp}.annotated.gtf | cut -d$'\t' -f9 | cut -d ";" -f1 | uniq | sed 's/$/;/' > ${spp}_known_list.txt
```
Remove the known transcripts by:
```bash
grep -vFf ${spp}_known_list.txt ${spp}.annotated.gtf > unknown_${spp}.annotated.gtf 
```
Check:
```bash
awk '$3=="transcript"' unknown_${spp}.annotated.gtf | wc -l
```
The file for lncRNAs identification is: _unknown_${spp}.annotated.gtf_

##  **Pathogen assembly**

Per-sample BAMs were assembled with StringTie in reference-free mode for the pathogen ([`scripts/05_1_stringtie_fc.sh`](scripts/05_1_stringtie_fc.sh). Per-sample GTFs were then merged into a non-redundant consensus using `stringtie --merge` [`scripts/05_2_stringtie_merge.sh`](scripts/05_2_stringtie_merge.sh) In order to get a fasta file with the sequences of the assembled transcripts of the pathogen, extract transcripts from the non-redundant GTF and convert to FASTA file using `gffread` and the reference genome: [`scripts/05_3_gffread_fc.sh`](scripts/05_3_gffread_fc.sh)

Check the number of transcripts in the FASTA files:

```bash
grep '^>' fc_transcripts.fa | wc -l
```

**>> Generate a protein-coding comparator GTF**

Build a curated GTF of protein-coding transcripts to (i) label/retain known coding loci, (ii) define the “known set” for downstream comparisons, and (iii) separate novel/unknown models (candidates for lncRNA).

**Inputs:**
- GTF of assembled pathogen transcripts ``fc_transcripts.gtf``
- FASTA of assembled pathogen transcripts ``fc_transcripts.fa``
- Pathogen genome FASTA ``Fusarium_circinatum_FC072V.fa``
- A reference known transcripts GTF for F. circinatum (e.g., curated/annotated set) → fc_known_transcripts.gtf

**Outputs:**
unknown_fusarium.annotated.gtf → GTF of novel/unknown transcripts (used later for lncRNA identification)

**1) Functional annotation with EnTAP**

Annotation of pathogen transcripts to identify those with functional hits/assignments → [`scripts/05_4_entap_fc.sh`](scripts/scripts/05_4_entap_fc.sh)

Get first functional view of what is likely coding/known vs unknown:

```bash
grep -c '^>' entap_outfiles/final_results/final_annotated.fnn
grep -c '^>' entap_outfiles/final_results/final_unannotated.fnn
```

Extract the transcripts that have been annotated:
```bash
grep "^>" final_annotated.fnn > list_annotated.txt
```
Change the format so that it is recognizable in the GTF: 
Remove the “>” symbol at the beginning of the row, and add ‘transcript_id’, quotation marks, and “;” at the end.
```bash
awk '{gsub(/^>/, ""); print "transcript_id \"" $0 "\";"}' list_annotated.txt > list_annotated_2.txt
```
Keep only these transcripts:
```bash
grep -Ff list_annotated_2.txt fc_transcripts.gtf > fc_known_transcripts.gtf
```
Check:
```bash
awk '$3=="transcript"' fc_known_transcripts.gtf | wc -l
```
Now we use this new GTF to compare with the previous one.


**2) Structural comparison against a “known transcripts” GTF**

Contrast the assembled GTF against a reference known set using gffcompare to classify each transcript structurally (match/novel, class codes) and generate .tracking and .loci files  →  [`scripts/05_5_gffcompare_fc.sh`](scripts/scripts/05_5_gffcompare_fc.sh)

Number of loci:

```bash
awk '{print $1}' fusarium.loci | wc -l 
```

Number of unique loci:

```bash
awk '{print $1}' fusarium.loci | sort -u | wc -l
```

The .tracking third column indicates closest reference match

```bash
awk '{print $4}' fusarium.tracking | sort | uniq -c
```

**3) Create a new GTF without the known transcripts ("=")**
```bash
grep 'class_code \"=\"' fusarium.annotated.gtf | wc -l
```
Usamos el archivo anterior: list_annotated_2.txt
```bash
grep -vFf list_annotated_2.txt fusarium.annotated.gtf > unknown_fusarium.annotated.gtf
```
Check:
```bash
awk '$3=="transcript"' unknown_fusarium.annotated.gtf | wc -l 
```
The file for lncRNAs identification is: _unknown_fusarium.annotated.gtf_

## **Identification of long non-coding RNAs**  

### Step 1:  → [`scripts/feelnc_filter.sh`](scripts/feelnc_filter.sh)  

Identify lncRNAs using the filter module of the `FEELnc` tool by removing the short (< 200 bp) and single-exon transcripts. After that, the sequences of the resulting transcripts (potential lncRNAs) were extracted `Gffread` ([`scripts/05_3_gffread_pine.sh`](scripts/05_3_gffread_pine.sh)).


   

Step 2:
  




## License
MIT — see `LICENSE`.

## Contact
Maintainer: Cristina Zamora Ballesteros (cristinazamoraballesteros@gmail.com) — issues and PRs welcome.
