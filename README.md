# Long non-coding RNAs (lncRNAs) in the pathosystem *Pinus* spp.–*Fusarium circinatum*

Long non-coding RNAs (lncRNAs) are emerging regulators of plant immunity, modulating gene expression through chromatin interactions, cis proximity to defense genes, RNA–protein decoys, and crosstalk with small RNAs. In forest trees, and particularly in conifers, their roles during pathogen challenge remain largely unresolved. This repository documents a reproducible dual RNA-seq workflow to profile defense-related lncRNAs in resistant *Pinus pinea* and susceptible *P. radiata* during infection by *Fusarium circinatum*, and to survey fungal lncRNAs expressed in planta during the pathogenesis. The code and configuration files provided here aim to enable transparent re-analysis and serve as a foundation for testable hypotheses on lncRNA-mediated regulation of pine resistance and fungal virulence.

## Data overview

This repository is **code-only**. Raw sequencing files and large intermediates are **not** distributed here due to size and licensing constraints. The goal of this document is to record provenance and make the analysis reproducible. 
Here `${spp}` is a placeholder for the organism code (pipi for _P. pinea_, pira for _P. radiata_).

## Raw reads (not included)

- Type: paired-end RNA-Seq (dual RNA-seq: *Pinus* host + *Fusarium circinatum* pathogen).
- Number of original files: 32 libraries.
- Source: Raw reads have been deposited in the NCBI SRA Database under accession numbers SRR13737940-53 (BioProject PRJNA702546).

## References (not included)

- **Host genome/annotation**
  - Organism: *Pinus taeda* (loblolly pine)
  - Assembly name/version: **Pita v2.01**
  - NCBI Assembly record: **GCA_000404065.3**
  - Assembly level: Scaffold; submitter: **TreeGenes Database**; TreeGenes FTP: `treegenesdb.org/FTP/Genomes/Pita/v2.01/`
  - Genome FASTA: **Pita.2_01.fa.gz**
  - Gene annotation (GTF): **Pita.2_01.gtf.gz** 
 
- **Pathogen genome/annotation**
  - Name: Fusarium circinatum strain FC072V
  - Organism: *Fusarium circinatum* (TaxID: 48490)
  - BioProject: **PRJNA716280** — Genome sequencing and assembly; scope: monoisolate; registered 2021-04-27
  - Assembly accession: **GCA_018163655.1**
  - WGS master: **JAGGEA000000000**
  - BioSample: **SAMN18415966**

# Bioinformatics pipeline

## Table of Contents
- [Quality control](#quality-control)
- [Trimming](#trimming)
- [Alignment (host / pathogen)](#alignment-host--pathogen)
- [SAM processing](#sam-processing)
- [Pine assembly](#pine-assembly)
- [Pathogen assembly](#pathogen-assembly)
- [Generate a protein-coding comparator GTF](#generate-a-protein-coding-comparator-gtf)
- [Identification of long non-coding RNAs](#identification-of-long-non-coding-rnas)
  - [Step 1: Identify lncRNAs with FEELnc](#step-1-identify-lncrnas-with-feelnc)
  - [Step 2: lncRNAs identification by coding potential assessment](#step-2-lncrnas-identification-by-coding-potential-assessment)
- [Generate a FASTA with the final lncRNAs](#generate-a-FASTA-with-the-final-lncrnas)
- [Expression](#expression)

## **Quality control**

Assess read quality with FastQC on all raw FASTQ files (see [`scripts/01_fastqc.sh`](scripts/01_fastqc.sh)). 

## **Trimming**  

Reads are trimmed with Trimmomatic 0.38 using Illumina adapter removal and light head cropping (`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:10`). Post-trim QC is re-run to confirm improvement ([`scripts/01_fastqc.sh`](scripts/01_fastqc.sh)). All auxiliary `*_2u.fastq.gz` (unpaired 2) files found empty might be removed. See [`scripts/02_trimming.sh`](scripts/02_trimming.sh)

## **Alignment (host / pathogen)**

Trimmed reads are aligned with HISAT2 (`--dta`, `stranded RF`) against the host genome (_Pinus taeda_ Ptaeda2.0) and, separately, against _Fusarium circinatum_ FC072V (see [`scripts/03_hisat2_pine.sh`](scripts/03_hisat2_pine.sh) and [`scripts/03_hisat2_fc.sh`](scripts/03_hisat2_fc.sh)). 

## **SAM processing**

Alignment outputs (SAM) are converted to BAM, coordinate-sorted, and indexed with samtools. Intermediate SAM files might be removed to save space. Host and pathogen BAMs are kept separately (one per sample).

   - SAM to BAM → [`scripts/04_1_samtools_stgtie.sh`](scripts/04_1_samtools_stgtie.sh) 
   - Sorting BAM by name → [`scripts/04_2_samtools_namesort.sh`](scripts/04_2_samtools_namesort.sh)
   - Fixing BAM by removing duplicates → [`scripts/04_3_samtools_fixmate.sh`](scripts/04_3_samtools_fixmate.sh)
   - Sorting by coordinates → [`scripts/04_4_samtools_sort.sh`](scripts/04_4_samtools_sort.sh)
   - Removing duplicates → [`scripts/04_5_samtools_rmdup.sh`](scripts/04_5_samtools_rmdup.sh)

##  **Pine assembly**

Per-sample BAMs are assembled with StringTie in reference-guided mode for the pine species (`Pita.2_01.gtf.gz`; [`scripts/05_1_stringtie_pita.sh`](scripts/05_1_stringtie_pita.sh)). Per-sample GTFs are then merged into a non-redundant consensus using `stringtie --merge` [`scripts/05_2_stringtie_merge.sh`](scripts/05_2_stringtie_merge.sh).

Number of assembled transcripts in each GTF:

```bash
awk '$3=="transcript"' ${file}_transcripts.gtf | wc -l
```

In order to get a FASTA file with the sequences of the assembled transcripts of each pine, extract transcripts from the non-redundant GTF and convert to FASTA file using `gffread` and the reference genome: [`scripts/05_3_gffread_pine.sh`](scripts/05_3_gffread_pine.sh)

Check the number of transcripts in the FASTA file:

```bash
grep '^>' transcripts_${spp}_consensus.fa | wc -l
```

Compare the assembled transcripts to the reference GTF using ``gffcompare`` to obtain class codes (e.g., =, u, x, i) for downstream lncRNA filtering. See [`scripts/05_4_gffcompare_pine.sh`](scripts/05_4_gffcompare_pine.sh)

Output:

```text
${spp}.missed_introns.gff
${spp}.R_missed.gff 
${spp}.tracking 
${spp}.${spp}_transcripts.gtf.tmap 
${spp}.${spp}_transcripts.gtf.refmap 
${spp}.loci 
${spp}.annotated.gtf 
${spp}.stats
```
To better understand the outputs, visit ``https://github.com/gpertea/gffcompare``

**Inspecting gffcompare outputs and defining the “known” vs “unknown” sets**

These quick checks help verify scale and composition before extracting novel models.

1) How many loci are represented?
Each line in ``${spp}.loci`` corresponds to a locus; counting unique locus IDs gives the total genomic loci covered by the assembly:

```bash
awk '{print $1}' ${spp}.loci | sort | uniq | wc -l
```

2) What do the tracking codes look like?
``${spp}.tracking`` links each assembled transcript to its closest reference counterpart. The third printed field carries the class code or the mapping summary. This distribution is a sanity check (e.g., many “u” = intergenic, some “i/x” = intronic/antisense, etc.).

```bash
awk '{print $4}' ${spp}.tracking | sort | uniq -c
```

3) How many transcripts were tracked overall?
A coarse count of records in ``*.tracking`` (one per assembled transcript).

```bash
wc -l ${spp}.tracking
```

4) How many transcripts are exact matches to the reference? ``(class_code "=")``
Exact matches are our conservative proxy for known protein-coding models when building the comparator set.

```bash
awk '$4=="="' ${spp}.tracking | wc -l
```

We now create a list of transcript IDs that gffcompare labeled as ``class_code "="`` in the annotated GTF, and remove them to retain only novel/unknown models (the input for lncRNA calling).

Count known (“=”) in the annotated GTF:

```bash
grep 'class_code \"=\"' ${spp}.annotated.gtf | wc -l
```

Extraction (assumes ``transcript_id`` is the first attribute in column 9):

```bash
grep 'class_code "=";' ${spp}.annotated.gtf \
  | cut -f9 \
  | cut -d';' -f1 \
  | uniq \
  | sed 's/$/;/' > ${spp}_known_list.txt
```

Remove the known (“=”) transcripts and verify the remaining count:

```bash
grep -vFf ${spp}_known_list.txt ${spp}.annotated.gtf > unknown_${spp}.annotated.gtf
awk '$3=="transcript"' unknown_${spp}.annotated.gtf | wc -l
```

The file for lncRNAs identification is: _unknown_${spp}.annotated.gtf_

## **Identification of long non-coding RNAs**  

### Step 1: Identify lncRNAs with FEELnc  

Identify lncRNAs using the filter module of the `FEELnc` tool by removing the short (< 200 bp) and single-exon transcripts (see [`scripts/06_feelnc_filter.sh`](scripts/06_feelnc_filter.sh)). After that, the sequences of the resulting transcripts (potential lncRNAs) must be extracted `Gffread` ([`scripts/06_gffread.sh`](scripts/06_gffread.sh)).

How many transcripts remain in the identification pipeline?

```bash
awk '$3=="transcript"' candidate_${spp}-lncRNA.gtf | wc -l
```
The file for coding potential assessment is: ``candidate_${spp}-lncRNA.fa``

### Step 2: lncRNAs identification by coding potential assessment
   
The filtered transcripts are screened for their respective coding potential using six different computational approaches:

  (1) **Coding Potential Calculator (CPC2)** → [`scripts/07_cpc2.sh`](scripts/07_cpc2.sh)

We extract those labelled as non-coding:

```bash
awk '$8 == "noncoding" {print $1}' CPC2_${spp}-lncRNA.txt > CPC2_${spp}-noncoding.txt
```
  
  (2) **Coding-Non-Coding Index (CNCI)** → [`scripts/07_cnci.sh`](scripts/07_cnci.sh)

The resulting file `CNCI.index` has a index column (2nd) where it is specified if noncoding or coding. We extract those labeled as noncoding:

```bash
awk '$2 == "noncoding" {print $1}' CNCI.index | wc -l
awk '$2 == "noncoding" {print $1}' CNCI.index > cnci_${spp}-noncoding.txt
```
  
  (3) **Coding-Potential Assessment Tool (CPAT)** → [`scripts/07_cpat.sh`](scripts/07_cpat.sh)

CPAT analysis is a method to judge transcript encoding ability by constructing logistic regression model, calculating coding probability based on ORF length and ORF coverage. When coding probability < 0.38, it is considered as noncoding RNA.

```bash
awk '$6 < 0.38 {print $1}' ${spp}_cpat | wc -l
awk '$6 > 0.38 {print $1}' ${spp}_cpat | wc -l
```

We extract the names of transcripts assigned as noncoding:

```bash
awk '$6 < 0.38 {print $1}' ${spp}_cpat | tail -n +2 > ${spp}_cpat_noncoding.txt
```

  (4) **PLEK** → [`scripts/07_plek.sh`](scripts/07_plek.sh)

Predicted positive samples are labeled as "Coding", and negative as "Non-coding". 
Output: ${spp}_plek_predicted --> first column label, third column name of transcripts

```bash
awk '$1 == "Non-coding" {print $3}' ${spp}_plek_predicted > ${spp}_plek_noncoding.txt
awk '$1 == "Non-coding" {print $3}' ${spp}_plek_predicted | wc -l
```
  
  (5) **FEELnc codpot module** → [`scripts/07_feelnc_codpot.sh`](scripts/07_feelnc_codpot.sh)

Outputs are written to the `--outdir`:

```text
- `${spp}-lncRNA.gtf`  (lncRNA candidates)
- `${spp}-mRNA.gtf`    (coding set for comparison)
- `${spp}_RF.txt`      (coding-potential scores/labels)
```

We extract those labelled as non-coding (label--> 1 = coding; 0 = noncoding):

```bash
awk '$11 == "0" {print $1}' ${spp}_RF.txt > feelnc_${spp}_noncoding.txt
```
  
  (6) **EnTAP** → [`scripts/07_entap_codpot.sh`](scripts/07_entap_codpot.sh)

The unannotated sequences (potential non-coding sequences) are located in these files:

```text
frame_selection/GeneMarkS-T/processed/sequences_removed.fn
final_results/final_unannotated.faa
```

To obtain the FASTA sequences of the unannotated sequences:
  
```bash
grep "^>" entap_outfiles/frame_selection/GeneMarkS-T/processed/sequences_removed.fnn | sed 's/^>//' > sequences_removed.txt
grep "^>" entap_outfiles/final_results/final_unannotated.faa | sed 's/^>//' > final_unannotated.txt
cat entap_outfiles/frame_selection/GeneMarkS-T/processed/sequences_removed.txt entap_outfiles/final_results/final_unannotated.txt > entap_outfiles/final_results/entap_${spp}_noncoding.txt
```

## **Build a consensus lncRNA set across predictors**

Each predictor (CNCI, CPAT, CPC2, PLEK, FEELnc codpot, EnTAP unannotated) has distinct biases. To be conservative, we define consensus lncRNAs as the intersection of the six non-coding sets per species. The resulting ID lists (``${spp}_consensus_noncoding.txt``) feed directly into the FASTA extraction step below.

Inputs (per species):

```text
cnci_${spp}-noncoding.txt

${spp}_cpat_noncoding.txt

CPC2_${spp}-noncoding.txt

${spp}_plek_noncoding.txt

feelnc_${spp}_noncoding.txt

entap_${spp}_noncoding.txt
```

```r
# ---- Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
})

# Helper: read 1-column files robustly and normalize IDs
read_ids <- function(path, strip_gt = FALSE) {
  x <- read.table(path, header = FALSE, sep = "", quote = "", comment.char = "", stringsAsFactors = FALSE)
  ids <- x[[1]]
  if (strip_gt) ids <- gsub("^>", "", ids)
  ids <- trimws(ids)
  ids[nzchar(ids)]
}

consensus_from_files <- function(tag) {
  # tag = "pira" or "pipi"
  files <- list(
    CNCI   = sprintf("cnci_%s-noncoding.txt", tag),
    CPAT   = sprintf("%s_cpat_noncoding.txt", ifelse(tag=="pira","Pira","Pipi")),
    CPC2   = sprintf("CPC2_%s-noncoding.txt", tag),
    PLEK   = sprintf("%s_plek_noncoding.txt", tag),
    FEELnc = sprintf("feelnc_%s_noncoding.txt", tag),
    EnTAP  = sprintf("entap_%s_noncoding.txt", tag)
  )

  sets <- list(
    CNCI   = read_ids(files$CNCI),
    CPAT   = read_ids(files$CPAT),
    CPC2   = read_ids(files$CPC2),
    PLEK   = read_ids(files$PLEK, strip_gt = TRUE),  # PLEK headers sometimes start with '>'
    FEELnc = read_ids(files$FEELnc),
    EnTAP  = read_ids(files$EnTAP)
  ) %>% lapply(unique)

  # Intersection (consensus)
  consensus <- Reduce(intersect, sets)
  message(sprintf("[%s] sizes: CNCI=%d, CPAT=%d, CPC2=%d, PLEK=%d, FEELnc=%d, EnTAP=%d; CONSENSUS=%d",
                  tag, lengths(sets)["CNCI"], lengths(sets)["CPAT"], lengths(sets)["CPC2"],
                  lengths(sets)["PLEK"], lengths(sets)["FEELnc"], lengths(sets)["EnTAP"], length(consensus)))

  out <- sprintf("%s_consensus_noncoding.txt", tag)
  write.table(consensus, out, quote = FALSE, col.names = FALSE, row.names = FALSE)
  invisible(list(consensus = consensus, sets = sets, out = out))
}

# ---- Run for both species ----
res_pira <- consensus_from_files("pira")  # Pinus radiata
res_pipi <- consensus_from_files("pipi")  # Pinus pinea
```

We next extract FASTA sequences only for transcripts present in the consensus lists created above (``${spp}_consensus_noncoding.txt``), using the species-specific transcript FASTA files (``transcripts_${spp}_consensus.fa``). Check that there are no printable characters in the text file and remove them:

```bash
cat -A ${spp}_consensus_noncoding.txt
dos2unix ${spp}_consensus_noncoding.txt
```

Convert the FASTA format of the transcripts into single-line FASTA:

```bash
awk '/^>/ {printf("%s%s\n",(NR>1?"\n":""),$0);next;} {printf("%s",$0);} END {printf("\n");}' transcripts_${spp}_consensus.fa > transcripts_${spp}_consensus.fa
```

Extract sequences and check:

```bash
grep -F -f ${spp}_consensus_noncoding.txt -A 1 transcripts_${spp}_consensus.fa | grep -v "^--$" > ${spp}_consensus_noncoding.fa
grep -c "^>" ${spp}_consensus_noncoding.fa
```

## **Expression**

To calculate the expression of transcripts, we need the BAM files used to generate the GTFs, the non-redundant GTF, and `StringTie -e`. We prepare the non-redundant GTF because it has a header that needs to be removed:

```bash
tail -n +3 ${spp}_transcripts.gtf > ${spp}_transcripts-mod.gtf
```

Now we calculate the expression for all identified transcripts. To determine the expression of lncRNAs, transcripts will be filtered later in the RStudio script. See [`scripts/08_stringtie_expression.sh`](scripts/08_stringtie_expression.sh)
To convert the GTFs with expression information into a quantification table, prepare the files that we need: `sample_lst_gtf.txt` and `prepDE.py` (visit https://github.com/gpertea/stringtie). After that, we run a script to run the python script:

```bash
python2 prepDE.py -i sample_lst_gtf.txt -s ${spp}
```

Now we export the final table `transcript_count_matrix.csv` to RStudio in order to carry out the differential expression analysis (DEA).

##  **Pathogen assembly**

Per-sample BAMs were assembled with StringTie in reference-free mode for the pathogen ([`scripts/05_1_stringtie_fc.sh`](scripts/05_1_stringtie_fc.sh). Per-sample GTFs were then merged into a non-redundant consensus using `stringtie --merge` [`scripts/05_2_stringtie_merge.sh`](scripts/05_2_stringtie_merge.sh) In order to get a FASTA file with the sequences of the assembled transcripts of the pathogen, extract transcripts from the non-redundant GTF and convert to FASTA file using `gffread` and the reference genome: [`scripts/05_3_gffread_fc.sh`](scripts/05_3_gffread_fc.sh)

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

Annotation of pathogen transcripts to identify those with functional hits/assignments → [`scripts/05_4_entap_fc.sh`](scripts/05_4_entap_fc.sh)

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
Keep only these transcripts and check:

```bash
grep -Ff list_annotated_2.txt fc_transcripts.gtf > fc_known_transcripts.gtf
awk '$3=="transcript"' fc_known_transcripts.gtf | wc -l
```
Now we use this new GTF to compare with the previous one.

**2) Structural comparison against a “known transcripts” GTF**

Contrast the assembled GTF against a reference known set using ``gffcompare`` to classify each transcript structurally (match/novel, class codes) and generate ``.tracking`` and ``.loci files``  →  [`scripts/05_5_gffcompare_fc.sh`](scripts/05_5_gffcompare_fc.sh)

Number of loci:

```bash
awk '{print $1}' fusarium.loci | wc -l 
```

Number of unique loci:

```bash
awk '{print $1}' fusarium.loci | sort -u | wc -l
```

The ``.tracking`` third column indicates closest reference match

```bash
awk '{print $4}' fusarium.tracking | sort | uniq -c
```

**3) Create a new GTF without the known/coding transcripts ("=")**

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

The process of pathogen lncRNA identification and expression assessment is similar to the one used for pine species, this time without using EnTAP as it was already used in previous step.

## How to cite
If you use this code, please cite this repository and the related manuscript (when available).  
*(A `CITATION.cff` file will be added for automated citation metadata.)*

## License
MIT — see `LICENSE`.

## Contact
Maintainer: Cristina Zamora Ballesteros (cristinazamoraballesteros@gmail.com) — issues and PRs welcome.
