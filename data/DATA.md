# Data overview

This repository is **code-only**. Raw sequencing files and large intermediates are **not** distributed here due to size and licensing constraints. The goal of this document is to record provenance and make the analysis reproducible.

---

## Raw reads (not included)

- Type: paired-end RNA-Seq (dual RNA-seq: *Pinus* host + *Fusarium circinatum* pathogen).
- Number of original files: `<e.g., 32 paired libraries>`; total size: `<e.g., ~135 GB>`.
- Example filenames: `SampleX_1.fastq.gz`, `SampleX_2.fastq.gz`.
- Source: Raw reads have been deposited in the NCBI SRA Database under accession numbers SRR13737940-53 (BioProject PRJNA702546).

### Preprocessing decisions (before analysis)
- Removed auxiliary/empty files: `<e.g., *_2u.fastq.gz found empty → deleted>`.
- Standard trimming and QC performed (see scripts and `docs/methods.md` for tools/parameters).

## References (not included)

- **Host genome/annotation**
  - Organism: *Pinus taeda* (loblolly pine) — TreeGenes code: **Pita**
  - Assembly name/version: **Ptaeda2.0**
  - Primary accession: **GCA_000404065.3** (NCBI Assembly)
  - Assembly level: Scaffold; submitter: **TreeGenes Database**; release date: **2017-01-09**
  - NCBI Assembly record: GCA_000404065.3 (Ptaeda2.0). 
  - TreeGenes FTP: `treegenesdb.org/FTP/Genomes/Pita/v2.0/` (subdir `genome/`).

- **Pathogen genome/annotation**
  - Name: Fusarium circinatum strain FC072V
  - Organism: *Fusarium circinatum* (TaxID: 48490)
  - BioProject: **PRJNA716280** — Genome sequencing and assembly; scope: monoisolate; registered 2021-04-27
  - Assembly accession: **GCA_018163655.1** (assembly level: Contig)
  - WGS master: **JAGGEA000000000**
  - BioSample: **SAMN18415966**
  - Source: NCBI BioProject / Assembly records
