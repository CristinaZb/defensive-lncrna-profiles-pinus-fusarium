# Pinus–Fusarium lncRNA (dual RNA-seq)

Reproducible dual RNA-seq analysis of defensive lncRNAs in *Pinus* challenged by *Fusarium circinatum*.  
**Scope:** code-only repository (no raw data).

## Overview
This repo contains the scripts and minimal configuration to reproduce:
QC → read processing → counting → differential expression (DESeq2) → co-expression (WGCNA) → lncRNA cis/trans context and module summaries for resistant vs susceptible *Pinus* spp.

## What’s here
- `scripts/` — Bash/R scripts used in the paper (cleaned and documented).
- `config/` — small text configs (e.g., parameters, contrasts).
- `results/` — lightweight exports only (tables/figures under 10–20 MB).
- `docs/` — brief methods and notes.

> **Note:** Raw FASTQs, genome indices, and large intermediates are intentionally **not** tracked.

## Quick start
1. Clone or download this repo.
2. Check `scripts/` headers for dependencies (R ≥ 4.3; packages listed at the top of each script).
3. Edit `config/` files (paths/contrasts) if you want to rerun parts of the analysis.
4. Run scripts in the numbered order (see filenames).

## Reproducibility
- Versions are tagged; releases will be archived in Zenodo to obtain a DOI.
- Where applicable, random seeds are set within scripts.
- If you find any non-determinism, open an issue.

## How to cite
If you use this code, please cite this repository and the related manuscript (when available).  
*(A `CITATION.cff` file will be added for automated citation metadata.)*

## License
MIT — see `LICENSE`.

## Contact
Maintainer: Cristina Zamora Ballesteros (cristinazamoraballesteros@gmail.com) — issues and PRs welcome.
