#!/bin/bash

# Run EnTAP

EnTAP --runP \
-i fc_transcripts.fa \
-d EnTAP/1.9.2/databases/diamond/nr.dmnd \
-d EnTAP/1.9.2/databases/diamond/refseq_v2.dmnd \
-d EnTAP/1.9.2/databases/diamond/uniprot_sprot.dmnd \
--level 0 \
--ontology 0 \
--ontology 1 \
--threads $SLURM_NTASKS \
-c bacteria \
--taxon Fusarium
