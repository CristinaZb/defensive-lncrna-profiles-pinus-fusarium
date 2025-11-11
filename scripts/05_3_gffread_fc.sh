#!/bin/bash

module load gffread_0.12.1

# Run gffread utility
gffread -w fc_transcripts.fa -g Fusarium_circinatum_FC072V.fa fc_transcripts.gtf
