#!/bin/bash

module load gffread_0.12.1

# Run gffread utility

gffread -w candidate_${spp}-lncRNA.fa -g Pita.1_2.fa candidate_${spp}-lncRNA.gtf
