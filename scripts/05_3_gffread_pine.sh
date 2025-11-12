#!/bin/bash

module load gffread_0.12.1

gffread -w transcripts_${spp}_consensus.fa -g Pita.2_01.fa ${spp}_transcripts.gtf
