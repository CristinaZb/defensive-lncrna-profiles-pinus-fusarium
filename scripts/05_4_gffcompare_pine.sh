#!/bin/bash

module load haswell/gffcompare_0.12.1_gcc8.2.0

gffcompare -o ${spp} -r $REFERENCE/Pita.2_01.gtf ${spp}_transcripts.gtf -p PINUS --debug
