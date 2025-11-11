#!/bin/bash

module load gffcompare_0.12.1_gcc8.2.0

# Run gffcompare utility
gffcompare -o fusarium -r fc_known_transcripts.gtf fc_transcripts.gtf -p FCIR --debug
