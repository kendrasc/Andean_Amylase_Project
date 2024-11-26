#!/bin/bash

READS="HG02106_allreads" # uncorrected reads
MIN_READ_LEN=10000 # min read length (bp), HERRO documentation suggests >10 kb
SP="HG02106" # species or sample prefix
CORR_READS="HG02106_10kb_HERRO.fasta" # name for corrected reads, must be .fasta
THREADS=32

fastp -i ${READS}.fastq -o ${SP}_L${MIN_READ_LEN}.fastq -h ${SP}_L${MIN_READ_LEN}.html -j ${SP}_L${MIN_READ_LEN}.json -l ${MIN_READ_LEN} --disable_quality_filtering --disable_adapter_trimming --disable_trim_poly_g -V -w 8

./dorado-0.7.1-linux-x64/bin/dorado correct -v  ${SP}_L${MIN_READ_LEN}.fastq > $CORR_READS


