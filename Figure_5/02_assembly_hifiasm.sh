#!/bin/bash
cor_reads="path/to/your/corrected/reads.fasta" #Reads corrected using herro
par1="your/path/to/parent1/short/reads.fastq.gz" #parent's 1 short reads for phasing
par1="your/path/to/parent2/short/reads.fastq.gz" #parent's 2 short reads for phasing
raw="raw_reads.fastq.gz" #path to your raw reads
sample="samplename"

#pipeline adapted from https://doi.org/10.1101/2024.05.18.594796
#hifiasm: https://github.com/chhylp123/hifiasm
#yak: https://github.com/lh3/yak. 

#If coverage is higher than 38x, then the reads need to be downsampled
#Corrected reads were first chopped to a length of 30,000 bp, and chunks shorter than 10,000 bp were filtered away using SeqKit v2.5.1, as suggested in Mile Šikić, 2024 (doi: https://doi.org/10.1101/2024.05.18.594796)
seqkit sliding -s 30000 -W 30000 -g ${cor_reads} > chopped.fasta
seqkit seq -m 10000 chopped.fasta > processed.fasta

#As ultra long reads, we are using only reads bigger than 50kb
seqkit seq -m 50000 $raw > raw_reads_${sample}_50kb.fastq

#please see to install yak: https://github.com/lh3/yak. 
/projects/academic/omergokc/Luane/softwares/yak/yak count -b37 -t 32 -o parent1.yak <(zcat ${par1})  
/projects/academic/omergokc/Luane/softwares/yak/yak count -b37 -t 32 -o parent2.yak <(zcat ${par2})   

./hifiasm -t 32 -o ${sample}_diploid_50kb --ul raw_reads_${sample}_50kb.fastq --ul-cut 50000 --dual-scaf -1 parent1.yak -2 parent2.yak  processed.fasta 

#For haploid samples:

#For Haploid samples
#./hifiasm -t 32 -o ${sample}_hifiasm_50kb --ul raw_reads_${sample}_50kb.fastq -l0 ${cor_reads}
