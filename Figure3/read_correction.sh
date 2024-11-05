#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --account=omergokc
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mem=500G
#SBATCH --gpus-per-node=1
#SBATCH --job-name="cgpu2106"
#SBATCH --output=Dorado_correct_HG02106_cluster.out
#SBATCH --error=Dorado_correct_HG02106_cluster.err
#SBATCH --export=NONE
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL

# ABOUT ################################################
# author: Dan MacGuigan
# designed to work on UB CCR cluster

# script to run Nanopore read correction using HERRO (via Dorado)
# requires basecalled reads fastq file
# outputs a fasta file of corrected reads

# example usage:
# first modify variables in the INPUTS section
# then run "sbatch AISO_Dorado_readCorrection.sh"

#### notes from Dorado Github: ############
# Dorado supports single-read error correction with the 
# integration of the HERRO algorithm. HERRO uses all-vs-all
# alignment followed by haplotype-aware correction using a 
# deep learning model to achieve higher single-read accuracies. 
# The corrected reads are primarily useful for generating de 
# novo assemblies of diploid organisms.

# Dorado correct only supports FASTX(.gz) as the input and 
# generates a FASTA file as output. An index file is generated 
# for the input FASTX file in the same folder unless one is already 
# present. Please ensure that the folder with the input file is 
# writeable by the dorado process and has sufficient disk space 
# (no more than 10GB should be necessary for a whole genome dataset).

# The error correction tool is both compute and memory intensive. 
# As a result, it is best run on a system with multiple high 
# performance CPU cores ( > 64 cores), large system memory ( > 256GB) 
# and a modern GPU with a large VRAM ( > 32GB).

# INPUTS ################################################
WD="/projects/academic/omergokc/Luane/HG02106/herro/ENA" # working directory
READS="HG02106_allreads_ENA_added" # uncorrected reads, must be located in WD
MIN_READ_LEN=10000 # min read length (bp), HERRO documentation suggests >10 kb
SP="HG02106" # species or sample prefix
CORR_READS="HG02106_ENA_10kb_HERRO.fasta" # name for corrected reads, must be .fasta
THREADS=32

# SCRIPT ################################################

module load gcc/11.2.0 fastp/0.23.2

#gunzip -c /projects/academic/omergokc/Luane/HG02106/HG02106_allreads_ENA_added.fastq.gz > /projects/academic/omergokc/Luane/HG02106/herro/ENA/HG02106_allreads_ENA_added.fastq

# filter by min read length
#fastp -i ${READS}.fastq -o ${SP}_L${MIN_READ_LEN}.fastq -h ${SP}_L${MIN_READ_LEN}.html -j ${SP}_L${MIN_READ_LEN}.json -l ${MIN_READ_LEN} --disable_quality_filtering --disable_adapter_trimming --disable_trim_poly_g -V -w 8

# convert to bgzip
#module load gcc/11.2.0 htslib/1.18
#gunzip -c ${SP}_L${MIN_READ_LEN}.fastq.gz > ${SP}_L${MIN_READ_LEN}.fastq
#bgzip -c ${SP}_L${MIN_READ_LEN}.fastq > ${SP}_L${MIN_READ_LEN}_2.fastq.gz


# run read correction
module load cuda/11.8.0
/projects/academic/tkrabben/software/Dorado/dorado-0.7.1-linux-x64/bin/dorado correct -v  HG02106_L10000.fastq > $CORR_READS

# OUTPUT ################################################
#echo "read correction complete"

