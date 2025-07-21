#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
 #The default is that the job will not be requeued.
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR
TMPDIR=/tmp
echo "TMPDIR="$TMPDIR
ulimit -s unlimited

module load gcc/11.2.0
module load clang
#module load bcftools
#module load vcftools

eval "$(/projects/academic/omergokc/charikleia/anaconda3/bin/conda shell.bash hook)"
conda activate gsl

export LD_LIBRARY_PATH=/projects/academic/omergokc/Kendra/software/2023.01/anaconda/envs/gsl/lib:$LD_LIBRARY_PATH
export CFLAGS="-I/projects/academic/omergokc/Kendra/software/2023.01/anaconda/envs/gsl/include"
export LDFLAGS="-L//projects/academic/omergokc/Kendra/software/2023.01/anaconda/envs/gsl/lib"

Info=$1

./selscan --xpehh \
          --vcf filtered_genetic_map_AllQuechua_no_EUR_PhasedChr1_no_indels_0.05_no_missing_unique.vcf.gz \
          --vcf-ref filtered_genetic_map_Maya_PhasedChr1_no_indels_0.05_no_missing_unique.vcf.gz \
          --out Quechua_Maya_${Info} \
          --max-extend 2000000 \
          --wagh \
          --map Quechua_Maya_chr1_bedtools_test_genetic_map_no_header.txt  \
          --threads 8

./norm --xpehh --files Quechua_Maya_${Info}.xpehh.out
