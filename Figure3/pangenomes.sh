#Chris might need to check this


#!/bin/bash
#SBATCH --qos=omergokc
#SBATCH --partition=omergokc
#SBATCH --cluster=faculty
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=500G
#SBATCH --job-name="pangenome_amy"
#SBATCH --output=pangenome_amy.out
#SBATCH --error=pangenome_amy.err
#SBATCH --export=NONE
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL

eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate herro 

pggb -i amylocus_haplotypes.fa \
     -o amylocus_5k \
     -n 46 \
     -t 32 \
     -p 99 \
     -s 5k \
