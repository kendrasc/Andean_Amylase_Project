#!/bin/bash

individual=$1
population=$2
directory=$3
export individual
echo "VAR = \"${individual}\""

# Sets the directory the population & current individual for this run is located (will change for HGDP files)
input_directory="${directory}/${population}/${individual}"
echo "Input directory = \"${input_directory}\""
cd $input_directory

# Sets the directory that all of the scripts will be located in + checks
kmer_direct="/path/to/Kmer/GeneToCN"
echo "Kmer directory = \"${kmer_direct}\""
echo "Current directory:" && pwd

# Have to make sure that just the fastq files end in ".gz" in the directory of interest
# There is probably a better way to do this but this is what I've got...
all=$(ls | wc -l | awk '{print $1}')
echo "Total length of directory: $all"
paste <(yes $input_directory/ | head -n ${all}) <(ls) > to_be_edited.txt
awk '{print $1$2}' to_be_edited.txt | grep ".gz" > second_file.txt
sed '/fastqs\.txt/d' ${population}_2.txt | tr '\n' ' ' > locations.txt
rm to_be_edited.txt
rm second_file.txt

# added to make sure other fastq files aren't grabbed by accident
ls | grep ".gz" | sed '/Fastq\.err/d' | sed '/Fastq\.out/d' | sed '/fastqs\.txt/d'  > fastqs.txt
fastq_files=$(awk '{print $1}' fastqs.txt)
echo $fastq_files

# Load Conda environment here
if [ $? -ne 0 ]; then
    echo "Error: Failed to start conda"
    exit 1
fi

conda activate GeneToCN
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate GeneToCN"
    exit 1
fi

file=$(<"${input_directory}/locations.txt")
if [ $? -ne 0 ]; then
    echo "Error: Failed to set file variable"
    exit 1
fi
echo "Fastq locations = \"$file\""

cd ${kmer_direct}
if [ $? -ne 0 ]; then
    echo "Error: Failed to go to original directory"
    exit 1
fi
echo "Current directory" && pwd

# Create directory for individual
mkdir -p ${kmer_direct}/Amylase_Counts/${population}/${individual}

echo "Began Python script"
python KmerToCN.py -db Amylase_Counts_all_kmers.db -kp kmer_db_locations.txt -s ${file} -d Amylase_Counts/${population}/${individual} -o ${individual} -gm /Path/to/Genometester4/src/ -r Amylase_Counts_region.txt -i

if [ $? -ne 0 ]; then
    echo "Error: Failed to run python script"
    exit 1
fi
