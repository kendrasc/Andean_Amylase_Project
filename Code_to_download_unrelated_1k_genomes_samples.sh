#!/bin/bash

# Sample Population
pop=$1
# Directory Location
directory=$2

grep "${pop}" All_1k_genomes_samples.txt | awk '{print $1}' | sort | uniq > ${pop}_all.txt
comm -12 ${pop}_all.txt 1k_Genomes_unrelated_individuals.txt > ${pop}_unrelated.txt

while read individual
do
    grep "${individual}" 1k_Genomes_unrelated.txt | awk '{print $2}' | awk -F 'uk/|;' '{print $2}' > sample_links.txt
    test_count1=$(wc -l sample_links.txt | awk '{print $1}')
    grep "${individual}" 1k_Genomes_unrelated.txt | awk '{print $2}' | awk -F 'uk/|;' '{print $4}'  >> sample_links.txt
    test_count2=$(wc -l sample_links.txt | awk '{print $1}')

# Checks to make sure that both the first and second runs have been pulled
    if (($test_count2/$test_count1 != 2)); then
        echo "Incorrect number of sample links"
        exit 1
    fi

    Number_of_files=$(wc -l sample_links.txt | awk '{print $1}')

    if [ -n "${pop}" ]; then
        Pop="${pop}/"
    else
        Pop=""
    fi

# Pulls the correct directory for the population
    yes ${directory}${Pop}${individual}"/" | head -n ${Number_of_files} > tail_location.txt

# Gets the ID number for each individual
    grep "${individual}" 1k_Genomes_unrelated.txt | awk '{print $1}' > Individuals.txt
    cat Individuals.txt Individuals.txt > Individuals_twice.txt

    awk -F '[/]' '{print $NF}' sample_links.txt > fastq_file.txt
    fastq_file_length=$(wc -l fastq_file.txt | awk '{print $1}')

    if (($test_count2 != $fastq_file_length)); then
        echo "Incorrect number of fastqs"
        exit 1
    fi

    paste -d'\0' tail_location.txt Individuals_twice.txt fastq_file.txt > adjusted_home_path.txt

    paste sample_links.txt adjusted_home_path.txt > ${individual}_to_pull.txt

    mkdir -p ${pop}_pulling
    mv ${individual}_to_pull.txt ${pop}_pulling/
    
    mkdir -p /vscratch${directory}${Pop}${individual}

done < ${pop}_unrelated.txt

cat  ${pop}_pulling/*_to_pull.txt >  ${pop}_to_pull.txt

globus transfer $ENA: $Buffalo: --label ${pop} --batch  ${pop}_to_pull.txt