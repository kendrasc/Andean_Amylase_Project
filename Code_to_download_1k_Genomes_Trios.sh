#!/bin/bash

batch_number=$1
output_dir=$2
file_directory=$3

while read samples
do
    grep "$samples" Trios_1k_genomes.txt >> Trios_batch_${batch_number}_all_info.txt
done < ${output_dir}/Trios_batch_${batch_number}.txt

grep "fastq.gz" Trios_batch_${batch_number}_all_info.txt | awk -F '[/\t;]' '{print $2 "/" $3 "/" $4 "/" $5 "/" $6 "/" $7}' > sample_links.txt
test_count1=$(wc -l sample_links.txt | awk '{print $1}')
grep "fastq.gz" Trios_batch_${batch_number}_all_info.txt  | awk -F '[/\t;]' '{print $9 "/" $10 "/" $11 "/" $12 "/" $13 "/" $14}' >> sample_links.txt
test_count2=$(wc -l sample_links.txt | awk '{print $1}')

if (($test_count2/$test_count1 != 2)); then
    echo "Incorrect number of sample links"
    exit 1
fi

Number_of_files=$(wc -l sample_links.txt | awk '{print $1}')
yes ${file_directory}"/" | head -n ${Number_of_files} > tail_location.txt

grep "fastq" Trios_batch_${batch_number}_all_info.txt | awk '{print $2"/"}' > Individuals.txt
cat Individuals.txt Individuals.txt > Individuals_twice.txt

awk -F '[;]' '{print $1}' Trios_batch_${batch_number}_all_info.txt | awk -F '[/]' '{print $NF}' > fastq_file.txt
awk '{print $1}' Trios_batch_${batch_number}_all_info.txt | awk -F '[/]' '{print $NF}' >> fastq_file.txt
fastq_file_length=$(wc -l fastq_file.txt | awk '{print $1}')

if (($test_count2 != $fastq_file_length)); then
    echo "Incorrect number of fastqs"
    exit 1
fi

# Combine the files
paste -d'\0' tail_location.txt Individuals_twice.txt fastq_file.txt > adjusted_home_path.txt
paste sample_links.txt adjusted_home_path.txt > ${batch_number}_to_download.txt

# To make directories
while read NAME ; do
    mkdir -p /vscratch/${file_directory}/${NAME};
done < ${output_dir}/Trios_batch_${batch_number}.txt

# Move download file to Trios directory
mv "${batch_number}_to_download.txt" "${output_dir}"

globus transfer $ENA: $Buffalo: --label ${batch_number} --batch ${output_dir}/${batch_number}_to_download.txt

rm tail_location.txt
rm Individuals.txt
rm Individuals_twice.txt
rm adjusted_home_path.txt
rm sample_links.txt
rm Trios_batch_${batch_number}_all_info.txt
rm fastq_file.txt
