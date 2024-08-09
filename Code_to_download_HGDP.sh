#!/bin/bash

pop=$1
directory=$2
mkdir /projects/${directory}/${pop}

awk 'NR==FNR { variables[$1]=$0; next } { for (var in variables) if (index($0, var) > 0) print $0 }' to_fix.txt ${pop}.tsv > ${pop}_to_fix.txt
awk '{print $1}' ${pop}_to_fix.txt | grep ".fastq.gz" | awk -F '[/]' '{print $4"/"$5"/"$6"/"$7"/"$8"/"$9}' | sed 's/\/$//' > sample_links.txt

Number_of_files=$(wc -l sample_links.txt | awk '{print $1}')
yes /${directory}/${pop}"/" | head -n ${Number_of_files} > tail_location.txt

grep ".fastq.gz" ${pop}_to_fix.txt | awk '{print $11"/"}' > Individuals.txt
awk '{print $1}' ${pop}_to_fix.txt | grep ".fastq.gz" | awk -F '[/]' '{print $NF}' > fastq_file.txt
paste -d'\0' tail_location.txt Individuals.txt fastq_file.txt > adjusted_home_path.txt
paste sample_links.txt adjusted_home_path.txt > ${pop}_to_pull.txt

awk '{print $11}' ${pop}_to_fix.txt | tail -n+2 | sort | uniq > ${pop}_samples.txt
while read NAME ; do mkdir /projects/${directory}/${pop}/${NAME}; done < ${pop}_samples.txt

globus transfer $ENA: $Buffalo: --label ${pop} --batch ${pop}_to_pull.txt
