#!/bin/bash

module load gcc
module load bcftools
module load gatk
module load samtools

test="biallelic_no_missingness_no_EUR_phased"
chrm="chr1"
gene="AMY"
location1="chr1:103347982-103830994"
info_file="IDs.txt"
file="AllQuechua_Maya_PhasedChr1_no_missing_biallelic_103347982_103830994_no_EUR.vcf.gz"

file_name="${gene}_region_${chrm}"
for sample in `bcftools query -l $file`; do
   bcftools view -Oz -s $sample $file -o ${gene}_region/${file_name}_${sample}.vcf.gz

 #Use awk to find the matching sample_id and extract Copy_number and Population
    Copy_number=$(awk -v id="$sample" '$1 == id {print $2}' $info_file)
    Population=$(awk -v id="$sample" '$1 == id {print $3}' $info_file)
    
        # remove the sites that are the same as the reference
    gatk IndexFeatureFile \
    -I ${gene}_region/${file_name}_${sample}.vcf.gz

   gatk SelectVariants \
    -R chr1_hg38.fasta \
    -V ${gene}_region/${file_name}_${sample}.vcf.gz \
    --select-type-to-include SNP \
    --exclude-non-variants true \
    --O ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz

        # For haplotype 1
    samtools faidx chr1_hg38.fasta ${location1} | bcftools consensus -H 1 ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz > ${gene}_region/${sample}_region_hap1.fasta

    # For haplotype 2
    samtools faidx chr1_hg38.fasta ${location1} | bcftools consensus -H 2 ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz > ${gene}_region/${sample}_region_hap2.fasta

    # Ensure the variables are correctly formatted
    formatted_label=$(printf "%s_%s_%s_hap1" "$sample_id" "$Copy_number" "$Population")
    # Relabel the header in each FASTA file for haplotype 1
    sed -e 's/\r//g' -e "s/^>.*/>${formatted_label}/" ${gene}_region/${sample_id}_region_hap1_${test}.fasta > ${gene}_region/${sample_id}_region_relabelled_hap1_${test}.fasta

    # For haplotype 2
    formatted_label=$(printf "%s_%s_%s_hap2" "$sample_id" "$Copy_number" "$Population")
    sed -e 's/\r//g' -e "s/^>.*/>${formatted_label}/" ${gene}_region/${sample_id}_region_hap2_${test}.fasta > ${gene}_region/${sample_id}_region_relabelled_hap2_${test}.fasta

    rm ${gene}_region/AMY_region_chr1_${sample_id}_${test}.vcf.gz.tbi
    rm ${gene}_region/AMY_region_chr1_${sample_id}_${test}.vcf.gz
    rm ${gene}_region_chr1_${sample_id}_excluded_nonvariants_${test}.vcf.gz
    rm ${gene}_region_chr1_${sample_id}_excluded_nonvariants_${test}.vcf.gz.tbi
    rm ${gene}_region/${sample_id}_region_hap1_${test}.fasta
    rm ${gene}_region/${sample_id}_region_hap1_${test}.fasta.fai
    rm ${gene}_region/${sample_id}_region_hap2_${test}.fasta
    rm ${gene}_region/${sample_id}_region_hap2_${test}.fasta.fai
done

cat AMY_region/*${test}.fasta > amy_trees/Quechua_Maya_${test}.fasta

conda activate iqtree
iqtree -s amy_trees/Quechua_Maya_${test}.fasta --seqtype DNA -m TEST --alrt 1000 -bb 1000 -bnni -allnni -nt AUTO -pre amy_trees/Quechua_Maya_haploid_amy_region_within_recombination_no_missing_${test}

