module load gcc
module load bcftools
module load gatk
module load samtools

chrm="chr1"
gene="AMY"
location1="chr1:103347982-103830994"
info_file="LWK.txt"
file="LWK_no_missing_103347982_103830994_biallelic_maf_5.vcf.gz"
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
    formatted_label=$(printf "%s_%s_%s_hap1" "$sample" "$Copy_number" "$Population")
    # Relabel the header in each FASTA file for haplotype 1
    sed -e 's/\r//g' -e "s/^>.*/>${formatted_label}/" ${gene}_region/${sample}_region_hap1.fasta > ${gene}_region/${sample}_region_relabelled_hap1.fasta

    # For haplotype 2
    formatted_label=$(printf "%s_%s_%s_hap2" "$sample" "$Copy_number" "$Population")
    sed -e 's/\r//g' -e "s/^>.*/>${formatted_label}/" ${gene}_region/${sample}_region_hap2.fasta > ${gene}_region/${sample}_region_relabelled_hap2.fasta

    rm ${gene}_region/AMY_region_chr1_${sample}.vcf.gz.tbi
    rm ${gene}_region/AMY_region_chr1_${sample}.vcf.gz
    rm${gene}_region/ ${gene}_region_chr1_${sample}_excluded_nonvariants.vcf.gz
    rm ${gene}_region/${gene}_region_chr1_${sample}_excluded_nonvariants.vcf.gz.tbi
    rm ${gene}_region/${sample}_region_hap1.fasta
    rm ${gene}_region/${sample}_region_hap1.fasta.fai
    rm ${gene}_region/${sample}_region_hap2.fasta
    rm ${gene}_region/${sample}_region_hap2.fasta.fai

done

cat AMY_region/*.fasta > amy_trees/LWK_Quechua_Maya_no_missing_biallelic.fasta

conda activate iqtree
iqtree -s amy_trees/LWK_Quechua_Maya_no_missing_biallelic.fasta --seqtype DNA -m TEST --alrt 1000 -bb 1000 -bnni -allnni -nt AUTO -pre amy_trees/LWK_Quechua_Maya_haploid_amy_region_within_recombination_no_missing

