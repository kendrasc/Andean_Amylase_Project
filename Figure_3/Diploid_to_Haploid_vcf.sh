chrm="chr1"
gene="AMY"
location1="chr1:103343429-103968264"
info_file="Andean.txt"

# Create Haploid fasta files for individuals from Quechua and Maya populations that have little to none other ancestry in amylase region
file="Quechua_Missing_EUR_ancestry_in_amylase_Maya_103343_10396_chr1_1_region_snp.vcf"
file_name="${gene}_region_${chrm}"
for sample in `bcftools query -l $file`; do
    bcftools view -Oz -s $sample $file -o ${gene}_region/${file_name}_${sample}.vcf.gz

    sample_id="${sample#*_}"

# Use awk to find the matching sample_id and extract Copy_number and Population
    Copy_number=$(awk -v id="$sample_id" '$1 == id {print $2}' $info_file)
    Population=$(awk -v id="$sample_id" '$1 == id {print $3}' $info_file)

# remove the sites that are the same as the reference
    gatk IndexFeatureFile \
    -I ${gene}_region/${file_name}_${sample}.vcf.gz

   gatk SelectVariants \
    -R ${chrm}_hg38.fasta \
    -V ${gene}_region/${file_name}_${sample}.vcf.gz \
    --select-type-to-include SNP \
    --exclude-non-variants true \
    --O ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz

    gatk IndexFeatureFile \
     -I ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz

    # For haplotype 1
    bcftools consensus -f ${chrm}_hg38_10334_10396.fasta -H 1 ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz > ${gene}_region/${sample_id}_region_hap1.fasta

    # For haplotype 2
    bcftools consensus -f ${chrm}_hg38_10334_10396.fasta -H 2 ${gene}_region/${file_name}_${sample}_excluded_nonvariants.vcf.gz > ${gene}_region/${sample_id}_region_hap2.fasta

    # Ensure the variables are correctly formatted
    formatted_label=$(printf "%s_%s_%s_hap1" "$sample_id" "$Copy_number" "$Population")
    # Relabel the header in each FASTA file for haplotype 1
    sed -e 's/\r//g' -e "s/^>.*/>${formatted_label}/" ${gene}_region/${sample_id}_region_hap1.fasta > ${gene}_region/${sample_id}_region_relabelled_hap1.fasta

    # For haplotype 2
    formatted_label=$(printf "%s_%s_%s_hap2" "$sample_id" "$Copy_number" "$Population")
    sed -e 's/\r//g' -e "s/^>.*/>${formatted_label}/" ${gene}_region/${sample_id}_region_hap2.fasta > ${gene}_region/${sample_id}_region_relabelled_hap2.fasta

    # remove unnecessary files
    rm ${gene}_region/AMY_region_chr1_${sample}.vcf.gz.tbi
    rm ${gene}_region/AMY_region_chr1_${sample}.vcf.gz
    rm ${gene}_region/${sample_id}_region_hap1.fasta
    rm ${gene}_region/${sample_id}_region_hap1.fasta.fai
    rm ${gene}_region/${sample_id}_region_hap2.fasta
    rm ${gene}_region/${sample_id}_region_hap2.fasta.fai

done

# Create a single file for all data
cat ${gene}_region/*.fasta > amy_trees/All_relabeled_all_snps.fasta

# Align fasta files
mafft --auto amy_trees/All_relabeled_all_snps.fasta > amy_trees/All_relabeled_Quechua_Maya_all_snps.fasta

# Create phylogeny
iqtree -s amy_trees/All_relabeled_Quechua_Maya_all_snps.fasta --seqtype DNA -m TEST --alrt 1000 -bb 1000 -nt AUTO -pre amy_trees/Quechua_Maya_haploid_amy_region_and_flanking_all_snps
