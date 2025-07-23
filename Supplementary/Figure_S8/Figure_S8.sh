module load gcc
module load bcftools
module load gatk
module load samtools

chrm="chr1"
gene="AMY_updated"
pop=$1
location1="chr1:103347982-103830994"
info_file="unrelated_${pop}.txt"
file="${pop}_no_missing_103347982_103830994_biallelic.vcf.gz"
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

    rm ${gene}_region/${gene}_region_chr1_${sample}.vcf.gz.tbi
    rm ${gene}_region/${gene}_region_chr1_${sample}.vcf.gz
    rm ${gene}_region/${gene}_region_chr1_${sample}_excluded_nonvariants.vcf.gz
    rm ${gene}_region/${gene}_region_chr1_${sample}_excluded_nonvariants.vcf.gz.tbi
    rm ${gene}_region/${sample}_region_hap1.fasta
    rm ${gene}_region/${sample}_region_hap1.fasta.fai
    rm ${gene}_region/${sample}_region_hap2.fasta
    rm ${gene}_region/${sample}_region_hap2.fasta.fai

done

cat AMY_updated/*.fasta > AMY_updated_EAS_tree/EAS_Quechua_Maya_no_missing_biallelic.fasta


conda activate iqtree

iqtree2  \
-s AMY_updated_EAS_tree/EAS_Quechua_Maya_no_missing_biallelic.fasta  \
-st DNA  \
-m MFP  \
-nt AUTO  \
-pre AMY_updated_EAS_tree/EAS_Quechua_Maya_no_missing_biallelic  \
-alrt 1000 -bb 1000  \
-allnni -bnni  \
-keep-ident

################# invariant sites only
# I got EAS_Quechua_Maya_no_missing_biallelic_corrected_for_ascertainment_bias.varsites.phy
# by running the above script with "+ASC". This fails but then creates this .phy file that only
# includes the invariant sites.

iqtree2 \
  -s AMY_updated_EAS_tree/EAS_Quechua_Maya_no_missing_biallelic_corrected_for_ascertainment_bias.varsites.phy \
  -st DNA \
  -m MFP+ASC \
  -nt AUTO \
  -pre AMY_updated_EAS_tree/EAS_Quechua_Maya_no_missing_biallelic_varsites_invariant_sites_only \
  -alrt 1000 -bb 1000 \
  -allnni -bnni \
  -keep-ident
