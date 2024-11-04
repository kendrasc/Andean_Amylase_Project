Pop1=$1
Pop2=$2

# Filter to just get populations of interest and remove SNPs with a minor allele frequency lower than 5% and missing positions
# VCF file name different for Maya and Quechua samples
bcftools view -S ${Pop1}_${Pop2}_unrelated.txt \
              1kGP_high_coverage_Illumina_chr1_amylase_filtered_SNV_INDEL_SV_phased_panel.vcf.gz \
              -Oz -o 1k_genomes_${Pop1}_${Pop2}_chr1.vcf.gz

bcftools view -i 'MAF>=0.05 & GT!="." & TYPE="snp"' \
              1k_genomes_${Pop1}_${Pop2}_chr1.vcf.gz \
              -Oz -o 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing.vcf.gz

# Step 2: Filter to keep only unique positions
bcftools norm -d all 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing.vcf.gz \
              -Oz -o 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing_unique.vcf.gz

# Step 3: Calculate Fst
vcftools --gzvcf 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing_unique.vcf.gz \
         --min-alleles 2 \
         --max-alleles 2 \
         --weir-fst-pop ${Pop1}_unrelated_sample_IDs.txt \
         --weir-fst-pop ${Pop2}_unrelated_sample_IDs.txt \
         --out ${Pop1}_${Pop2}_amylase_all_chr1_biallelic_0.05_1_region_fst

# Step 4: Convert to Genotype Format
vcftools --gzvcf 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing_unique.vcf.gz \
         --min-alleles 2 \
         --max-alleles 2 \
         --extract-FORMAT-info GT \
         --out 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing_unique_Genotype

mv 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing_unique_Genotype.GT.FORMAT 1k_genomes_${Pop1}_${Pop2}_chr1_no_indels_0.05_no_missing_unique_Genotype.txt
