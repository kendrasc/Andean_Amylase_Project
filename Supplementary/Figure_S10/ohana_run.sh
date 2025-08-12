# Remove sites with missing data (keep only complete sites)
vcftools --gzvcf merged_Quechua_IBSYRI_v42.vcf.gz \
  --max-missing 1.0 \
  --recode --recode-INFO-all \
  --out no_missing_merged_Quechua_IBSYRI

# Convert VCF to PLINK PED format (no missing genotypes allowed)
plink --vcf no_missing_merged_Quechua_IBSYRI.recode.vcf \
      --recode 12 tab \
      --geno 0.0 \
      --out merged_chr1_ibsyri_clean

# Convert PED to Ohana DGM format
convert ped2dgm merged_chr1_ibsyri_clean.ped g_full_yriibs_clean.dgm

# Randomly sample 5% of SNPs from the full dataset for Q matrix inference
python ./sample-sites.py ./g_full_yriibs_clean.dgm 5 ./g_yriibs_5percent_clean.dgm

# Estimate Q and F matrices from the full dataset (K=3), using 5% SNP subset
qpas ./g_full_yriibs_clean.dgm -k 3 -qo ./yriibs_5percent_clean_Q.matrix -fo ./yriibs_5percent_clean_F.matrix -e 1e-5 -mi 200

# Generate covariance (C) matrix from the full dataset using F matrix from 5% subset (null hypothesis)
nemeco ./g_full_yriibs_clean.dgm ./yriibs_5percent_clean_F.matrix -e 0.0 -co ./yriibs_5percent_C.matrix -mi 5

# Generate F matrix for full dataset using fixed Q from 5% subset
qpas ./g_full_yriibs_clean.dgm -qi ./yriibs_5percent_clean_Q.matrix -fo ./f_full_ibsyri_clean_5per.matrix -fq -e 1e-5

# Perform selection scan using full dataset, full F matrix, null C matrix, and contrast-specific C matrix
selscan g_full_yriibs_clean.dgm f_full_ibsyri_clean_5per.matrix yriibs_5percent_C.matrix -cs cs.matrix_ibsyri_clean_5per > lle-ratios_ibsyri_clean_5per.txt 

