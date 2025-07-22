module load gcc/11.2.0
module load clang

conda activate gsl

Info=$1

./selscan --xpehh \
          --vcf filtered_genetic_map_AllQuechua_no_EUR_PhasedChr1_no_indels_0.05_no_missing_unique.vcf.gz \
          --vcf-ref filtered_genetic_map_Maya_PhasedChr1_no_indels_0.05_no_missing_unique.vcf.gz \
          --out Quechua_Maya_${Info} \
          --max-extend 2000000 \ # extend because the SNPs in LD with each other in the amylase locus are pretty far away from each other
          --wagh \
          --map Quechua_Maya_chr1_bedtools_test_genetic_map_no_header.txt  \ # genetic map created from "chr1_avg.bedgraph" file which is from Halldorsson et al 2019 recombination graphs
          --threads 8

./norm --xpehh --files Quechua_Maya_${Info}.xpehh.out
