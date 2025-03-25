module load gcc/11.2.0
module load openmpi/4.1.1
module load r/4.2.0

PATH_TO_RELATE="/projects/academic/omergokc/Kendra/Relate/relate"
PATH_TO_FILES="/projects/academic/omergokc/Kendra/Relate/vcf_files"
test="1k_map"

${PATH_TO_RELATE}/bin/RelateFileFormats \
                 --mode ConvertFromVcf \
                 --haps AllQuechua_Maya_Chr1_no_EUR_no_missing.haps \
                 --sample AllQuechua_Maya_Chr1_no_EUR_no_missing.sample \
                 -i AllQuechua_Maya_Chr1_no_EUR_no_missing

${PATH_TO_RELATE}/bin/Relate --mode All \
	--haps AllQuechua_Maya_Chr1_no_EUR_no_missing.haps \
	--sample AllQuechua_Maya_Chr1_no_EUR_no_missing.sample \
	--map ${PATH_TO_FILES}/genetic_map_chr1.txt \
	-N 30000 \
	-m 1.25e-8 \
	-o AllQuechua_Maya_Chr1_no_EUR_no_missing_${test} \
	--seed 1

${PATH_TO_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
	-i AllQuechua_Maya_Chr1_no_EUR_no_missing_${test} \
	-o AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_${test} \
	-m 1.25e-8 \
	--poplabels ${PATH_TO_FILES}/Quechua_Maya_no_EUR_shared_snps_no_indels_0.05_no_missing_unique.poplabels \
	--years_per_gen 28 \
        --bins 3,7,0.2 \
        --threads 2 \
	--seed 1


${PATH_TO_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
  -i AllQuechua_Maya_Chr1_no_EUR_no_missing_${test} \
  -o AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_mrate_${test} \
  -m 1.25e-8 \
  --poplabels ${PATH_TO_FILES}/Quechua_Maya_no_EUR_shared_snps_no_indels_0.05_no_missing_unique.poplabels \
  --norm_mutrate 1 \
  --years_per_gen 28 \
  --bins 3,7,0.2 \
  --threads 2 \
  --seed 1

${PATH_TO_RELATE}/scripts/SampleBranchLengths/ReEstimateBranchLengths.sh \
	-i AllQuechua_Maya_Chr1_no_EUR_no_missing_${test} \
	-o AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_${test} \
	-m 1.25e-8 \
	--coal AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop.coal \
	--seed 1

${PATH_TO_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
	-i AllQuechua_Maya_Chr1_no_EUR_no_missing_${test} \
	-o AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_sampled_${test} \
	-m 1.25e-8 \
	--coal AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop.coal \
	--first_bp 103347080 \
	--last_bp 103831082 \
	--seed 1 \
	--num_samples 100

${PATH_TO_RELATE}/scripts/TreeView/TreeViewSample.sh \
	--haps AllQuechua_Maya_Chr1_no_EUR_no_missing.haps \
	--sample AllQuechua_Maya_Chr1_no_EUR_no_missing.sample \
	--anc AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_sampled_${test}.anc \
	--mut AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_sampled_${test}.mut \
	--bp_of_interest 103614521 \
	--years_per_gen 28 \
	--poplabels ${PATH_TO_FILES}/Quechua_Maya_no_EUR_shared_snps_no_indels_0.05_no_missing_unique.poplabels \
	--dist AllQuechua_Maya_Chr1_no_EUR_no_missing_bypop_sampled_${test}.dist \
	-o 103614521_plot_${test}
