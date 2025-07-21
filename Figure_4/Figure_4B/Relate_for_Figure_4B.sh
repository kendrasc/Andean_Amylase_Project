module load gcc/11.2.0
module load openmpi/4.1.1
module load r/4.2.0

test="CLUES"

${PATH_TO_RELATE}/bin/RelateFileFormats \
    --mode ConvertFromVcf \
    --haps Quechua_PhasedChr1_no_missing_biallelic_no_EUR.haps \
    --sample Quechua_PhasedChr1_no_missing_biallelic_no_EUR.sample \
    -i Quechua_PhasedChr1_no_missing_biallelic_no_EUR


${PATH_TO_RELATE}/bin/Relate --mode All \
	--haps Quechua_PhasedChr1_no_missing_biallelic_no_EUR.haps \
	--sample Quechua_PhasedChr1_no_missing_biallelic_no_EUR.sample \
	--map ${PATH_TO_FILES}/genetic_map_chr1.txt \
	-N 30000 \
	-m 1.25e-8 \
	-o Quechua_PhasedChr1_no_missing_biallelic_no_EUR_${test} \
	--seed 1

${PATH_TO_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
	-i Quechua_PhasedChr1_no_missing_biallelic_no_EUR_${test} \
	-o Quechua_PhasedChr1_no_missing_biallelic_no_EUR_bypop_${test} \
	-m 1.25e-8 \
	--poplabels Quechua_no_EUR.poplabels \
	--years_per_gen 28 \
	--bins 3,7,0.2 \
	--threads 2 \
	--seed 1


${PATH_TO_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i Quechua_PhasedChr1_no_missing_biallelic_no_EUR_${test} \
        -o Quechua_PhasedChr1_no_missing_biallelic_no_EUR_bypop_sampled_${test}_single_snp \
        -m 1.25e-8 \
	--coal Quechua_PhasedChr1_no_missing_biallelic_no_EUR_bypop_${test}.coal \
        --format b \
        --num_samples 100 \
        --first_bp 103614521 \
	--last_bp 103614521



