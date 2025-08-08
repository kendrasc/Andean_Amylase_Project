We used the same output files as Figure 4B and Figure S12 for Quechua.

${PATH_TO_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i Quechua_PhasedChr1_no_missing_biallelic_no_EUR_${test} \
        -o Quechua_PhasedChr1_no_missing_biallelic_no_EUR_bypop_sampled_${test}_EGLN1_rs1769792_coal \
        -m 1.25e-8 \
        --coal Quechua_PhasedChr1_no_missing_biallelic_no_EUR_bypop_${test}.coal \
        --format b \
        --num_samples 100 \
        --first_bp 231462872 \
        --last_bp 231462872
