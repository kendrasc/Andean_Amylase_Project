################################# STEP 1: MAPPING #################################

#!/bin/bash
#SBATCH --qos=omergokc
#SBATCH --partition=omergokc
#SBATCH --cluster=faculty
#SBATCH --account=omergokc
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=250G
#SBATCH --job-name=map_reads_all_haplotypes
#SBATCH --output=map_allreads_all_haplotypes_%A_%a.out
#SBATCH --error=map_allreads_all_haplotypes_%A_%a.err
#SBATCH --export=NONE
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-2 # job array index

# Read the assembly name and path from the file
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" reads.txt)
outdir=$(echo "$line" | awk '{print $1}')
reads_hap1=$(echo "$line" | awk '{print $2}')
reads_hap2=$(echo "$line" | awk '{print $3}')

# Load conda environment
eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate herro
haplotype_dir="/projects/academic/omergokc/Luane/amy_hap/haplotypes"
THREADS=32

# Create output directories
mkdir -p ${outdir}_hap1 ${outdir}_hap2

# Loop through each haplotype reference file
for ref in ${haplotype_dir}/*.fa*; do
    # Extract the base name of the reference file (e.g., amy_H10 from amy_H10.fasta)
    sample=$(basename "${ref}" | sed 's/\.[^.]*$//')

    echo "Processing haplotype reference: ${sample}"

    # Map reads for hap1
    minimap2 -ax map-ont -t ${THREADS} --secondary=no "${ref}" "${reads_hap1}" | \
    samtools sort -@ ${THREADS} -m 4G -o "${outdir}_hap1/${sample}_${outdir}_mapped.bam"
    samtools index -@ ${THREADS} "${outdir}_hap1/${sample}_hap1_mapped.bam"

    # Map reads for hap2
    minimap2 -ax map-ont -t ${THREADS} --secondary=no "${ref}" "${reads_hap2}" | \
    samtools sort -@ ${THREADS} -m 4G -o "${sample}_hap2/${outdir}_hap2_mapped.bam"
    samtools index -@ ${THREADS} "${sample}_hap2/${outdir}_${outdir}_mapped.bam"
done

echo "Mapping completed for all haplotypes."


################################# STEP 2: CALCULATING READ DEPTH ################################################
#First, I am filtering for reads with alignments that are bigger than 5kb and generating filtered bam files
for bam in *_hap1/*.bam *_hap2/*.bam; do
    # Get the directory of the current BAM file
    dir=$(dirname "${bam}")
    # Get the base name of the BAM file without the extension
    sample=$(basename "${bam}" .bam)
    echo "Filtering ${sample} for alignments ≥ 5kb..."

    # Define the output BAM file path in the same directory
    filtered_bam="${dir}/filtered_${sample}.bam"

    # Filter alignments using samtools and awk
    samtools view -h "${bam}" | \
    awk 'BEGIN {OFS="\t"}
        /^@/ {print}
        !/^@/ {
            cigar = $6;
            total_match = 0;
            while (match(cigar, /([0-9]+)M/, arr)) {
                total_match += arr[1];
                cigar = substr(cigar, RSTART + RLENGTH);
            }
            if (total_match >= 5000) print
        }' | \
    samtools view -b -o "${filtered_bam}"

    # Index the filtered BAM file
    samtools index "${filtered_bam}"
done

echo "Filtering completed for all BAM files."

########################Calculating read depth for the filtered bam files:
#!/bin/bash

# Enable nullglob to handle cases where no files match the pattern
shopt -s nullglob

# Define the output file
output_file="haplotype_depth_stats.txt"

# Initialize the output file with a header
echo -e "Sample\tMean_Depth\tStd_Dev\tHaplotype\tReference" > "${output_file}"

# Loop through all *_mapped_depth.txt files in the specified directories
for file in *_hap1/*_mapped_depth.txt *_hap2/*_mapped_depth.txt; do
    # Get the directory and extract the haplotype information
    dir=$(dirname "${file}")
    dir_basename=$(basename "${dir}")

    # Extract haplotype from directory name
    if [[ "${dir_basename}" == "hap1" || "${dir_basename}" == "hap2" ]]; then
        haplotype="${dir_basename}"
    elif [[ "${dir_basename}" == *"_hap1" ]]; then
        haplotype="hap1"
    elif [[ "${dir_basename}" == *"_hap2" ]]; then
        haplotype="hap2"
    else
        haplotype="unknown"
    fi

    # Get the filename without the directory path and extension
    filename=$(basename "${file}" _mapped_depth.txt)
    # Remove the "filtered_" prefix
    filename=${filename#filtered_}

    # Extract the reference and sample from the filename
    # Assuming the format is "ref_sample", where ref can contain underscores
    ref="${filename%_*}"    # Reference is everything before the last underscore
    sample="${filename##*_}" # Sample is everything after the last underscore

    echo "Processing Sample: ${sample}, Haplotype: ${haplotype}, Reference: ${ref}..."

    # Use awk to calculate the mean and standard deviation
    awk '{
        sum += $3;         # Sum of depth values (column 3)
        sumsq += $3*$3;    # Sum of squares of depth values
        n++;
    }
    END {
        if (n > 0) {
            mean = sum / n;
            stddev = sqrt(sumsq / n - (mean * mean));
            printf "%s\t%.2f\t%.2f\t%s\t%s\n", "'${sample}'", mean, stddev, "'${haplotype}'", "'${ref}'";
        } else {
            printf "%s\tNA\tNA\t%s\t%s\n", "'${sample}'", "'${haplotype}'", "'${ref}'";
        }
    }' "${file}" >> "${output_file}"
done

echo "Mean read depth and standard deviation calculated for all haplotypes."
echo "Results saved in ${output_file}."

