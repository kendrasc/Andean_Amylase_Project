import msprime
import numpy as np
import csv
import collections
import sys

########################################
# 1) Define generation time & convert years -> generations
########################################
GENERATION_TIME = 25  # years per generation

# Read command-line arguments
Quechua = int(sys.argv[1])
Maya_drop = int(sys.argv[2])
Maya_end = int(sys.argv[3])

# Convert times into generations
T_60000 = 60000 // GENERATION_TIME   # 60,000 yrs -> ~2400 gens
T_16000 = 16000 // GENERATION_TIME   # 16,000 yrs -> ~640 gens
T_SPLIT = 10000 // GENERATION_TIME   # 10,000 yrs -> ~400 gens
T_6000  = 6000 // GENERATION_TIME    # 6,000 yrs  -> ~240 gens
T_1000  = 1000 // GENERATION_TIME    # 1,000 yrs  -> ~40 gens

print(f"T_SPLIT: {T_SPLIT}, T_1000: {T_1000}, T_6000: {T_6000}, T_16000: {T_16000}, T_60000: {T_60000}")

# Effective population sizes
ANCESTRAL_NE = 13600  # Original effective size before bottlenecks
BOTTLENECK1_NE = 1750  # Ne at 60,000 years ago
BOTTLENECK2_NE = 2125  # Ne at 16,000 years ago
Quechua_Maya_Bottleneck_NE = 3000 # Ne at 10,000 years ago

########################################
# 2) Build the Demography Model
########################################
dem = msprime.Demography()

# Add ancestral population with its initial size
dem.add_population(name="Ancestral", initial_size=ANCESTRAL_NE)

# Add bottleneck 1 (60,000 ya -> Ne = 1750)
dem.add_population_parameters_change(
    time=T_60000,
    population="Ancestral",
    initial_size=BOTTLENECK1_NE
)

# Add bottleneck 2 (16,000 ya -> Ne = 2125)
dem.add_population_parameters_change(
    time=T_16000,
    population="Ancestral",
    initial_size=BOTTLENECK2_NE
)

# Add split populations (Quechua and Maya) at T_SPLIT
dem.add_population(name="Quechua", initial_size=Quechua_Maya_Bottleneck_NE)
dem.add_population(name="Maya", initial_size=Quechua_Maya_Bottleneck_NE)
dem.add_population_split(
    time=T_SPLIT,
    derived=["Quechua", "Maya"],
    ancestral="Ancestral"
)

# Quechua: Exponential decline to 1500 from [T_SPLIT, T_1000]
r_quechua = np.log(Quechua / Quechua_Maya_Bottleneck_NE) / (T_1000 - T_SPLIT)
dem.add_population_parameters_change(time=T_SPLIT - 1, population="Quechua", initial_size=Quechua_Maya_Bottleneck_NE, growth_rate=r_quechua)
dem.add_population_parameters_change(time=T_1000, population="Quechua", initial_size=Quechua, growth_rate=0)

# Maya: Two growth phases
r_maya_1 = np.log(Maya_drop / Quechua_Maya_Bottleneck_NE) / (T_6000 - T_SPLIT)
dem.add_population_parameters_change(time=T_SPLIT - 1, population="Maya", initial_size=Quechua_Maya_Bottleneck_NE, growth_rate=r_maya_1)
dem.add_population_parameters_change(time=T_6000, population="Maya", initial_size=Maya_drop)
r_maya_2 = np.log(Maya_end / Maya_drop) / (T_1000 - T_6000)
dem.add_population_parameters_change(time=T_6000, population="Maya", growth_rate=r_maya_2)
dem.add_population_parameters_change(time=T_1000, population="Maya", initial_size=Maya_end, growth_rate=0)

# Ensure all events are sorted
dem.sort_events()

########################################
# 3) Simulate ancestry
########################################
samples = {"Quechua": 78, "Maya": 48}
ts_anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=3e9, random_seed=42)

########################################
# 4) Overlay neutral mutations
########################################
ts_mut = msprime.sim_mutations(ts_anc, rate=1.165e-8, random_seed=43)

########################################
# 5) Frequency Analysis & CSV Export
########################################
output_file = "simulation_mutations.csv"

# Map sample node IDs to indices in the genotype array
sample_id_to_index = {node: idx for idx, node in enumerate(ts_mut.samples())}
quechua_samples = [sample_id_to_index[node] for node in ts_mut.samples() if ts_mut.node(node).population == 1]
maya_samples = [sample_id_to_index[node] for node in ts_mut.samples() if ts_mut.node(node).population == 2]

for i, variant in enumerate(ts_mut.variants()):
    print(f"Variant {i+1}: Position={variant.site.position}, Genotypes={variant.genotypes[:10]}")
    if i == 5:  # Stop after 5 variants
        break

output_file = f"simulation_{Quechua}_{Maya_drop}_{Maya_end}.csv"

with open(output_file, mode="w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow([
        "Position",
        "AlleleIndex",
        "Count_Quechua",
        "Count_Maya",
        "Number_Quechua_Samples",
        "Number_Maya_Samples",
        "Frequency_Quechua",
        "Frequency_Maya"
    ])  # Header

    for variant in ts_mut.variants():
        genotypes = variant.genotypes
        q_counter = collections.Counter(genotypes[i] for i in quechua_samples)
        m_counter = collections.Counter(genotypes[i] for i in maya_samples)

        total_Q = len(quechua_samples)
        total_M = len(maya_samples)

        # Gather all possible allele indices that appear in either population
        all_alleles = set(q_counter.keys()).union(m_counter.keys())
        for allele_index in sorted(all_alleles):
            # Skip ancestral allele (index 0)
            if allele_index == 0:
                continue

            count_Q = q_counter[allele_index]
            count_M = m_counter[allele_index]

            freq_Q = count_Q / total_Q
            freq_M = count_M / total_M

            csv_writer.writerow([
                variant.site.position,
                allele_index,
                count_Q,
                count_M,
                total_Q,
                total_M,
                freq_Q,
                freq_M
            ])

print(f"Per-allele mutation data written to {output_file}")
