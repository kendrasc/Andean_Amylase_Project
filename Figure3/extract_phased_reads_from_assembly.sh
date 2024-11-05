#!/bin/bash
######## INPUTS
reads="/projects/academic/omergokc/Luane/HG01946/herro/HG01946_10kbp_HERRO.fasta" #path to corrected reads
raw_reads="/projects/academic/omergokc/Luane/HG01946/HG01946_allreads_Apr16.fastq" #path to raw reads

Ref="/projects/academic/omergokc/hg38.fa" #Reference genome to generate paf file
Ref_cut="/projects/academic/omergokc/Luane/amy_hap/hg38_103570000-103760000.fa" #Fasta file with the region of the genome you're aiming to align to

assembly_hap1="/projects/academic/omergokc/Luane/HG01946/hifiasm/diploid_raw_long_reads/HG01946_diploid_50kb.dip.hap1.p_ctg.fasta" #path to hifiasm assembly, its not the final one! - hap1 
assembly_hap2="/projects/academic/omergokc/Luane/HG01946/hifiasm/diploid_raw_long_reads/HG01946_diploid_50kb.dip.hap2.p_ctg.fasta" #path to hifiasm assembly, its not the final one! - hap2

rn="hg38" #referece genome name
id="HG01946" #sample id/organism id

###################Create paf files and dotplots for entire chr and assembly##############
mkdir -p "${id}_hap1_50kb"
mkdir -p "${id}_hap2_50kb"

minimap2 -x asm5 ${Ref} ${assembly_hap1}.fasta > "${id}_hap1"/${rn}_vs_${id}_hap1.paf
minimap2 -x asm5 ${Ref} ${assembly_hap2}.fasta > "${id}_hap2"/${rn}_vs_${id}_hap2.paf

################################# HAP1 ###################################################
cd "${id}_hap1"
n="1"

#Extracting only chromosome 1
awk -F'\t' '$6 == "chr1"' ${rn}_vs_${id}_hap1.paf > ${rn}_${id}_chr1_only.paf

#What I need to to here is first extract all contigs that align to that region
awk '($8 <= 103760000) && ($9 >= 103570000)' ${rn}_${id}_chr1_only.paf > contigs_chr1_103570000-103760000_hap${n}_${id}.paf

#get only the names of the contigs
cat contigs_chr1_103570000-103760000_hap${n}_${id}.paf | cut -f1 > contigs_chr1_103570000-103760000_hap${n}_${id}_names.txt

#this command will extract only the information for each contig (which contains which reads made up to each contig) without the actual sequence of the contig (which makes the file easier to work with)
awk 'NR==FNR {contigs[$1]; next} 
     ($2 in contigs) && ($1 == "S") {print $1, $2, "."} 
     ($2 in contigs) && ($1 != "S") {print}' contigs_chr1_103570000-103760000_hap${n}_${id}_names.txt $assembly_hap1.gfa > selected_contigs_info_no_seq.gfa

#get rid of the contig names that are in each line with "S"
cut -f5 selected_contigs_info_no_seq.gfa | grep -v "^S" > contigs_chr1_103570000-103760000_hap${n}_${id}_info_no_seq_onlyreads.txt

#get the read name, start, and end in one file from the gfa file you just created
awk -F'[_:]' '{
    split($0, parts, /[_:]/)
    read_name = parts[1]
    for(i=1;i<=NF;i++) {
        if($i ~ /sliding/) {
            split($(i+1), coords, "-")
            start = coords[1]
            end = coords[2]
            print read_name, start, end
            break
        }
    }
}' OFS='\t' contigs_chr1_103570000-103760000_hap${n}_${id}_info_no_seq_onlyreads.txt > read_start_end_contigs_chr1_103570000-103760000_hap${n}_${id}.tsv

#Now I need to extract the reads

#Getting the names of the reads, sorting them and getting only unique read names.
cut -f1 read_start_end_contigs_chr1_103570000-103760000_hap${n}_${id}.tsv | sort | uniq > read_names_contigs_chr1_103570000-103760000_hap${n}_${id}_uniq.txt

#doing the same thing for the contigs
cat contigs_chr1_103570000-103760000_hap${n}_${id}_names.txt | sort | uniq > contigs_chr1_103570000-103760000_hap${n}_${id}_names_uniq.txt

module load gcc seqtk

#Extracting the reads that align to the contigs from the allreads file
seqtk subseq ${reads} read_names_contigs_chr1_103570000-103760000_hap${n}_${id}_uniq.txt > reads_contigs_chr1_103570000-103760000_hap${n}_${id}.fq

#Extracting the contigs that align to the amylase region
seqtk subseq ${assembly_hap1}.fasta contigs_chr1_103570000-103760000_hap${n}_${id}_names_uniq.txt > contigs_chr1_103570000-103760000_hap${n}_${id}.fq

########## Mapping the reads back into the contigs, so I can extract only the reads that make up the region that I want. This step is necessary because some of these contigs actually map to longer regions, and therefore I am not sure if the reads are actually from my region of interest. 

#mapping the reads to my contigs and creating a paf file.
minimap2 -x asm5 contigs_chr1_103570000-103760000_hap${n}_${id}.fq reads_contigs_chr1_103570000-103760000_hap${n}_${id}.fq > ${id}_hap${n}_contig_vs_reads_amyregion.paf

#From this paf file I need to get the reads that map to the coordinates that I actually mapped to the region of interest of hg38. 

#Transforming my paf to bed file only with coordinates for both paf files. 
cut -f 1,3,4 contigs_chr1_103570000-103760000_hap${n}_${id}.paf > contigs_chr1_103570000-103760000_hap${n}_${id}_to_hg38_coord.bed

awk '{print $6"\t"$8"\t"$9"\t"$1}' ${id}_hap${n}_contig_vs_reads_amyregion.paf > reads_to_contigs_chr1_103570000-103760000_hap${n}_${id}_coord.bed

#using bedtools to intersect files. 
module load gcc bedtools

#-wa: Write the original entry from file a (reads) that has any overlap with b.
#-a reads_on_contigs.bed: Our reads mapped to contigs.
#-b contig_regions.bed: Contig regions of interest.
bedtools intersect -a reads_to_contigs_chr1_103570000-103760000_hap${n}_${id}_coord.bed -b contigs_chr1_103570000-103760000_hap${n}_${id}_to_hg38_coord.bed -wa > overlapping_reads.bed

#Getting only the read names that I want to extract
awk '{print $4}' overlapping_reads.bed | sort | uniq > read_names.txt

#Extract the reads I want from the original fastq file. 
seqtk subseq $reads read_names.txt > reads_amylocus_${id}_hap${n}_chr1_103570000-103760000.fq

#Extract the reads I want from the original fastq file. 
seqtk subseq $raw_reads read_names.txt > raw_reads_amylocus_${id}_hap${n}_chr1_103570000-103760000.fq


################################# HAP2 ###################################################
cd "${id}_hap2"
n="2"

#Extracting only chromosome 1
awk -F'\t' '$6 == "chr1"' ${rn}_vs_${id}_hap${n}.paf > ${rn}_${id}_chr1_only.paf

#What I need to to here is first extract all contigs that align to that region
awk '($8 <= 103760000) && ($9 >= 103570000)' ${rn}_${id}_chr1_only.paf > contigs_chr1_103570000-103760000_hap${n}_${id}.paf

#get only the names of the contigs
cat contigs_chr1_103570000-103760000_hap${n}_${id}.paf | cut -f1 > contigs_chr1_103570000-103760000_hap${n}_${id}_names.txt

#this command will extract only the information for each contig (which contains which reads made up to each contig) without the actual sequence of the contig (which makes the file easier to work with)
awk 'NR==FNR {contigs[$1]; next} 
     ($2 in contigs) && ($1 == "S") {print $1, $2, "."} 
     ($2 in contigs) && ($1 != "S") {print}' contigs_chr1_103570000-103760000_hap${n}_${id}_names.txt $assembly_hap2.gfa > selected_contigs_info_no_seq.gfa

#get rid of the contig names that are in each line with "S"
cut -f5 selected_contigs_info_no_seq.gfa | grep -v "^S" > contigs_chr1_103570000-103760000_hap${n}_${id}_info_no_seq_onlyreads.txt

#get the read name, start, and end in one file from the gfa file you just created
awk -F'[_:]' '{
    split($0, parts, /[_:]/)
    read_name = parts[1]
    for(i=1;i<=NF;i++) {
        if($i ~ /sliding/) {
            split($(i+1), coords, "-")
            start = coords[1]
            end = coords[2]
            print read_name, start, end
            break
        }
    }
}' OFS='\t' contigs_chr1_103570000-103760000_hap${n}_${id}_info_no_seq_onlyreads.txt > read_start_end_contigs_chr1_103570000-103760000_hap${n}_${id}.tsv

#Now I need to extract the reads

#Getting the names of the reads, sorting them and getting only unique read names.
cut -f1 read_start_end_contigs_chr1_103570000-103760000_hap${n}_${id}.tsv | sort | uniq > read_names_contigs_chr1_103570000-103760000_hap${n}_${id}_uniq.txt

#doing the same thing for the contigs
cat contigs_chr1_103570000-103760000_hap${n}_${id}_names.txt | sort | uniq > contigs_chr1_103570000-103760000_hap${n}_${id}_names_uniq.txt

module load gcc seqtk

#Extracting the reads that align to the contigs from the allreads file
seqtk subseq ${reads} read_names_contigs_chr1_103570000-103760000_hap${n}_${id}_uniq.txt > reads_contigs_chr1_103570000-103760000_hap${n}_${id}.fq

#Extracting the contigs that align to the amylase region
seqtk subseq ${assembly_hap2}.fasta contigs_chr1_103570000-103760000_hap${n}_${id}_names_uniq.txt > contigs_chr1_103570000-103760000_hap${n}_${id}.fq

########## Mapping the reads back into the contigs, so I can extract only the reads that make up the region that I want. This step is necessary because some of these contigs actually map to longer regions, and therefore I am not sure if the reads are actually from my region of interest. 

#mapping the reads to my contigs and creating a paf file.
minimap2 -x asm5 contigs_chr1_103570000-103760000_hap${n}_${id}.fq reads_contigs_chr1_103570000-103760000_hap${n}_${id}.fq > ${id}_hap${n}_contig_vs_reads_amyregion.paf

#From this paf file I need to get the reads that map to the coordinates that I actually mapped to the region of interest of hg38. 

#Transforming my paf to bed file only with coordinates for both paf files. 
cut -f 1,3,4 contigs_chr1_103570000-103760000_hap${n}_${id}.paf > contigs_chr1_103570000-103760000_hap${n}_${id}_to_hg38_coord.bed

awk '{print $6"\t"$8"\t"$9"\t"$1}' ${id}_hap${n}_contig_vs_reads_amyregion.paf > reads_to_contigs_chr1_103570000-103760000_hap${n}_${id}_coord.bed

#using bedtools to intersect files. 
module load gcc bedtools

#-wa: Write the original entry from file a (reads) that has any overlap with b.
#-a reads_on_contigs.bed: Our reads mapped to contigs.
#-b contig_regions.bed: Contig regions of interest.
bedtools intersect -a reads_to_contigs_chr1_103570000-103760000_hap${n}_${id}_coord.bed -b contigs_chr1_103570000-103760000_hap${n}_${id}_to_hg38_coord.bed -wa > overlapping_reads.bed

#Getting only the read names that I want to extract
awk '{print $4}' overlapping_reads.bed | sort | uniq > read_names.txt

#Extract the reads I want from the original fastq file. 
seqtk subseq $reads read_names.txt > reads_amylocus_${id}_hap${n}_chr1_103570000-103760000.fq

#Extract the reads I want from the original fastq file. 
seqtk subseq $raw_reads read_names.txt > raw_reads_amylocus_${id}_hap${n}_chr1_103570000-103760000.fq

conda deactivate
