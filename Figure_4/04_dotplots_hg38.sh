#!/bin/bash

#This script uses the paf files you have generated with the dotplots scrips. 
#The purpose is to extract certain regions of the contigs from the assembly fasta file, and
#regionally align and visualize it using a reference genome.
#You will need the coordinates of the region of interest from the reference genome. I got it 
#from UCSC. You will also need to extract the fasta file of the region you are interested in. 
#nucmer and mummerplot are part of this package: https://github.com/mummer4/mummer

################################################    Defining variables
Ref="/projects/academic/omergokc/hg38.fa" #Reference genome to generate paf file
Ref_cut="/projects/academic/omergokc/Luane/amy_hap/hg38_103570000-103760000.fa" #Fasta file with the region of the genome you're aiming to align to
assembly_hap1="/projects/academic/omergokc/Luane/HG01922/hifiasm/hifiasm_diploid_ULR/HG01922_diploid_50kb.dip.hap1.p_ctg.fasta" #haplotype 1 of your assembly
assembly_hap2="/projects/academic/omergokc/Luane/HG01922/hifiasm/hifiasm_diploid_ULR/HG01922_diploid_50kb.dip.hap2.p_ctg.fasta" #haplotype 2 of your assembly
rn="hg38" #reference name identifier
#rn="t2t"
outdir="HG01922_hap1_50kb" #outdir for hap1
outdir2="HG01922_hap2_50kb" #outdir for hap2
id="HG01922_50kb" #sample name/identifier

##################Create paf files and dotplots for entire chr and assembly##############
mkdir -p ${outdir}
mkdir -p ${outdir2}

minimap2 -x asm5 ${Ref} ${assembly_hap1} > ${outdir}/${rn}_vs_${id}_hap1.paf
minimap2 -x asm5 ${Ref} ${assembly_hap2} > ${outdir2}/${rn}_vs_${id}_hap2.paf

#to create one dotplot for the entire paf file
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_vs_${id}_hap1.paf #haplotype 1
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir2}/${rn}_vs_${id}_hap2.paf #haplotype 2

#only for chr1, because if not it'll be grouped with all the chr that start with 'chr1'
#I am interested in chr1
awk -F'\t' '$6 == "chr1"' ${outdir}/${rn}_vs_${id}_hap1.paf > ${outdir}/${rn}_${id}_chr1_only.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_${id}_chr1_only.paf

awk -F'\t' '$6 == "chr1"' ${outdir2}/${rn}_vs_${id}_hap2.paf > ${outdir2}/${rn}_${id}_chr1_only.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir2}/${rn}_${id}_chr1_only.paf

##########################################################################################
cd $outdir

#This gets all the contigs that cover the region of interest. In a paf file, the 8th column is the start in the reference 
#and the 9th column is the end. So I am extracting from the whole paf file only the region that aligns to my region of interest. 
#The following steps is to filter the dotplots, because many contigs can be pulled from this. So I developped a way to filter that. 
#Beware that the contigs can start way earlier in the chr1 and end way later. Or they can start/end inside the region -- very likely.
awk '($8 <= 103760000) && ($9 >= 103570000)' ${rn}_${id}_chr1_only.paf > all_contigs_103570000-103760000.paf

#sorting first by the start position, then by contig name. So I have all contigs aligned in which the first ones to show will be the ones that start earlier
#This step is optional in case you want to check the paf file late, this is easier to visualize.
sort -k3 all_contigs_103570000-103760000.paf | sort -k1 > all_contigs_103570000-103760000_sorted.paf

#getting only the start and end of all contigs and the type of strand alignment (+ or -).
cut -f1,3,4,5 all_contigs_103570000-103760000_sorted.paf > all_contigs_103570000-103760000_sorted_start_end.paf

#This step is extracting the parts of the contigs I have on the paf file and saving them to a fasta file.
while read -r i start end strand; do
  samtools faidx $assembly_hap1 "${i}:${start}-${end}" > "${i}_${start}_${end}_amylocus.fa"
done < all_contigs_103570000-103760000_sorted_start_end.paf

concatenating all contigs in one fasta file
cat *_amylocus.fa > allcontigs_amylocus.fa

#nucmer realigns only the selected contigs (extracted) to the extracted reference fasta.
nucmer  --prefix amylocus_${id}_hap1 $Ref_cut allcontigs_amylocus.fa

#Delta-filter filters the alignment in many ways. In this case, the "-q" filters for best alignment by contig. 
#This means for all the possible alignments, it will filter by best identity and larger lenght of the alignment, 
#giving you only one alignment per contig. The -l 1000 will give you alignments of at least 1000bp. Modify this as 
#wished based on the region you're capturing.
delta-filter -q -l 5000 amylocus_${id}_hap1.delta > amylocus_${id}_hap1.filter

#This next step is necessary because you will re-xtract the contigs based on the alignments made by nucmer. 
#This command shows the coordinates of your alignment, so you can get the filtered alignment and re-extract the contigs.
show-coords amylocus_${id}_hap1.filter > amylocus_${id}_hap1.filter.coords.txt

#This command is made to get the specific coordinates on the contigs for the best alignments.
#It will get it from the delta file you just filtered, etc.
delta="amylocus_${id}_hap1.filter.coords.txt"
awk 'BEGIN { FS="\\|"; OFS="\t"; } /^\s*$/ || /^={5,}/ || /^NUCMER/ || /^\s*\[S1\]/ || $0 !~ /\|/ || NF < 5 { next; } { for (i = 1; i <= NF; i++) { gsub(/^ +| +$/, "", $i); } split($2, pos2_fields, " +"); S2 = pos2_fields[1] + 0; E2 = pos2_fields[2] + 0; gsub(/^ +| +$/, "", $5); split($5, tag_fields, "\t"); contig_tag = tag_fields[2]; if (contig_tag == "" || contig_tag == "-") { next; } split(contig_tag, contig_parts, ":"); contig_name = contig_parts[1]; coords = contig_parts[2]; split(coords, coord_arr, "-"); contig_start = coord_arr[1] + 0; start = contig_start + (S2 - 1); end = contig_start + (E2 - 1); if (start > end) { temp = start; start = end; end = temp; } print contig_name, start, end; }' "${delta}" > contig_coordinates_alignments_nucmer.txt

#Now I am reextracting the contigs again, so I get only the parts of the contigs I know better align
#to my region of interest. 
while read -r i start end strand; do
  samtools faidx $assembly_hap1 "${i}:${start}-${end}" > "${i}_${start}_${end}_amylocus_nucmer.fa"
done < contig_coordinates_alignments_nucmer.txt

#Same process again. Concatenate them all in one fasta file. 
cat *_nucmer.fa > allcontigs_amylocus_nucmer.fa

#Realign them using nucmer.
nucmer  --prefix amylocus_${id}_hap1_nucmer $Ref_cut allcontigs_amylocus_nucmer.fa
#This step is for only keeping the best alingments and not showing secondary alignments.
delta-filter -q -l 10000 amylocus_${id}_hap1_nucmer.delta > amylocus_${id}_hap1_nucmer.filter

#plotting the final dotplots.
module load gcc gnuplot
mummerplot -prefix amylocus_${id}_hap1_nucmer.filter --png --large -R $Ref_cut -Q allcontigs_amylocus_nucmer.fa amylocus_${id}_hap1_nucmer.filter

#########################################################################################
#This step is doing the same thing for haplotype 2.

cd ../$outdir2

awk '($8 <= 103760000) && ($9 >= 103570000)' ${rn}_${id}_chr1_only.paf > all_contigs_103570000-103760000.paf

sort -k3 all_contigs_103570000-103760000.paf | sort -k1 > all_contigs_103570000-103760000_sorted.paf

cut -f1,3,4,5 all_contigs_103570000-103760000_sorted.paf > all_contigs_103570000-103760000_sorted_start_end.paf

while read -r i start end strand; do
  samtools faidx $assembly_hap2 "${i}:${start}-${end}" > "${i}_${start}_${end}_amylocus.fa"
done < all_contigs_103570000-103760000_sorted_start_end.paf

cat *_amylocus.fa > allcontigs_amylocus.fa

nucmer  --prefix amylocus_${id}_hap2 $Ref_cut allcontigs_amylocus.fa

delta-filter -q -l 5000 amylocus_${id}_hap2.delta > amylocus_${id}_hap2.filter

show-coords amylocus_${id}_hap2.filter > amylocus_${id}_hap2.filter.coords.txt

delta="amylocus_${id}_hap2.filter.coords.txt"
awk 'BEGIN { FS="\\|"; OFS="\t"; } /^\s*$/ || /^={5,}/ || /^NUCMER/ || /^\s*\[S1\]/ || $0 !~ /\|/ || NF < 5 { next; } { for (i = 1; i <= NF; i++) { gsub(/^ +| +$/, "", $i); } split($2, pos2_fields, " +"); S2 = pos2_fields[1] + 0; E2 = pos2_fields[2] + 0; gsub(/^ +| +$/, "", $5); split($5, tag_fields, "\t"); contig_tag = tag_fields[2]; if (contig_tag == "" || contig_tag == "-") { next; } split(contig_tag, contig_parts, ":"); contig_name = contig_parts[1]; coords = contig_parts[2]; split(coords, coord_arr, "-"); contig_start = coord_arr[1] + 0; start = contig_start + (S2 - 1); end = contig_start + (E2 - 1); if (start > end) { temp = start; start = end; end = temp; } print contig_name, start, end; }' "${delta}" > contig_coordinates_alignments_nucmer.txt

while read -r i start end strand; do
  samtools faidx $assembly_hap2 "${i}:${start}-${end}" > "${i}_${start}_${end}_amylocus_nucmer.fa"
done < contig_coordinates_alignments_nucmer.txt

cat *_nucmer.fa > allcontigs_amylocus_nucmer.fa

nucmer  --prefix amylocus_${id}_hap2_nucmer $Ref_cut allcontigs_amylocus_nucmer.fa
delta-filter -q -l 10000 amylocus_${id}_hap2_nucmer.delta > amylocus_${id}_hap2_nucmer.filter

module load gcc gnuplot
mummerplot -prefix amylocus_${id}_hap2_nucmer.filter --png --large -R $Ref_cut -Q allcontigs_amylocus_nucmer.fa amylocus_${id}_hap2_nucmer.filter
###############################Best commands^############################################

###Keep in mind that after these steps, sometimes one contig might have duplicates on the dotplots. So I manually removed the duplicates, and left only the best alignments.
