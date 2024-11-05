#Still editing


prefix="HG01946_50kb_contig_hap1"
contig="HG01946_hap1_100kb_contig_h1tg000119l:101760897-101966413_amylocus.fa"
partitions="/projects/academic/omergokc/Luane/amy_hap/partitions.fa"

nucmer  --prefix $prefix $contig $partitions 

delta-filter -q -l 10000 ${prefix}.delta > ${prefix}.filter

mummerplot -prefix $prefix --png --large -R $contig  -Q $partitions ${prefix}.filter
