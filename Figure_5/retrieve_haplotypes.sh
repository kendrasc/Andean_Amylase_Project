#!/bin/bash
samples="/projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02004/GCA_042027455.1_HG02004_pat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02004/GCA_042027105.1_HG02004_mat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG01978/GCA_018472865.2_HG01978_pat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG01978/GCA_018472845.2_HG01978_mat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02293/GCA_042026455.1_HG02293_mat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02293/GCA_042027495.1_HG02293_pat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02300/GCA_042028475.1_HG02300_pat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02300/GCA_042027635.1_HG02300_mat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02148/GCA_018471535.2_HG02148_mat_hprc_f2_genomic.fna /projects/academic/omergokc/Luane/PEL_data_HGSV_Sudmant/HG02148/GCA_018471525.2_HG02148_pat_hprc_f2_genomic.fna"

haplotypes="hap1.fa hap2.fa"

for a in $samples; do
    blastn -query up_down.fa -subject $a -outfmt 6 -out ${a}_HPRC_vs_up_down.blastn.out -max_hsps 1 -max_target_seqs 1
    
    awk -v coordfile="${a}_coords.txt" '
    {
      if ($1 == "upstream_amy") {
        up_contig = $2;
        up_start = $9;
        up_end = $10;
      } else if ($1 == "downstream_amy") {
        down_contig = $2;
        down_start = $9;
        down_end = $10;
      }

      if (up_contig == down_contig) {
        contig = up_contig;

        start = (up_start < down_start) ? up_start : down_start;
        end = (up_end > down_end) ? up_end : down_end;

        # reverse coordinates when necessary
        if (start > end) {
          temp = start;
          start = end;
          end = temp;
        }

        printf("%s\t%s\t%s\n", contig, start, end) > coordfile;
      }
    }
    ' ${a}_HPRC_vs_up_down.blastn.out
    
    while IFS=$'\t' read -r contig start end; do
      if [ "$start" -gt "$end" ]; then
        temp=$start
        start=$end
        end=$temp
      fi

      samtools faidx $a "${contig}:${start}-${end}" >> ${a}_HPRC_extracted_sequences.fa
    done < ${a}_coords.txt
  done
done
