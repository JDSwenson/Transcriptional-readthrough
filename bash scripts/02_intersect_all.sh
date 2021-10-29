#!/bin/bash
#intersect reads

#Load required modules
module load bedtools/2.29.2

#Set overlap requirement - default = 30% of read must overlap specified region
over=0.3

for file in 01_mapped_reads/*.bed
do

sample=${file%%.bed} #Store path and basename of file minus read info
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

#Intersect with upstream gene region - stranded
bedtools intersect -a references/expressed_genes_upstream_region.bed -b $file -c -S -F $over > 03_bedtools_out/opposite_strand/${out_name}.all_expressed_genes_upstream_genic_region.intersect.bed

#Intersect with downstream readthrough region - stranded
bedtools intersect -a references/expressed_genes_downstream_region.bed -b $file -c -S -F $over > 03_bedtools_out/opposite_strand/${out_name}.all_expressed_genes_downstream_readthrough_region.intersect.bed

done

#Process all files at once instead of in loop. Would need to adapt pipeline for this, so it's commented out, but want to keep it here because it would be cool to do sometime and would ultimately make things more straightforward ... 
# -C is count all files, -names means give them names instead of file numbers
#bedtools intersect -a references/expressed_gene_boundaries_stranded.bed -b 02_mapped_reads/genome_alignments/*sorted.bed -C -s -F 0.3 -names > 04_bedtools_out/all_genes_stranded.intersect.bed