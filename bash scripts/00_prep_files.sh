#!/bin/bash

#Load required modules
module load gcc/8.1.0
module load gffread/0.9.8
module load RSEM/1.3.0
module load star/2.7.0e
module load bedops/2.4.14-x86_64
module load fastqc/0.11.5
module load python3/3.5.0_packages/multiqc/1.4

#I.2.1) Download arabidopsis genome & gff file
curl ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas -o references/TAIR10_chr_all.fas

wait

##Change names of chromosomes to match annotation file
sed -i 's/1 CHROMOSOME/Chr1/; s/2 CHROMOSOME/Chr2/; s/3 CHROMOSOME/Chr3/; s/4 CHROMOSOME/Chr4/; s/5 CHROMOSOME/Chr5/; s/mitochondria CHROMOSOME/ChrM/; s/chloroplast CHROMOSOME/ChrC/' references/TAIR10_chr_all.fas

##Gene annotation
curl https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.Mar92021.gff -o references/araport11_Chr.gff

wait

#I.2.2) Index genome file
hisat2-build TAIR10_ChrAll.fas 00_genome_index/Tair10idx

#I.2.3) Prepare annotation file for downstream use with bedtools
##Convert gff annotation file to bed file
gff2bed < ./references/araport11_Chr.gff > ./references/araport11_genes.bed

wait

##Extract boundaries for full genes
grep "gene" ./references/araport11_genes.bed > ./references/gene_boundaries.bed

#I.2.4) Concatenate reads
#### Will need to manually adjust input below to adjust for number of lanes. The present script assumes four lanes with the same prefix and single-end reads.
#### Setup: check if files created in the script below already exist. If they do, delete previous versions.

if test -f ./metrics/line_counts.txt; then
rm ./metrics/line_counts.txt
fi

##Concatenate reads from each sample####
for lane in raw_reads/*L001_R1*.fastq.gz
do

sample=${lane%%_L001*} #Store path and basename of file minus read info
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

cat ${sample}_L001_R1_001.fastq.gz ${sample}_L002_R1_001.fastq.gz ${sample}_L003_R1_001.fastq.gz ${sample}_L004_R1_001.fastq.gz > 00_cat_reads/${out_name}.fastq.gz

echo "Finished with sample $sample"

done

wait

# I.2.5) Count lines in raw files and concatenated files to ensure they match (ie no data loss)
touch ./metrics/line_counts.txt #Create text file to store line counts

for file in raw_reads/*L001_R1*.fastq.gz

do

sample1=${file%%_L001*} #Store path and basename of file minus read info
out_name1=${sample1##*/} #Remove path, but keep basename top be used as a prefix to output

echo "Line count for Sample $out_name1" >> ./metrics/line_counts.txt

####Make a variable for each lane
f1=${sample1}_L001_R1*
f2=${sample1}_L002_R1*
f3=${sample1}_L003_R1*
f4=${sample1}_L004_R1*

####Add each of the variables specified above to the wc -l command
linesF=$(wc -l $f1 $f2 $f3 $f4 | awk '/total/ {print $1}')

c1=00_cat_reads/${out_name1}.fastq.gz

linesF_cat=$(wc -l $c1 | awk '{print $1}')

echo "raw reads: $linesF" >> ./metrics/line_counts.txt
echo "concatenated reads: $linesF_cat" >> ./metrics/line_counts.txt

done

wait

# I.2.6) QC individual files with FastQC
for file in 00_cat_reads/*.fastq.gz

do

echo $file

fastqc $file -o logs/fastqc_raw -t 4

done

wait

# I.2.7) Combine QC reports for all files with MultiQC
multiqc logs/fastqc_raw/ -o logs/multiqc_raw -f