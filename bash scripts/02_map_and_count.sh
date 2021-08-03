#!/bin/bash
#Map reads to indexed reference genome and produce matrix of counts

#BSUB -q long
#BSUB -W 40:00
#BSUB -R rusage[mem=32000]
#BSUB -n 8
#BSUB -R span[hosts=1]
#BSUB -e logs/02_map_and_count.err
#BSUB -oo logs/02_map_and_count.log

#Load required modules
module load star/2.7.0e
module load gcc/8.1.0
module load samtools/1.9
module load RSEM/1.3.0

#II.2: Map each sample to indexed reference genome using STAR

#### Check if unmapped reads file exists; if it does, delete. If it doesn't, create it.
if test -f ./metrics/unmapped_reads.txt; then
rm ./metrics/unmapped_reads.txt
fi

touch metrics/unmapped_reads.txt

#### II.2.1a) Run STAR to align reads to indexed genome using data filtered for rRNA
for file in 01_rRNA_filtered/*.fastq
do
sample=${file%%.fastq} #Store path and basename of file minus read info
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

STAR --runThreadN 8 --genomeDir 00_genome_index --readFilesIn $file --outFileNamePrefix ./02_mapped_reads/$out_name --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

#### II.2.1b) Run STAR to align reads to indexed genome using data NOT filtered for rRNA
#for file in 00_cat_reads/*.fastq.gz
#do
#sample=${file%%.fastq.gz} #Store path and basename of file minus read info
#out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

#STAR --runThreadN 8 --genomeDir 00_genome_index --readFilesIn $file --outFileNamePrefix ./02_mapped_reads/$out_name --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

#### II.2.2) Report unmapped reads for all samples in one file
echo "sample $out_name" >> ./metrics/unmapped_reads.txt
grep "unmapped" ./02_mapped_reads/${out_name}Log.final.out >> ./metrics\
/unmapped_reads.txt

done

wait

#### II.2.3) Move alignments to appropriate folders
mv ./02_mapped_reads/*.sortedByCoord* ./02_mapped_reads/genome_alignments
mv ./02_mapped_reads/*toTranscriptome* ./02_mapped_reads/transcriptome_alignments


#II.3: Convert files to formats that are more user-friendly and/or can be used with IGV

for file2 in 02_mapped_reads/genome_alignments/*.sortedByCoord.out.bam
do

sample=${file2%%Aligned*} #Store path and basename of file minus read info
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

####II.3.1) Index bam files for use with IGV
samtools index $file2

####II.3.2) Convert bam files to bed files so they're human readable
bedtools bamtobed -i $file2 > ${sample}_sorted.bed

echo "finished with $file2"

done

wait

# II.4: Quantify expression with RSEM

for file3 in 02_mapped_reads/transcriptome_alignments/*.bam

do

sample=${file3%%Aligned.toTranscriptome.out.bam} #Store path and basename of file minus read info
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output


#### II.4.1) Quantify expression for each sample individually
#####Arguments are (after paired-end): bam file, location/prefix of rsem index output files, location/prefix of rsem calulate expression output files (from this run)

rsem-calculate-expression --bam --no-bam-output -p 12 --strandedness reverse $file3 00_genome_index/rsem_TAIR10 03_transcript_quant/$out_name

done

wait

#### II.4.2) Combine counts from all individual samples into a single counts matrix
rsem-generate-data-matrix 03_transcript_quant/*.genes.results > all.counts.matrix

#### Re-zip original files to conserve space
gzip 00_cat_reads/*.fastq