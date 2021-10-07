#!/bin/bash

#II.2: Map reads to indexed genome using ht-seq
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/h1-4-H1.fastq -S h1-4-H1.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/h1-4-H2.fastq -S h1-4-H2.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/h1-4-H3.fastq -S h1-4-H3.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/h1-4-RT1.fastq -S h1-4-RT1.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/h1-4-RT2.fastq -S h1-4-RT2.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/h1-4-RT3.fastq -S h1-4-RT3.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/sh2-H1.fastq -S sh2-H1.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/sh2-H2.fastq -S sh2-H2.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/sh2-H3.fastq -S sh2-H3.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/sh2-RT1.fastq -S sh2-RT1.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/sh2-RT2.fastq -S sh2-RT2.sam
hisat2 --dta -q -x Tair10idx -U /project/uma_m2m/SHOT2_RNAseq/sh2-RT3.fastq -S sh2-RT3.sam

wait

#II.3: Convert sam files to bam and bed format
for file in *.sam
do

sample=${file%%.sam} #Store path and basename of file
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

####II.3.1) Convert sam to bam
samtools view -b $file -o ${out_name}.bam

####II.3.2) Index bam files for use with IGV
samtools index ${out_name}.bam

####II.3.3) Convert bam files to bed files so they're human readable
bedtools bamtobed -i ${out_name}.bam > ${sample}.bed

echo "finished with $file2"

done

wait

#II.4: Count reads with ht-seq count
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene h1-4-RT1.bam araport11_Chr.gff > count_hRT1.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene h1-4-RT2.bam araport11_Chr.gff > count_hRT2.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene h1-4-RT3.bam araport11_Chr.gff > count_hRT3.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene sh2-RT1.bam araport11_Chr.gff > count_s2RT1.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene sh2-RT2.bam araport11_Chr.gff > count_s2RT2.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene sh2-RT3.bam araport11_Chr.gff > count_s2RT3.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene h1-4-H1.bam araport11_Chr.gff > count_hH1.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene h1-4-H2.bam araport11_Chr.gff > count_hH2.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene h1-4-H3.bam araport11_Chr.gff > count_hH3.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene sh2-H1.bam araport11_Chr.gff > count_s2H1.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene sh2-H2.bam araport11_Chr.gff > count_s2H2.txt
htseq-count -m intersection-nonempty --stranded reverse -f bam -i ID -t gene sh2-H3.bam araport11_Chr.gff > count_s2H3.txt