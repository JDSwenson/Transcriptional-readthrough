#BSUB -q long
#BSUB -W 20:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e logs/01_sortrna.err
#BSUB -oo logs/01_sortrna.log

module load sortmerna/4.2.0

#Step II.1: Filter rRNA
#### II.1.1) Unzip files (sortmerna does not accept zipped files as input)
gunzip 00_cat_reads/*

wait

#### II.1.2) Align each sample against a database of eukaryotic 18s and 28s rRNA sequences
### First, download sortmerna and the associated eukaryotic rRNA databases
###Check if the program has already been downloaded and remove all associated files if it has
if test -f ./references/sortmerna; then
rm -r ./references/sortmerna
fi

git clone https://github.com/biocore/sortmerna.git references/sortmerna

wait

#Check if sortmerna output already exists and remove all associated files if it has
if test -f ./z_rRNA_reads/idx; then
rm ./z_rRNA_reads/idx/* ./z_rRNA_reads/kvdb/* ./z_rRNA_reads/out/*
rmdir ./z_rRNA_reads/idx/ ./z_rRNA_reads/kvdb/ ./z_rRNA_reads/out/
fi

for file in 00_cat_reads/*.fastq

do

sample=${file%%.fastq} #Store path and basename of file minus read info
out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output

sortmerna -ref references/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta -ref references/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta -reads $file -fastx -other -workdir z_rRNA_reads

wait

#### II.1.3) Move non-rRNA reads and rRNA reads to appropriate folders
mv z_rRNA_reads/out/aligned.log ./metrics/${out_name}.rRNA.Log
mv z_rRNA_reads/out/other.fastq ./01_rRNA_filtered/${out_name}.fastq
mv z_rRNA_reads/out/aligned.fastq ./z_rRNA_reads/${out_name}.rRNA.Aligned.fastq
rm ~/sortmerna/run/*/*

done

wait

#### II.1.4) Create a single summary file with percent of reads that did and did not align to rRNA

##### Make sure file of metrics doesn't already exist. If it does, create it
if test -f ./metrics/line_counts.txt; then
rm ./metrics/line_counts.txt
fi

touch metrics/all_rRNA_metrics.txt

for f in metrics/*.rRNA.Log
do

echo "rRNA aligned reads for sample ${f%%.rRNA.Log}" >> metrics/all_rRNA_metrics.txt
grep "E-value threshold" $f | awk -F" " '{print $3, $8}' | cut -d ')' -f1 | awk -F"(" '{print $1, $2}' | awk -F" " '{print $1, $2}' | sed "s/passing/aligned_to_rRNA/; s/failing/did_NOT_align_to_rRNA/; s/ /:\t/" >> metrics/all_rRNA_metrics.txt

done