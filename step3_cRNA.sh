#!/bin/bash

##### Description
# This small shell script wraps the data generation for R in Albertsen et. al., 2013

##### Needed input files
# Fasta file with all assembled scaffolds (keep the naming as >1, >2 etc): assembly.fa

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git

# touch step3_cRNA.log
clear
echo "---Metagenomics workflow step 3 script v.1.0---" > step3_cRNA.log

echo "" >> step3_cRNA.log
echo "Calculating read length 'n'" >> step3_cRNA.log
echo "n=$(<'read_length.txt')" >> step3_cRNA.log
n=$(<'read_length.txt') 
echo "$n" 

echo "" >> step3_cRNA.log
echo "Building bowtie2 index for cRNA contigs" >> step3_cRNA.log
echo "/home/maglau/tools/bowtie2-2.2.5/bowtie2-build -f cRNA/Trinity.fasta cRNA/reads_cRNA_trinity" >> step3_cRNA.log
/home/maglau/tools/bowtie2-2.2.5/bowtie2-build -f cRNA/Trinity.fasta cRNA/reads_cRNA_trinity 

echo "" >> step3_cRNA.log
echo "Mapping cRNA reads to cRNA contigs (use very sensitive mode, -D 20 -R 3 -N 0 -L 20 -i S,1,0.50; --un <path> to keep the unaligned reads)" >> step3_cRNA.log
echo "/home/maglau/tools/bowtie2-2.2.5/bowtie2 --end-to-end --very-sensitive -p 10 --norc -t -x cRNA/reads_cRNA_trinity -f reads_cRNA.fasta -S cRNA/reads_trinity_mappedReads.sam --un cRNA/reads_trinity_unmapped.fasta" >> step3_cRNA.log
/home/maglau/tools/bowtie2-2.2.5/bowtie2 --end-to-end --very-sensitive -p 10 --norc -t -x cRNA/reads_cRNA_trinity -f reads_cRNA.fasta -S cRNA/reads_trinity_mappedReads.sam --un cRNA/reads_trinity_unmapped.fasta

echo "" >> step3_cRNA.log
echo "Calculating statistics for cRNA contigs" >> step3_cRNA.log
echo "/home/maglau/tools/samtools-1.2/samtools view -bS cRNA/reads_trinity_mappedReads.sam > cRNA/reads_trinity_mappedReads.bam" >> step3_cRNA.log
echo "/home/maglau/tools/samtools-1.2/samtools sort cRNA/reads_trinity_mappedReads.bam cRNA/reads_trinity_mappedReads.sorted" >> step3_cRNA.log
echo "/home/maglau/tools/samtools-1.2/samtools index cRNA/reads_trinity_mappedReads.sorted.bam" >> step3_cRNA.log
echo "/home/maglau/tools/bedtools-2.17.0/bin/bedtools bamtobed -i cRNA/reads_trinity_mappedReads.sorted.bam > cRNA/reads_trinity_mappedReads.bed" >> step3_cRNA.log
echo "python Calculate_stats_v3.py cRNA/Trinity.fasta $n cRNA/reads_trinity" >> step3_cRNA.log
/home/maglau/tools/samtools-1.2/samtools view -bS cRNA/reads_trinity_mappedReads.sam > cRNA/reads_trinity_mappedReads.bam
/home/maglau/tools/samtools-1.2/samtools sort cRNA/reads_trinity_mappedReads.bam cRNA/reads_trinity_mappedReads.sorted
/home/maglau/tools/samtools-1.2/samtools index cRNA/reads_trinity_mappedReads.sorted.bam
/home/maglau/tools/bedtools-2.17.0/bin/bedtools bamtobed -i cRNA/reads_trinity_mappedReads.sorted.bam > cRNA/reads_trinity_mappedReads.bed
python Calculate_stats_v3.py cRNA/Trinity.fasta $n cRNA/reads_trinity


