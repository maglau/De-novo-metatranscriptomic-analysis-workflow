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

# touch step4_ncRNA.log
clear
echo "---Metagenomics workflow step 3 script v.1.0---" > step4_ncRNA.log

echo "" >> step4_ncRNA.log
echo "Calculating read length 'n'" >> step4_ncRNA.log
echo "n=$(<'read_length.txt')" >> step4_ncRNA.log
n=$(<'read_length.txt')
echo "$n"

echo "" >> step4_ncRNA.log
echo "Calculating statistics" >> step4_ncRNA.log
echo "/home/maglau/tools/samtools-1.2/samtools view -bS ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sam > ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.bam" >> step4_ncRNA.log
echo "/home/maglau/tools/samtools-1.2/samtools sort ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.bam ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sorted" >> step4_ncRNA.log
echo "/home/maglau/tools/samtools-1.2/samtools index ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sorted.bam" >> step4_ncRNA.log
echo "/home/maglau/tools/bedtools-2.17.0/bin/bedtools bamtobed -i ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sorted.bam > ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_BT2vs_mappedReads.bed" >> step4_ncRNA.log
echo "python Calculate_stats_v3.py ncRNA/transabyss-merged-reads_nc_new.fa $n ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_BT2vs_" >> step4_ncRNA.log

/home/maglau/tools/samtools-1.2/samtools view -bS ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sam > ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.bam
/home/maglau/tools/samtools-1.2/samtools sort ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.bam ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sorted
/home/maglau/tools/samtools-1.2/samtools index ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sorted.bam
/home/maglau/tools/bedtools-2.17.0/bin/bedtools bamtobed -i ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sorted.bam > ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_BT2vs_mappedReads.bed
python Calculate_stats_v3.py ncRNA/transabyss-merged-reads_nc_new.fa $n ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_BT2vs_

echo "" >> step4_ncRNA.log
echo "Running CD-Hit" >> step4_ncRNA.log
# for ncRNA nucleotide sequences
echo "/home/maglau/tools/cd-hit-v4.6.4-2015-0603/cd-hit-est -i ncRNA/transabyss-merged-reads_nc_new.fa -o cd-hit/reads_ncRNA_95G -c 0.95 -n 10 -g 1 -G 0 -aS 0.8 -d 0 -T 16 > cd-hit/reads_ncRNA_95G.log" >> step4_ncRNA.log
/home/maglau/tools/cd-hit-v4.6.4-2015-0603/cd-hit-est -i ncRNA/transabyss-merged-reads_nc_new.fa -o cd-hit/reads_ncRNA_95G -c 0.95 -n 10 -g 1 -G 0 -aS 0.8 -d 0 -T 16 > cd-hit/reads_ncRNA_95G.log



