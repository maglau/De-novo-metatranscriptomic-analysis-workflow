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

touch step3_ncRNA.log

clear 
echo "---Metatransciptomic workflow step 2 script v.1.0---" > step3_ncRNA.log

echo "" >> step3_ncRNA.log
echo "module load openmpi" >> step3_ncRNA.log
module load openmpi

echo "" >> step3_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss-merge ncRNA/reads_nc_61/transabyss-final.fa ncRNA/reads_nc_65/transabyss-final.fa ncRNA/reads_nc_69/transabyss-final.fa ncRNA/reads_nc_73/transabyss-final.fa ncRNA/reads_nc_77/transabyss-final.fa ncRNA/reads_nc_81/transabyss-final.fa ncRNA/reads_nc_85/transabyss-final.fa ncRNA/reads_nc_89/transabyss-final.fa ncRNA/reads_nc_93/transabyss-final.fa ncRNA/reads_nc_97/transabyss-final.fa ncRNA/reads_nc_101/transabyss-final.fa --mink 61 --maxk 101 --length 200 --threads 20 --out ncRNA/transabyss-merged-reads_nc_new.fa --prefix k61. k65. k69. k73. k77. k81. k85. k89. k93. k97. k101" >> step3_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss-merge ncRNA/reads_nc_61/transabyss-final.fa ncRNA/reads_nc_65/transabyss-final.fa ncRNA/reads_nc_69/transabyss-final.fa ncRNA/reads_nc_73/transabyss-final.fa ncRNA/reads_nc_77/transabyss-final.fa ncRNA/reads_nc_81/transabyss-final.fa ncRNA/reads_nc_85/transabyss-final.fa ncRNA/reads_nc_89/transabyss-final.fa ncRNA/reads_nc_93/transabyss-final.fa ncRNA/reads_nc_97/transabyss-final.fa ncRNA/reads_nc_101/transabyss-final.fa --mink 61 --maxk 101 --length 200 --threads 20 --out ncRNA/transabyss-merged-reads_nc_new.fa --prefix k61. k65. k69. k73. k77. k81. k85. k89. k93. k97. k101


echo "" >> step3_ncRNA.log
echo "python Calculate_avReadLen.py reads_ncRNA.fasta" >> step3_ncRNA.log
python Calculate_avReadLen.py reads_ncRNA.fasta

echo "" >> step3_ncRNA.log
echo "Calculating read length 'n'" >> step3_ncRNA.log
echo "n=$(<'read_length.txt')" >> step3_ncRNA.log
n=$(<'read_length.txt')
echo "$n"

echo "" >> step3_ncRNA.log
echo "Building bowtie2 index for cRNA contigs" >> step3_ncRNA.log
echo "/home/maglau/tools/bowtie2-2.2.5/bowtie2-build -f ncRNA/transabyss-merged-reads_nc_new.fa ncRNA/reads_ncRNA_transabyssBT2" >> step3_ncRNA.log
/home/maglau/tools/bowtie2-2.2.5/bowtie2-build -f ncRNA/transabyss-merged-reads_nc_new.fa ncRNA/reads_ncRNA_transabyssBT2

echo "" >> step3_ncRNA.log
echo "/home/maglau/tools/bowtie2-2.2.5/bowtie2 --end-to-end --very-sensitive -p 10 --norc -t -x ncRNA/reads_ncRNA_transabyssBT2 -f reads_ncRNA.fasta -S ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sam --un ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_unmapped.fasta"
/home/maglau/tools/bowtie2-2.2.5/bowtie2 --end-to-end --very-sensitive -p 10 --norc -t -x ncRNA/reads_ncRNA_transabyssBT2 -f reads_ncRNA.fasta -S ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_mappedReads_BT2vs.sam --un ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_unmapped.fasta 
