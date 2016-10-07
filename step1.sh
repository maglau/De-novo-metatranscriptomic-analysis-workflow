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

# touch step1.log

clear
echo "---Metagenomics workflow script v.1.0---" > step1.log

echo "" >> step1.log
echo "Usearch SSU1" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_1.udb -id 0.8 -blast6out reads_SSU1.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_1.udb -id 0.8 -blast6out reads_SSU1.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU2" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_2.udb -id 0.8 -blast6out reads_SSU2.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_2.udb -id 0.8 -blast6out reads_SSU2.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU3" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_3.udb -id 0.8 -blast6out reads_SSU3.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_3.udb -id 0.8 -blast6out reads_SSU3.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU4" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_4.udb -id 0.8 -blast6out reads_SSU4.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_4.udb -id 0.8 -blast6out reads_SSU4.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU5" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_5.udb -id 0.8 -blast6out reads_SSU5.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_5.udb -id 0.8 -blast6out reads_SSU5.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU6" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_6.udb -id 0.8 -blast6out reads_SSU6.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_6.udb -id 0.8 -blast6out reads_SSU6.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU7" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_7.udb -id 0.8 -blast6out reads_SSU7.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_7.udb -id 0.8 -blast6out reads_SSU7.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch SSU8" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_8.udb -id 0.8 -blast6out reads_SSU8.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_8.udb -id 0.8 -blast6out reads_SSU8.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch LSU" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/LSUParc119_ATCG.udb -id 0.8 -blast6out reads_LSU.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/LSUParc119_ATCG.udb -id 0.8 -blast6out reads_LSU.m8 -strand both -threads 12

echo "" >> step1.log
echo "Usearch tRNA and 5SrRNA" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/5SrRNA.udb -id 0.8 -blast6out reads_5SrRNA.m8 -strand both -threads 12" >> step1.log
echo "/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/tRNA.udb -id 0.8 -blast6out reads_tRNA.m8 -strand both -threads 12" >> step1.log
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/5SrRNA.udb -id 0.8 -blast6out reads_5SrRNA.m8 -strand both -threads 12
/home/maglau/tools/usearch -usearch_global $1 -db /scratch/gpfs/maglau/DB/tRNA.udb -id 0.8 -blast6out reads_tRNA.m8 -strand both -threads 12

echo "" >> step1.log
echo "Combine all SSU files into a single file" >> step1.log
echo "cat reads_SSU1.m8 reads_SSU2.m8 reads_SSU3.m8 reads_SSU4.m8 reads_SSU5.m8 reads_SSU6.m8 reads_SSU7.m8 reads_SSU8.m8 > reads_SSU.m8" >> step1.log
cat reads_SSU1.m8 reads_SSU2.m8 reads_SSU3.m8 reads_SSU4.m8 reads_SSU5.m8 reads_SSU6.m8 reads_SSU7.m8 reads_SSU8.m8 > reads_SSU.m8

echo "" >> step1.log
echo "Summarize search results to determine which sequences are related to LSU, SSU, 5SrRNA and tRNA" >> step1.log
echo "python UsearchStat_v3.py reads_LSU.m8 reads_SSU.m8 reads_5SrRNA.m8 reads_tRNA.m8 > reads_usearchsummary.txt" >> step1.log
python UsearchStat_v3.py reads_LSU.m8 reads_SSU.m8 reads_5SrRNA.m8 reads_tRNA.m8 > reads_usearchsummary.txt

echo "" >> step1.log
echo "Separate the original fasta into _cRNA.fasta and _ncRNA.fasta" >> step1.log
echo "python ParseSeqsIntoGroups_v4.py reads_usearchsummary.txt reads.fasta reads" >> step1.log
python ParseSeqsIntoGroups_v4.py reads_usearchsummary.txt reads.fasta reads






