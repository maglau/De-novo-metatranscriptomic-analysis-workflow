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

# touch step2_cRNA.log

clear
echo "---Metatranscriptomic workflow step 2 script v.1.0---" > step2_cRNA.log

echo "" >> step2_cRNA.log
echo "Execute Trinity" >> step2_cRNA.log
echo "/home/maglau/tools/trinityrnaseq_r20131110/Trinity.pl --seqType fa --single reads_cRNA.fasta --output reads_cRNA --JM 1G --SS_lib_type F --CPU 16" >> step2_cRNA.log
/home/maglau/tools/trinityrnaseq_r20131110/Trinity.pl --seqType fa --single reads_cRNA.fasta --output cRNA --JM 1G --SS_lib_type F --CPU 16

echo "" >> step2_cRNA.log
echo "python Calculate_avReadLen.py reads_cRNA.fasta" >> step2_cRNA.log
python Calculate_avReadLen.py reads_cRNA.fasta

