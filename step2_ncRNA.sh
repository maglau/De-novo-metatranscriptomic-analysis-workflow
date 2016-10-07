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

touch step2_ncRNA.log

clear
echo "---Metatransciptomic workflow step 1 script v.1.0---" > step2_ncRNA.log

echo "" >> step2_ncRNA.log
echo "module load openmpi" >> step2_ncRNA.log
module load openmpi

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_61 --kmer 61 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_61 --kmer 61 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_65 --kmer 65 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_65 --kmer 65 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_69 --kmer 69 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_69 --kmer 69 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_73 --kmer 73 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_73 --kmer 73 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_77 --kmer 77 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_77 --kmer 77 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_81 --kmer 81 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se reads_ncRNA.fasta --outdir ncRNA/reads_nc_81 --kmer 81 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_85 --kmer 85 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_85 --kmer 85 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_89 --kmer 89 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_89 --kmer 89 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_93 --kmer 93 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_93 --kmer 93 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_97 --kmer 97 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_97 --kmer 97 --threads 16 --SS

echo "" >> step2_ncRNA.log
echo "/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_101 --kmer 101 --threads 16 --SS" >> step2_ncRNA.log
/home/youmio/tools/transabyss-1.5.2/transabyss --se $1 --outdir ncRNA/reads_nc_101 --kmer 101 --threads 16 --SS
