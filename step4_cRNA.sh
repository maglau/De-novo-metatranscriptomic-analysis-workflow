#!/bin/bash

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git

# touch step4_cRNA.log

clear
echo "---Metatranscriptomic workflow step 4 script v.1.0---" > step4_cRNA.log

echo "" >> step4_cRNA.log
echo "Running Prodigal on cRNA contigs" >> step4_cRNA.log
echo "/home/maglau/tools/prodigal -i cRNA/Trinity.fasta -a cRNA/Trinity.fasta.prodigal.faa -p meta -g 11 -d cRNA/Trinity.fasta.prodigal.fna -o cRNA/Prodigal_coords.gbk" >> step4_cRNA.log
echo "python /scratch/gpfs/abehmard/PickSenseORFs.py cRNA/Trinity.fasta.prodigal.faa cRNA/Trinity.fasta.prodigal" >> step4_cRNA.log
/home/maglau/tools/prodigal -i cRNA/Trinity.fasta -a cRNA/Trinity.fasta.prodigal.faa -p meta -g 11 -d cRNA/Trinity.fasta.prodigal.fna -o cRNA/Prodigal_coords.gbk
python /scratch/gpfs/abehmard/PickSenseORFs.py cRNA/Trinity.fasta.prodigal.faa cRNA/Trinity.fasta.prodigal

echo "" >> step4_cRNA.log
echo "Running CD-hit on cRNA contigs for cDNA ORFs" >> step4_cRNA.log
echo "/home/maglau/tools/cd-hit-v4.6.4-2015-0603/cd-hit -i ../cRNA/Trinity.fasta.prodigal_sense.faa -o cd-hit/reads_cRNA_sense90G -c 0.9 -n 5 -g 1 -G 0 -aS 0.8 -d 0 -T 16 > cd-hit/reads_cRNA_sense90G.log" >> step4_cRNA.log
/home/maglau/tools/cd-hit-v4.6.4-2015-0603/cd-hit -i cRNA/Trinity.fasta.prodigal_sense.faa -o cd-hit/reads_cRNA_sense90G -c 0.9 -n 5 -g 1 -G 0 -aS 0.8 -d 0 -T 16 > cd-hit/reads_cRNA_sense90G.log


echo "" >> step4_cRNA.log
echo "Processing CD-hit output files for cRNA" >> step4_cRNA.log
echo "python reformat_cdhit_clstr_v2_cRNA_singlesample.py cd-hit/reads_cRNA_sense90G.clstr cd-hit/reads_cRNA_sense90G" >> step4_cRNA.log
python reformat_cdhit_clstr_v2_cRNA_singlesample.py cd-hit/reads_cRNA_sense90G.clstr cd-hit/reads_cRNA_sense90G

echo "" >> step4_cRNA.log
echo "Splits the PEG cluster data into subfiles with 2000 lines each for easier BLASTing"
echo "python Split_FASTA.py cd-hit/reads_cRNA_sense90G.rename 2000"
python Split_FASTA.py cd-hit/reads_cRNA_sense90G.rename 2000
