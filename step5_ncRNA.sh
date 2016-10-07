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

# touch step5_ncRNA.log
clear
echo "---Metagenomics workflow step 4 script v.1.0---" > step5_ncRNA.log

echo "" >> step5_ncRNA.log
echo "UCHIME: Compares ncRNA contigs with non-coding sequences in databases" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/LSUParc119_ATCG.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeLSU -chimeras ncRNA/reads_ncRNA_transabyss.chimerasLSU" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/LSUParc119_ATCG.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeLSU -chimeras ncRNA/reads_ncRNA_transabyss.chimerasLSU

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_1.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU1 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU1" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_1.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU1 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU1

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_2.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU2 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU2" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_2.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU2 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU2

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_3.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU3 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU3" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_3.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU3 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU3

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_4.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU4 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU4" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_4.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU4 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU4

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_5.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU5 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU5" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_5.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU5 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU5

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_6.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU6 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU6" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_6.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU6 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU6

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_7.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU7 -chimerasn ncRNA/reads_ncRNA_transabyss.chimerasSSU7" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_7.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU7 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU7

echo "" >> step5_ncRNA.log
echo "usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_8.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU8 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU8" >> step5_ncRNA.log
usearch -uchime_ref ncRNA/transabyss-merged-reads_nc_new.fa -db /scratch/gpfs/maglau/DB/SSUParc119_ATCG_8.udb -strand plus -threads 8 -uchimeout ncRNA/reads_ncRNA_transabyss.uchimeSSU8 -chimeras ncRNA/reads_ncRNA_transabyss.chimerasSSU8

echo "" >> step5_ncRNA.log
echo "Concatenates all chimera files" >> step5_ncRNA.log
echo "cat reads_ncRNA_transabyss.chimerasLSU reads_ncRNA_transabyss.chimerasSSU1 reads_ncRNA_transabyss.chimerasSSU2 reads_ncRNA_transabyss.chimerasSSU3 reads_ncRNA_transabyss.chimerasSSU4 reads_ncRNA_transabyss.chimerasSSU5 reads_ncRNA_transabyss.chimerasSSU6 reads_ncRNA_transabyss.chimerasSSU7 reads_ncRNA_transabyss.chimerasSSU8 | grep ">" > reads_ncRNA_transabyss.chimerastitles" >> step5_ncRNA.log
cat ncRNA/reads_ncRNA_transabyss.chimerasLSU ncRNA/reads_ncRNA_transabyss.chimerasSSU1 ncRNA/reads_ncRNA_transabyss.chimerasSSU2 ncRNA/reads_ncRNA_transabyss.chimerasSSU3 ncRNA/reads_ncRNA_transabyss.chimerasSSU4 ncRNA/reads_ncRNA_transabyss.chimerasSSU5 ncRNA/reads_ncRNA_transabyss.chimerasSSU6 ncRNA/reads_ncRNA_transabyss.chimerasSSU7 ncRNA/reads_ncRNA_transabyss.chimerasSSU8 | grep ">" > ncRNA/reads_ncRNA_transabyss.chimerastitles
rm " > reads_ncRNA_transabyss.chimerastitles"

echo "" >> step5_ncRNA.log
echo "Removes chimera sequences from CD-HIT output files" >> step5_ncRNA.log
echo "python reformat_cdhit_clstr_v3_ncRNA_singlesample.py cd-hit/reads_ncRNA_95G.clstr cd-hit/reads_ncRNA_95G ncRNA/reads_ncRNA_transabyss.chimerastitles" >> step5_ncRNA.log
python reformat_cdhit_clstr_v3_ncRNA_singlesample.py cd-hit/reads_ncRNA_95G.clstr cd-hit/reads_ncRNA_95G ncRNA/reads_ncRNA_transabyss.chimerastitles

echo "" >> step5_ncRNA.log
echo "blastn -task megablast -db SILVA_119_rRNA -query cd-hit/reads_ncRNA_95G.renameTrimmed -max_target_seqs 10 -num_threads 16 -outfmt 11 -out cd-hit/reads_ncRNA_95G_SILVAbest10.asn" >> step5_ncRNA.log
blastn -task megablast -db SILVA_119_rRNA -query cd-hit/reads_ncRNA_95G.renameTrimmed -max_target_seqs 10 -num_threads 16 -outfmt 11 -out cd-hit/reads_ncRNA_95G_SILVAbest10.asn

echo "" >> step5_ncRNA.log
echo "blast_formatter -archive cd-hit/reads_ncRNA_95G_SILVAbest10.asn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out cd-hit/reads_ncRNA_95G_SILVAbest10.out" >> step5_ncRNA.log
blast_formatter -archive cd-hit/reads_ncRNA_95G_SILVAbest10.asn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out cd-hit/reads_ncRNA_95G_SILVAbest10.out

echo "" >> step5_ncRNA.log
echo "python Split_FASTA.py cd-hit/reads_ncRNA_95G.renameTrimmed 2000" >> step5_ncRNA.log
python Split_FASTA.py cd-hit/reads_ncRNA_95G.renameTrimmed 2000
