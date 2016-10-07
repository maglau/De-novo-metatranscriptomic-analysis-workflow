#!/bin/bash

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git

touch step6_ncRNA.log

clear
echo "---Metatranscriptomic workflow step 5 script v.1.0---" > step6_ncRNA.log

# -------
# Based on how many lines the input file contained, determine how many subfiles will be created
# If the default number as given in step4_ncRNA is used, each subfile will contain 2000 lines
# Each subfile will have a line similar to the one below, with 'X' replaced by '1....n' for n files

# EDIT THIS LINE:
echo "" >> step6_ncRNA.log
echo "blastn -task megablast -db SILVA_119_rRNA -query cd-hit/reads_ncRNA_95G_X.fasta -max_target_seqs 10 -num_threads 16 -outfmt 11 -out cd-hit/reads_ncRNA_95G_SILVAbest10_X.asn" >> step6_ncRNA.log
blastn -task megablast -db SILVA_119_rRNA -query cd-hit/reads_ncRNA_95G_X.fasta -max_target_seqs 10 -num_threads 16 -outfmt 11 -out cd-hit/reads_ncRNA_95G_SILVAbest10_X.asn

# Next, combine all subfiles as follows
# EDIT THIS LINE:
echo "" >> step6_ncRNA.log
echo "cat cd-hit/reads_ncRNA_95G_SILVAbest10_1.asn ... cd-hit/reads_ncRNA_95G_SILVAbest10_n.asn > cd-hit/reads_ncRNA_95G_SILVAbest10.asn" >> step6_ncRNA.log
cat cd-hit/reads_ncRNA_95G_SILVAbest10_1.asn ... cd-hit/reads_ncRNA_95G_SILVAbest10_n.asn > cd-hit/reads_ncRNA_95G_SILVAbest10.asn

echo "" >> step6_ncRNA.log
echo "blast_formatter -archive cd-hit/reads_ncRNA_95G_SILVAbest10.asn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out cd-hit/reads_ncRNA_95G_SILVAbest10.out" >> step6_ncRNA.log
blast_formatter -archive cd-hit/reads_ncRNA_95G_SILVAbest10.asn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out cd-hit/reads_ncRNA_95G_SILVAbest10.out


# EDIT THIS LINE:
echo "" >> step6_ncRNA.log
echo "rm cd-hit/reads_ncRNA_95G_SILVAbest10_1.asn .... cd-hit/reads_ncRNA_95G_SILVAbest10_n.asn cd-hit/reads_ncRNA_95G_1.fasta cd-hit/reads_ncRNA_95G_n.fasta" >> step6_ncRNA.log
rm cd-hit/reads_ncRNA_95G_SILVAbest10_1.asn .... cd-hit/reads_ncRNA_95G_SILVAbest10_n.asn cd-hit/reads_ncRNA_95G_1.fasta cd-hit/reads_ncRNA_95G_n.fasta
