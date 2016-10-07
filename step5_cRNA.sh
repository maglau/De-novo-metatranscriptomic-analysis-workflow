#!/bin/bash

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git

touch step5_cRNA.log

clear
echo "---Metatranscriptomic workflow step 5 script v.1.0---" > step5_cRNA.log

# Runs BLASTp on representative sequences of PEG clusters to obtain consensus protein IDs

# -------
# EDIT THIS LINE
# Based on how many lines the input file contained, determine how many subfiles will be created 
# If the default number as given in step4_cRNA is used, each subfile will contain 2000 lines
# Each subfile will have a line similar to the one below, with 'X' replaced by '1....n' for n files  
echo "blastp -db nr -query cd-hit/reads_cRNA_sense90G_X.fasta -max_target_seqs 10 -num_threads 16 -outfmt 11 -out cd-hit/reads_cRNA_sense90G_NRbest10_X.asn" >> step5_cRNA.log
blastp -db nr -query cd-hit/reads_cRNA_sense90G.rename -max_target_seqs 10 -num_threads 16 -outfmt 11 -out cd-hit/reads_cRNA_sense90G_NRbest10_X.asn

# EDIT THIS LINE:
echo "" >> step5_cRNA.log
echo "Combines all subfiles" >> step5_cRNA.log
echo "cat cd-hit/reads_cRNA_sense90G_NRbest10_1.asn ... cd-hit/reads_cRNA_sense90G_NRbest10_n.asn > cd-hit/reads_cRNA_sense90G_NRbest10.asn" >> step5_cRNA.log
cat cd-hit/reads_cRNA_sense90G_NRbest10_1.asn ... cd-hit/reads_cRNA_sense90G_NRbest10_n.asn > cd-hit/reads_cRNA_sense90G_NRbest10.asn

echo "" >> step5_cRNA.log
echo "blast_formatter -archive cd-hit/reads_cRNA_sense90G_NRbest10.asn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out cd-hit/reads_cRNA_sense90G_NRbest10.out" >> step5_cRNA.log
blast_formatter -archive cd-hit/reads_cRNA_sense90G_NRbest10.asn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out cd-hit/reads_cRNA_sense90G_NRbest10.out


# EDIT THIS LINE:
echo "" >> step5_cRNA.log
echo "Removes all subfiles" >> step5_cRNA.log
echo "rm cd-hit/reads_cRNA_sense90G_NRbest10_1.asn .... cd-hit/reads_cRNA_sense90G_NRbest10_n.asn" >> step5_cRNA.log
rm cd-hit/reads_cRNA_sense90G_NRbest10_1.asn .... cd-hit/reads_cRNA_sense90G_NRbest10_n.asn cd-hit/reads_cRNA_sense90G_1.fasta .... cd-hit/reads_cRNA_sense90G_n.fasta


