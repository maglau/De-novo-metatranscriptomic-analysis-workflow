#!/bin/bash

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git

# touch step7_ncRNA.log

clear
echo "---Metatranscriptomic workflow step 6 script v.1.0---" > step7_ncRNA.log

echo "" >> step7_ncRNA.log
echo "Find taxonomic classification for each of the query sequences" >> step7_ncRNA.log
echo "python blast-basedLCA80.py cd-hit/reads_ncRNA_95G_SILVAbest10.out cd-hit/reads_ncRNA_95G.tableTrimmed" >> step7_ncRNA.sh
python blast-basedLCA80.py cd-hit/reads_ncRNA_95G_SILVAbest10.out cd-hit/reads_ncRNA_95G.tableTrimmed

echo "" >> step7_ncRNA.log
echo "Merge the blast results to the list of cluster and rep sequences" >> step7_ncRNA.log
echo "python merge_LCAdiversity_clstrCOL_v2_singlesample.py cd-hit/reads_ncRNA_95G_SILVAbest10_LCAdiversity.txt cd-hit/reads_ncRNA_95G.clstrCOLtrimmed" >> step7_cRNA.log
python merge_LCAdiversity_clstrCOL_v2_singlesample.py cd-hit/reads_ncRNA_95G_SILVAbest10_LCAdiversity.txt cd-hit/reads_ncRNA_95G.clstrCOLtrimmed

echo "" >> step7_ncRNA.log
echo "Add the read counts and coverage to the above output files" >> step7_ncRNA.log
echo "python add_ncRNAs_ReadCounts_Coverage_v2_singlesample.py ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_BT2vs__coverage.txt cd-hit/reads_ncRNA_95G.clstrCOLtrimmed_LCAdiversity cd-hit/reads_ncRNA_95G.tableTrimmed.updated_wConsensusLength" >> step7_ncRNA.log
python add_ncRNAs_ReadCounts_Coverage_v2_singlesample.py ncRNA/reads_ncRNAreadsXncRNAcontigs_transabyss_BT2vs__coverage.txt cd-hit/reads_ncRNA_95G.clstrCOLtrimmed_LCAdiversity cd-hit/reads_ncRNA_95G.tableTrimmed.updated_wConsensusLength
python blast-majority_v3.py cd-hit/reads_ncRNA_95G_SILVAbest10.out cd-hit/reads_ncRNA_95G.tableTrimmed
python blast-basedLCA80.py cd-hit/reads_ncRNA_95G_SILVAbest10.out cd-hit/reads_ncRNA_95G.tableTrimmed
python blast-basedLCA80.py cd-hit/reads_ncRNA_95G_SILVAbest10.out cd-hit/reads_ncRNA_95G.tableTrimmed
python blast-basedLCA80.py cd-hit/reads_ncRNA_95G_SILVAbest10.out cd-hit/reads_ncRNA_95G.tableTrimmed
