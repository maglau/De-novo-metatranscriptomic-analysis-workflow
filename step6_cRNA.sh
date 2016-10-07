#!/bin/bash

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git

# touch step6_cRNA.log

clear
echo "---Metatranscriptomic workflow step 6 script v.1.0---" > step6_cRNA.log

echo "" >> step6_cRNA.log
echo "Find consensus protein ID for each of the query sequences" >> step6_cRNA.log
echo "python blast-majority_v3.py cd-hit/reads_cRNA_sense90G_NRbest10.out cd-hit/reads_cRNA_sense90G.table" >> step6_cRNA.sh
python blast-majority_v3.py cd-hit/reads_cRNA_sense90G_NRbest10.out cd-hit/reads_cRNA_sense90G.table

echo "" >> step6_cRNA.log
echo "Merge the blast results to the list of cluster and rep sequences" >> step6_cRNA.log
echo "python merge_BestRep_clstrCOL_v2_singlesample.py cd-hit/reads_cRNA_sense90G_NRbest10_newBestRep.txt cd-hit/reads_cRNA_sense90G.clstrCOL" >> step6_cRNA.log
python merge_BestRep_clstrCOL_v2_singlesample.py cd-hit/reads_cRNA_sense90G_NRbest10_newBestRep.txt cd-hit/reads_cRNA_sense90G.clstrCOL

echo "" >> step6_cRNA.log
echo "Add the read counts and coverage to the above output files" >> step6_cRNA.log
echo "python add_cRNAs_ReadCounts_Coverage_v3_singlesample.py cRNA/reads_trinity__coverage.txt cd-hit/reads_cRNA_sense90G.clstrCOL_newBestRep cd-hit/reads_cRNA_sense90G.table.updated_wConsensusLength" >> step6_cRNA.log
python add_cRNAs_ReadCounts_Coverage_v3_singlesample.py cRNA/reads_trinity_coverage.txt cd-hit/reads_cRNA_sense90G.clstrCOL_newBestRep cd-hit/reads_cRNA_sense90G.table.updated_wConsensusLength

python blast-majority_v3.py cd-hit/reads_cRNA_sense90G_NRbest10.out cd-hit/reads_cRNA_sense90G.table
python blast-majority_v3.py cd-hit/reads_cRNA_sense90G_NRbest10.out cd-hit/reads_cRNA_sense90G.table
python blast-majority_v3.py cd-hit/reads_cRNA_sense90G_NRbest10.out cd-hit/reads_cRNA_sense90G.table
