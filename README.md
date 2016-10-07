# De-novo-metatranscriptomic-analysis-workflow
RNA-sequencing is an important tool to reveal the identity and metabolic capability of the active members. Total RNA (i.e. ribosomal and messenger RNAs) can be used to profile the taxonomic and functional diversity in an environmental sample. This workflow first separates rRNA from mRNA and assembles them independently. Each assembly is followed by annotation and statistics calculation.

------

Purpose: This set of shell scripts functions as a metatranscriptomics study of coding (cRNA) and non-coding RNA (ncRNA) sequence reads. They utilize the Trinity and Trans-ABySS assemblers to produce contig scaffolds for the cRNA and ncRNA reads, respectively (M. Grabber et al. 2011, G. Robertson et al. 2010). Protein-encoding genes (PEGs) and taxonomic classifications are outputted with corresponding coverage values by matching resulting contigs to the SILVA RNA database (http://www.arb-silva.de/).

------

Required Algorithms:

- USEARCH (http://www.drive5.com/usearch/download.html)
- BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- Prodigal (http://prodigal.ornl.gov/downloads.php)
- Trinity (https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- Trans-ABySS (https://github.com/bcgsc/transabyss)
- ABySS (http://www.bcgsc.ca/platform/bioinfo/software/abyss)
- cd-hit (http://weizhongli-lab.org/cd-hit/download.php)

Required Databases:
- SILVA (https://www.arb-silva.de/download/archive/)
- BLAST 

------

Make the following directories:

- cRNA
- ncRNA
- cd-hit

------ 

Required files at the first directory level (in addition to slurm and shell scripts):
- [input xxx.fasta file of all sequence reads]
- Calculate_avReadLen.py
- Calculate_stats_v3.py
- ParseSeqsIntoGroups_v4.py
- reformat_cdhit_clstr_v2_cRNA_singlesample.py
- reformat_cdhit_clstr_v3_ncRNA_singlesample.py
- Split_FASTA.py
- UsearchStat_v3.py
- merge_BestRep_clstrCOL_v2_singlesample.py
- merge_LCAdiversity_clstrCOL_v2_singlesample.py
- add_cRNAs_ReadCounts_Coverage_v3_singlesample.py
- add_ncRNAs_ReadCounts_Coverage_v2_singlesample.py
- blast-basedLCA80.py
- blast-majority_v3.py

------

Shell Scripts (Shell scripts can be run independently of slurm scripts if sample is small and a supercomputer is not needed. *****MUST BE EDITED*****: Paths to all algorithms and databases within shell scripts will need editing):

- step1.sh: Compares the user-provided RNA sequence reads against four SILVA databses (small subunit rRNA and large subunit rRNA), and transfer RNA (tRNA) and 5S ribosomal RNA (rRNA), through USEARCH (R.C. Edgar. 2010), then separates cRNA and ncRNA reads and dumps the contents into two .fasta files. *****MUST BE EDITED*****: database names will need editing according to downloaded version. Add xxx.fasta file containing all sequence reads. If the user desires, each line can be separated into individual slurm scripts made to run in parallel in the interest of speeding up the process.


~ For processing cRNA reads:

- step2_cRNA.sh: Runs the Trinity assembler on the cRNA reads to generate contig scaffolds. Should be edited if a different assmebler is used. 

- step3_cRNA.sh: Maps the cRNA dataset created in step1.sh to the cRNA assembly created in step2_cRNA.sh using Bowtie 2 version 2.2.5 (B. Langmead & S.L. Salzberg. 2012) and calculates contig coverage and statistics. 

- step4_cRNA.sh: Uses Prodigal version 2.6.1 (D. Hyatt et al. 2010) to generate PEGs from cRNA transcript contigs, applies the CD-HIT algorithm version 4.6.4 (L. Fu et al. 2012) to cluster the resulting PEGs, and splits resulting PEG cluster file into subfiles of 2000 lines (can be edited) to create files of a more manageable size for BLAST search.

- step5_cRNA.sh: Searches the representative (i.e. longest) sequence in each PEG cluster against the NCBI non-redundant protein (nr) database (http://www.ncbi.nlm.nih.gov/refseq/) for the best ten hits using BLASTp (C. Camacho et al. 2009) to obtain consensus protein identities. *****MUST BE EDITED*****: User must add lines to process each subfile; the number of lines will depend on the size of the user’s input file as each subfile is automatically set to be 2000 lines long.

- step6_cRNA.sh: Constructs a table of protein IDs corresponding to each PEG query, along with read counts and coverage.


~ For processing ncRNA reads:

- step2_ncRNA.sh: Runs the Trans-ABySS assembler on the ncRNA reads to generate contig scaffolds.

- step3_ncRNA.sh: Maps the ncRNA dataset created in step1.sh to the ncRNA assembly also created in step2_ncRNA.sh using Bowtie 2 version 2.2.5 (B. Langmead & S.L. Salzberg. 2012).

- step4_ncRNA.sh: Calculates the ncRNA transcript contig coverage and statistics, and applies the CD-HIT-EST (L. Weizhong & A. Godzik. 2006) algorithm to cluster the contigs.

- step5_ncRNA.sh: Identifies chimeric gene sequences using UCHIME (R.C. Edgar et al. 2011) and splits resulting contig cluster file into subfiles of 2000 lines each for easier BLASTing. The ncRNA transcript contigs in clusters that were represented by a chimeric sequence are omitted from further analysis.

- step6_ncRNA.sh: Searches the representative sequence in each contig cluster against the rRNA database comprised of SILVA SSU and LSU sequences for the best ten hits using BLASTn (C. Camacho et al. 2009) to obtain consensus rRNA identities by majority rule, and taxonomic rankings by the lowest common ancestor (LCA) principle.*****MUST BE EDITED*****: User must add lines to process each subfile; the number of lines will depend on the size of the user’s input file as each subfile is automatically set to be 2000 lines long.

- step7_ncRNA.sh: Constructs a table of taxonomic classifications corresponding to each contig cluster sequence query, along with read counts and coverage.

------

Slurm Scripts (*****MUST BE EDITED*****: [user email] and [path to working directory] in all slurm script headings must be replaced with desired email and path: 

- SlurmScript_step1: *****MUST BE EDITED*****: user must add the desired input file name.

- Each script can be edited to specify the desired number of CPU cores and allotted run time.

------

Bibliography:

- Grabherr, M. G. et al. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat. Biotechnol. 29, 644–652 (2011).

- Robertson, G. et al. De novo assembly and analysis of RNA-seq data. Nat. Methods 7, 909–912 (2010).

- Edgar, R. C. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26, 2460–2461 (2010).

- Langmead, B. & Salzberg, S. L. Fast gapped-read alignment with Bowtie 2. Nat. Methods 9, 357–359 (2012).

- Hyatt, D. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010).

- Fu, L., Niu, B., Zhu, Z., Wu, S. & Li, W. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics 28, 3150–3152 (2012).

- Camacho, C. et al. BLAST+: architecture and applications. BMC Bioinformatics 10, 421 (2009).

- Weizhong, L & Godzik, A.  Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics 22, 1658-1659 (2006).

- Edgar, R. C., Haas, B. J., Clemente, J. C., Quince, C. & Knight, R. UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27, 2194–2200 (2011).
