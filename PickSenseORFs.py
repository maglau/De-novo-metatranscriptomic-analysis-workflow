
#!/usr/bin python

# PickSenseORFs.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Nov 04, 2014

# Protein-coding genes were called using Prodigal which search for all 6 translation frames. Since the transcript reads and contigs were known to be sense strands, this script is to remove those reverse strand genes based on the header info in the prodigal.faa file. The removed sequences will be written into a new file as 'prodigal_removed.faa', while those from predicted from the sense strand are retained and written into a 'prodigal_sense.faa'. 

# define usage message
Usage =  """Usage: PickSenseORFs.py [Trinity.fasta.prodigal.faa] [outfile prefix]"""

import sys
from Bio import SeqIO

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) != 3:
    print Usage
    exit(1)

prodigalORF = open(sys.argv[1], 'rU') # open the file as a text file
fileprefix = sys.argv[2]
sense = open(fileprefix+"_sense.faa", 'w') # change to "_sense.fna" for nucleotide sequences
remove = open(fileprefix+"_removed.faa", 'w') # change to "_removed.fna" for nucleotide sequences

for each in SeqIO.parse(prodigalORF, "fasta"):
	col = each.description.split('#')
#	print '%d' % int(col[3])
	if int(col[3]) > 0:
		SeqIO.write(each, sense, 'fasta')
	else:
		SeqIO.write(each, remove, 'fasta')

prodigalORF.close()
sense.close()
remove.close()
