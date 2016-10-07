#!/usr/bin python

# merge_BestRep_clstrCOL_v2_singlesample.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Apr 18, 2016

# add the consensus protein results in "_newBestRep.txt" to the ".clstrCOL"

# define usage message
Usage =  """Usage: python merge_BestRep_clstrCOL_v2_singlesample.py [_newBestRep.txt] [.clsterCOL]"""

import sys

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) != 3:
    print Usage
    exit(1)

BestRep = open(sys.argv[1], 'rU') # open the file as a text file
clstr = open(sys.argv[2], 'rU')
merged = open(sys.argv[2]+"_newBestRep", 'w')

BestRepdict = {}

for line in BestRep:
	if line.startswith('#') == False:
		line = line.strip()
		col = line.split("\t")
		title = col[0].split("_")
		seqID = '_'.join(title[1:])
		info = '\t'.join(col[1:])
		BestRepdict[seqID] = info	


# some representative sequences yielded no blast results (! e.g. Cluster6_Trinity|comp8954_c0_seeq3_1) using the stand-a-lone blast engine, and online BLAST site tells that "no putative conserved domains have been detected)


for line in clstr:
	line = line.strip()
	col = line.split("\t")
	SeqID = col[1]
	if SeqID in BestRepdict:
		merged.write('%s\t%s\n' % (line, BestRepdict[SeqID]))
	
	else:
		merged.write('%s\t%s\t%d\t%d\n' % (line, 'blast detected no putative conserved domains', 0, 0))

BestRep.close()
clstr.close()
merged.close()







