#!/usr/bin python

# merge_LCAdiversity_clstrCOL_v2_singlesample.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Apr 18, 2016

# add the consensus rRNA results in "_LCAdiversity.txt" to the ".clstrCOL"

# define usage message
Usage =  """Usage: merge_BestRep_clstrCOL_v2_singlesample.py [_LCAdiversity.txt] [.clsterCOL]"""

import sys

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) != 3:
    print Usage
    exit(1)

BestRep = open(sys.argv[1], 'rU') # open the file as a text file
clstr = open(sys.argv[2], 'rU')
merged = open(sys.argv[2]+"_LCAdiversity", 'w')

LCAdict = {}

for line in BestRep:
	if line.startswith('#') == False:
		line = line.strip()
		col = line.split("\t")
		title = col[0].split("_")
		seqID = '_'.join(title[1:])
		info = '\t'.join(col[1:])
		LCAdict[seqID] = info	


# some representative sequences yielded no blast results


for line in clstr:
	line = line.strip()
	col = line.split("\t")
	SeqID = col[1]
	if SeqID in LCAdict:
		merged.write('%s\t%s\n' % (line, LCAdict[SeqID]))
	
	else:
		merged.write('%s\t%s\n' % (line, 'blast detected no hits'))

BestRep.close()
clstr.close()
merged.close()







