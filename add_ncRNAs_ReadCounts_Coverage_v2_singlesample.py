#!/usr/bin python

# add_ncRNAs_ReadCounts_Coverage_v2_singlesample.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Apr 18, 2016

# this script adds the read counts (RC) and coverage in "__coverage.txt" to the corresponding transcript contig.


# define usage message
Usage =  """Usage: python add_ncRNAs_ReadCounts_Coverage_v2_singlsample.py [transabyss__coverage.txt ** BT2 version] [xxx.clstrCOL_LCAdiversity] [.table.updated_wConsensusLength]"""


import sys
import re

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) <= 2:
    print Usage
    exit(1)

transabyss = open(sys.argv[1], 'rU')

### extract the Read Counts (RC) and coverage from the "__coverage.txt" (BT2 version) 

def extractRC(x):
	RCdict = {}
	for line in x:
		line = line.strip()
		col = line.split('\t')
		transcriptID = col[0]
		RC = col[2]
		cov = col[3]
		RCdict[transcriptID] = (RC, cov)
	return RCdict

transabyssRCdict = extractRC(transabyss)

### add the Read Counts (RC) and coverage to the ".clstrCOL_LCAdiversity"

clstrCOL = open(sys.argv[2], 'rU')
newclstrCOL = open(sys.argv[2]+"_ReadCountsCoverage", 'w')

for line in clstrCOL:
	line = line.strip()
	col = line.split('\t')
	TranscriptID = col[2]
	if TranscriptID in transabyssRCdict:
		newclstrCOL.write('%s\t%s\t%s\n' % ('\t'.join(col[0:4]), '\t'.join(transabyssRCdict[TranscriptID]), '\t'.join(col[4:])))
	else:
		newclstrCOL.write('%s\t%d\t%d\t%s\n' % ('\t'.join(col[0:4]), 0, 0, '\t'.join(col[4:])))

### in the "_ReadCounts" file, the columns are [0] Cluster [1] RepSeq [2] Assembler [3] transcript contig [4] transcript contig length [5] read count [6] coverage [7] aln to the RepSeq [8] strand [9] % identity to the RepSeq [10] consensus rRNA molecule [11] max length of the consensus rRNA [12] taxonomy



### add the Read Counts (RC) and coverage to the ".table.updated_wConsensusLength"
newclstrCOL = open(sys.argv[2]+"_ReadCountsCoverage", 'rU')

RCdict = {}
RCsum = 0
covsum =0

for line in newclstrCOL:
	line = line.strip()
	col = line.split('\t')
	RefSeq = col[1]
	RC = int(col[4])
	cov = float(col[5])
        if RefSeq not in RCdict:
                RCdict[RefSeq] = (RC, cov)
        else:
                RCsum = RCdict[RefSeq][0] + RC
                covsum = RCdict[RefSeq][1] + cov
                RCdict[RefSeq] = (RCsum, covsum)
                RCsum = 0
                covsum = 0


table = open(sys.argv[3], 'rU')
newtable = open(sys.argv[3]+"_ReadCountsCoverage", 'w')

newtable.write('# cluster\tcluster size\trepresentative rRNA seq\ttranscript contig length\tread count\tcoverage\tconsensus rRNA molecule\tmax length of the consensus rRNA\ttaxonomy\n')

for line in table:
	if line.startswith('#') == False:
		line = line.strip()
		col = line.split('\t')
		RefSeqID = col[2]
		if RefSeqID not in RCdict:
                        RCdict[RefSeqID] = (0, 0)
		newtable.write('%s\t%s\t%s\n' % ('\t'.join(col[0:4]), '\t'.join(tuple(str(tr) for tr in RCdict[RefSeqID])), '\t'.join(col[4:])))

transabyss.close()
newclstrCOL.close()
table.close()
newtable.close()


