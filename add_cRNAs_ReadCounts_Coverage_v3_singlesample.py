#!/usr/bin python

# add_cRNAs_ReadCounts_Coverage_v2_singlesample.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Apr 18, 2016

# this script adds the read counts (RC) and coverage in "__coverage.txt" to the corresponding parent contig of each PEG.



# this version is cleaner and computationally less demanding (i.e. faster)

# define usage message
Usage =  """Usage: python add_cRNAs_ReadCounts_Coverage.py [trinity_coverage.txt ** BT2 version] [xxx.clstrCOL_newBestRep] [.table.updated_wConsensusLength]"""


import sys
import re

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) <= 2:
    print Usage
    exit(1)

trinity = open(sys.argv[1], 'rU') # open the file as a text file

### extract the Read Counts (RC) and coverage from the "_coverage.txt" (BT2 version) 

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

trinityRCdict = extractRC(trinity)

### add the Read Counts (RC) and coverage to the ".clstrCOL_newBestRep"

clstrCOL = open(sys.argv[2], 'rU')
newclstrCOL = open(sys.argv[2]+"_ReadCountsCoverage", 'w')

for line in clstrCOL:
	line = line.strip()
	col = line.split('\t')
	TranscriptID = col[4]
	if TranscriptID in trinityRCdict:
		newclstrCOL.write('%s\t%s\t%s\n' % ('\t'.join(col[0:5]), '\t'.join(trinityRCdict[TranscriptID]), '\t'.join(col[5:])))
	else:
		newclstrCOL.write('%s\t%d\t%d\t%s\n' % ('\t'.join(col[0:5]), 0, 0, '\t'.join(col[5:])))

### in the "_ReadCounts" file, the columns are [0] Cluster [1] RepSeq [2] Assembler [3] PEG [4] PEGlength [5] Parent transcript contig [6] Parent read count [7] Parent coverage [8] aln to the RepSeq [9] % identity to the RepSeq [10] consensus protein ID [11] consensus protein length [12] number of hit (out of 10) passed quality thresholds



### add the Read Counts (RC) to the ".table.updated_wConsensusLength"
newclstrCOL = open(sys.argv[2]+"_ReadCountsCoverage", 'rU')

RCdict = {}
RCsum = 0
covsum = 0

for line in newclstrCOL:
	line = line.strip()
	col = line.split('\t')
	RefSeq = col[1]
	RC = int(col[5])
	cov = float(col[6])
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

newtable.write('# cluster\tcluster size\trepresentative PEG seq\tPEG length\tread count\tcoverage\tconsensus protein ID\tmax length of the consensus protein\tnumber of qualified matches\n')

for line in table:
	if line.startswith('#') == False:
		line = line.strip()
		col = line.split('\t')
		RefSeqID = col[2]
		if RefSeqID not in RCdict:
			RCdict[RefSeqID] = (0, 0)
		newtable.write('%s\t%s\t%s\n' % ('\t'.join(col[0:4]), '\t'.join(tuple(str(t) for t in RCdict[RefSeqID])), '\t'.join(col[4:])))


trinity.close()
newclstrCOL.close()
table.close()
newtable.close()


