#!/usr/bin python

# blast-basedLCA3.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Oct 10, 2015

# ncRNA contigs searched against the SILVA SSU NR99+LSU database with the best 10 hits recorded. based on the taxonomic lineage information given in the description, a consensus lineage is consolidated based on 80% of confidence and down to the lowest common taxonomic rank possible.

# define usage message
Usage =  """Usage: blast-basedLCA3.py [blast result table] [.table]"""

import sys

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) != 3:
    print Usage
    exit(1)

SILVA = open(sys.argv[1], 'rU') # open the file as a text file
fileprefix = sys.argv[1][:-4]
LCA = open (fileprefix+"_LCAdiversity.txt", 'w')

prevID = None
prevMol = None
hitnum = 1
hitlenlist = []

LCA.write("# qseqid\tconsensus rRNA\tmax length of the consensus rRNA\ttaxonomy\n")


for line in SILVA:
        line = line.strip()
	col = line.split("\t")
	seqID = col[0]

	if prevID == None or prevID != seqID:
		if prevID != None:
			if hitlenlist != []:
				LCA.write('%s\t%s\t%d\t' % (prevID, prevMol, max(hitlenlist)))
			else:
				LCA.write('%s\t%s\t%d\t' % (prevID, 'unknown', 0))
			consensus = [] 
			confidences = []	
			for rankDict in tax:
				goodhitsum = sum(rankDict.values())
				for rank in rankDict:
                                        confidence = rankDict[rank]/float(goodhitsum)
                                        confidences.append(confidence)
					if confidence >= 0.8:
						LCA.write( '%s\t' % (rank) )
				if all(i < 0.8 for i in confidences):			
					LCA.write('unknown\t')
				confidences = []
			LCA.write ('\n')
			hitlenlist = []

		
		tax = [dict() for x in range(10)]

        cov = (float(col[8])/float(col[1]))*100
        e = float(col[12])
        bit = float(col[13])
        if int(cov) >= 50  and e <= 1e-05 and bit >=50:
		ranks = col[-1].split(";")
		numshort= 10 - len(ranks)
		ranks.extend(['unknown']*numshort)
		for rank, rankDict in zip(ranks, tax):
			ct = rankDict.setdefault(rank, 0)
			rankDict[rank] = ct + 1	
		hitlen = int(col[3])
		hitlenlist.append(hitlen)

	prevID = seqID
	prevMol = col[2][0:3]

if hitlenlist != []:
	LCA.write('%s\t%s\t%d\t' % (prevID, prevMol, max(hitlenlist)))
else:
	LCA.write('%s\t%s\t%d\t' % (prevID, 'unknown', 0))
consensus = []
confidences = []
for rankDict in tax:
	goodhitsum = sum(rankDict.values())
        for rank in rankDict:
		confidence = rankDict[rank]/float(goodhitsum)
                confidences.append(confidence)
                if confidence >= 0.8:
                       	LCA.write( '%s\t' % (rank) )
	if all(i < 0.8 for i in confidences):
        	LCA.write('unknown\t')
        confidences = []
LCA.write ('\n')


SILVA.close()
LCA.close()


#######

# this part of script adds the consensus rRNA identity to the summary file - xxx.table.

fileprefix = sys.argv[1][:-4]
infile = open(fileprefix+"_LCAdiversity.txt", 'rU')
oldsummary = open(sys.argv[2], 'rU')
newsummary = open(sys.argv[2]+".updated_wConsensusLength", 'w')

clstrIDdict = {}

for line in infile:
        line = line.strip()
        col = line.split("\t")
        key = col[0].split("_")
        clstr = key[0]
        clstrIDdict[clstr] = '\t'.join(col[1:])

# some representative sequences yielded no blast results

newsummary.write("# cluster\tcluster size\trepresentative seq\tseq length\tconsensus rRNA molecule\tmax length of the consensus protein\ttaxonomy\n")

for line in oldsummary:
        if line.startswith("#") == False:
                line = line.strip()
                col = line.split("\t")
                key = col[0]
                if key in clstrIDdict:
                        newsummary.write("%s\t%s\n" % (line, clstrIDdict[key]))
                else:
                        newsummary.write("%s\tblast detected no hits\n" % (line))

infile.close()
oldsummary.close()
newsummary.close()




