#!/usr/bin python

# blast-majority_v3.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on Oct 06, 2015

# with the blastp output of the top 10 bits, this script will determine, on the basis of majority, the best hit among those that have met the criteria: alignment covers > 50% of the query sequence, e-value < E-05and bitscore > 50. The column titles in the blastp output table: qseqid [0] qlen sseqid slen qstart qend [5] sstart send length pident nident [10] mismatch evalue bitscore staxids saccver [15] stitle 

# this version is cleaner than v2, and it also outputs the length of the longest consensus protein

# define usage message
Usage =  """Usage: python blast-majority_v3.py [blastp result table] [xxx.table]"""

import sys
import re
from collections import Counter

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) <= 2:
    print Usage
    exit(1)

blastp = open(sys.argv[1], 'rU') # open the file as a text file
fileprefix = sys.argv[1][:-4]
outfile = open(fileprefix+"_newBestRep.txt", 'w')

hitlist = []
hitlen = 0
prevpeptide = None
consensusID = None

outfile.write("# qseqid\tconsensus protein ID\tmax length of the consensus protein\tnumber of matches\n")

for line in blastp:
        line = line.strip()
	col = line.split("\t")
        cov = (float(col[8])/float(col[1]))*100
        e = float(col[12])
        bit = float(col[13])
	peptide = col[0]
	stitle = col[16].split(" [")
	match = stitle[0].replace("MULTISPECIES: ", "") 
	hitlen = int(col[3])
	if prevpeptide != None and peptide != prevpeptide:
		if hitlist == []:
                        outfile.write("%s\t%s\t%d\t%d\n" % (prevpeptide, None, 0, 0))
		else:
			proteinCnt = Counter([x[0] for x in hitlist]) # set up a counter for the protein matches that stored in the hitlist
			consensusID,hitCnt = proteinCnt.most_common(1)[0]
			consensuslen = dict(hitlist)[consensusID]
			outfile.write("%s\t%s\t%d\t%d\n" % (prevpeptide, consensusID, consensuslen, hitCnt))
		hitlist = []
		prevpeptide = peptide
			
	if int(cov) >= 50 and e <= 1e-05 and bit >=50:
		hitlist.append((match,hitlen))
		prevpeptide = peptide
	
if hitlist == []:
	outfile.write("%s\t%s\t%d\t%d\n" % (prevpeptide, None, 0, 0))
else:
	proteinCnt = Counter([x[0] for x in hitlist]) # set up a counter for the protein matches that stored in the hitlist
	consensusID,hitCnt = proteinCnt.most_common(1)[0]
	consensuslen = dict(hitlist)[consensusID]
	outfile.write("%s\t%s\t%d\t%d\n" % (prevpeptide, consensusID, consensuslen, hitCnt))

blastp.close()
outfile.close()

##############

# this part of script adds the consensus protein identity to the summary file - xxx.table.

fileprefix = sys.argv[1][:-4]
infile = open(fileprefix+"_newBestRep.txt", 'rU')
oldsummary = open(sys.argv[2], 'rU')
newsummary = open(sys.argv[2]+".updated_wConsensusLength", 'w')

clstrIDdict = {}

for line in infile:
	line = line.strip()
	col = line.split("\t")
	key = col[0].split("_")
	clstr = key[0]
	clstrIDdict[clstr] = col[1] + "\t" + col[2] + "\t" + col[3]

# some representative sequences yielded no blast results (! e.g. Cluster6_Trinity|comp8954_c0_seeq3_1) using the stand-a-lone blast engine, and online BLAST site tells that "no putative conserved domains have been detected)

newsummary.write("# cluster\tcluster size\trepresentative seq\tseq length\tcontig count(Trinity)\tcontig count(Rockhopper)\tcontig count(IDBA)\tcontig count(Oases)\tcontig count(TransABYSS)\tconsensus protein ID\tmax length of the consensus protein\tnumber of matches\n")

for line in oldsummary:
	if line.startswith("#") == False:
		line = line.strip()
		col = line.split("\t")
		key = col[0]
		if key in clstrIDdict:
			newsummary.write("%s\t%s\n" % (line, clstrIDdict[key]))
		else:
			newsummary.write("%s\tblastp detected no putative conserved domains\n" % (line))

infile.close()
oldsummary.close()
newsummary.close()




