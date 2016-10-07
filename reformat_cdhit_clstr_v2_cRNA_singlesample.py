#!/usr/bin python

# reformat_cdhit_clstr_v2_cRNA_singlesample.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on January 30, 2016



import sys
import re
#from collections import Counter


# define usage message
Usage =  """Usage: reformat_cdhit_clstr_v2_cRNA_singlesample.py [.clstr] [Representative sequences]"""


#   Check to make sure that the correct number of arguments are given
if len(sys.argv) != 3:
    print Usage
    exit(1)


####### PART 1 #######

# with the output from cd-hit (.clstr) as input, information will be summarized to show the Cluster, corresponding representative seq title, size of each cluster

clstr = open(sys.argv[1], 'rU') # open the file as a text file
fileprefix = sys.argv[1][:-6]
outtable = open(fileprefix+".table", 'w')
oldfasta = open(sys.argv[2], 'rU')
newfasta = open(sys.argv[2]+".rename", 'w')

prevcluster = None
cluster = None
clstrsize = 0
RepSeq = None
RepSeqLength = 0
clstrdict = {}

outtable.write("# cluster\tcluster size\trepresentative seq\tseq length\n")


for line in clstr:
        line = line.strip()
	if line.startswith(">"):
		title = line[1:].split()
		cluster = title[0] + title[1]
		if cluster != prevcluster and prevcluster != None:
			outtable.write("%s\t%d\t%s\t%s\n" % (prevcluster, clstrsize, RepSeq, RepSeqLength))
		clstrsize = 0
		prevcluster = cluster
		
	else:
		clstrsize += 1
		col = line.split()
                seqtitle = col[2][1:-3]
              	if '*' in line:
                        RepSeq = seqtitle
                        RepSeqLength = col[1][:-3]
			clstrdict[RepSeq] = cluster

# this is to print out the data for the last cluster
outtable.write("%s\t%d\t%s\t%s\n" % (cluster, clstrsize, RepSeq, RepSeqLength))


for l in oldfasta:
	l = l.strip()
	if l.startswith(">"):
		col = l.split()
		oldseqtitle = col[0][1:]
		newseqtitle = clstrdict[oldseqtitle] + "_" + oldseqtitle
		newfasta.write(l.replace(oldseqtitle, newseqtitle))
		newfasta.write("\n")
	else:
		newfasta.write("%s\n" % (l))
		
clstr.close()
outtable.close()
oldfasta.close()
newfasta.close()


####### PART 2 #######

# with the clstr file remains open

# reformat it to produce a look-up table for each sequence name and its corresponding representative sequence name, and perhaps the % of identity

clstr = open(sys.argv[1], 'rU')
columntable = open(fileprefix+".clstrCOL", 'w') # clster info in column format

prevclstrname = None
infolist = []

for line in clstr:
        line = line.strip()
        if line.startswith(">"):
                col = line.split()
                clstrname = ''.join(col)
		if infolist != None:
                	for i in infolist:
                        	columntable.write('%s\t%s\t%s\n' % (prevclstrname, refseq, '\t'.join(i)))
                	infolist = []
                prevclstrname = clstrname
        else:
                col = line.split()
                PEGtitle = col[2][1:-3]
                PEGlen = col[1][:-3]
                parentcontig = col[2][1:-5]
                if '*' not in line:
                        aln = col[4].split("/")
                        alnpos = aln[0]
                        percent = aln[1]
                        info = (PEGtitle, PEGlen, parentcontig, alnpos, percent)
                        infolist.append(info)
                else:
			refseq = PEGtitle
                        info = (PEGtitle, PEGlen, parentcontig, 'RefSeq', 'RefSeq')
                        infolist.append(info)

for i in infolist:
        columntable.write('%s\t%s\t%s\n' % (prevclstrname, refseq, '\t'.join(i)))



clstr.close()
columntable.close()

