#!/usr/bin python

# reformat_cdhit_clstr_v3_ncRNA_singlesample.py
# Created by Maggie Lau (maglau@princeton.edu)
# Created on January 30, 2016



import sys
from Bio import SeqIO

# define usage message
Usage =  """Usage: reformat_cdhit_clstr_v3_ncRNA_singlesample.py [.clstr] [Representative sequences] [chimeraslist]"""


#   Check to make sure that the correct number of arguments are given
if len(sys.argv) <= 3:
    print Usage
    exit(1)


####### PART 1 #######

# with the output from cd-hit-est (.clstr) as input, information will be summarized to show the Cluster, corresponding representative seq title, size of each cluster, constituent samples for each cluster

clstr = open(sys.argv[1], 'rU') # open the file as a text file
fileprefix = sys.argv[1][:-6]
outtable = open(fileprefix+".tableTrimmed", 'w')
oldfasta = open(sys.argv[2], 'rU')
newfasta = open(sys.argv[2]+".renameTrimmed", 'w')
chimeras = open(sys.argv[3], 'rU')

prevcluster = None
cluster = None
clstrsize = 0
RepSeq = None
RepSeqLength = 0
clstrdict = {}

outtable.write("# cluster\tcluster size\trepresentative seq\tseq length\t\n")

chiDict ={}

for line in chimeras:
	line=line.strip()
	if line.startswith(">"):
		title = line[1:]
		chiDict[title] = 1

for line in clstr:
        line = line.strip()
	if line.startswith(">"):
		title = line[1:].split()
		cluster = title[0] + title[1]
		if cluster != prevcluster and prevcluster != None and RepSeq not in chiDict:
			outtable.write("%s\t%d\t%s\t%s\n" % (prevcluster, clstrsize, RepSeq, RepSeqLength))
		clstrsize = 0
		prevcluster = cluster
	else:
		col = line.split()
		seqtitle = col[2][1:-3]
		if seqtitle not in chiDict:
			clstrsize += 1
              	if '*' in line:
                        RepSeq = seqtitle
                        RepSeqLength = col[1][:-3]
			clstrdict[RepSeq] = cluster

# this is to print out the data for the last cluster
outtable.write("%s\t%d\t%s\t%s\n" % (cluster, clstrsize, RepSeq, RepSeqLength))


for each in SeqIO.parse(oldfasta, 'fasta'):
	if each.id not in chiDict:
		oldseqtitle = each.id
		newseqtitle = clstrdict[oldseqtitle] + "_" + oldseqtitle
		seq = each.seq
		newfasta.write(">%s\n" % (newseqtitle))	
		newfasta.write("%s\n" % (seq))
		
clstr.close()
outtable.close()
oldfasta.close()
newfasta.close()
chimeras.close()

####### PART 2 #######

# with the clstr file remains open

# reformat it to produce a look-up table for each sequence name and its corresponding represence sequence name, and perhaps the % of identity

clstr = open(sys.argv[1], 'rU')
columntable = open(fileprefix+".clstrCOLtrimmed", 'w') # clster info in column format

repseq = None
prevclstrname = None
infolist = []

for line in clstr:
        line = line.strip()
        if line.startswith(">"):
		col = line.split()
                clstrname = ''.join(col)
                if infolist != None and repseq not in chiDict:
			for i in infolist:
                        	columntable.write('%s\t%s\t%s\n' % (prevclstrname, repseq, '\t'.join(i)))
                infolist = []
                prevclstrname = clstrname
        else:
                col = line.split()
                seqlen = col[1][:-3]
                title = col[2][1:-3]
                if '*' not in line:
                        aln = col[4].split("/")
                        alnpos = aln[0]
			strand = aln[1]
                        percent = aln[2]
			if title in chiDict:
                        	info = (title, seqlen, alnpos, strand, percent, 'CHIMERA')
                        	infolist.append(info)
			else:
				info = (title, seqlen, alnpos, strand, percent)
                                infolist.append(info)
                else:
                        	repseq = title
                        	info = (repseq, seqlen, 'RefSeq', 'RefSeq', 'RefSeq')
                        	infolist.append(info)

for i in infolist:
        columntable.write('%s\t%s\t%s\n' % (prevclstrname, repseq, '\t'.join(i)))



clstr.close()
columntable.close()

