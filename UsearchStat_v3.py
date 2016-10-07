#!/usr/bin/python2.6

# UsearchStat_v3.py

# Seqs were searched against LSU, SSU, 5S and tRNA databases using USEARCH. Using the default search criteria and id=0.8, it is returned with tables listing the query ID, the subject ID and the matching info (e.g. % id). The query IDs are believed to be related to LSU, SSU, 5S and tRNA. This script is to summarize how many reads are suggested to be LSU, SSU, 5S and tRNA, and which are they. 
# There are multiple hits per each query read!

# define usage message
Usage = """UsearchStat_v3.py [USEARCH 4LSU result table] [USEARCH 4SSU result table] [USEARCH 5S result table] [USEARCH 4tRNA result table] > [readsfilename_usearchsummary.csv]"""

import sys

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) <= 2:
    print Usage
    exit(1)

LSUtable = open(sys.argv[1], 'rU') # open the file as a text file
SSUtable = open(sys.argv[2], 'rU') # open the file as a text file
SrRNAtable = open(sys.argv[3], 'rU') # open the file as a text file
tRNAtable = open(sys.argv[4], 'rU') # open the file as a text file

list = []
LSUlist = []
SSUlist = []
SrRNAlist = []
tRNAlist = []

def hitscount(x):
	listdict = {}
	for line in x:
		col = line.split("\t")
		key = col[0]
		if key in listdict:
			listdict[key] +=1 # same as listdict[key] = listdict[key] + 1
		else:
			listdict[key] = 1
	return listdict

LSUdict = hitscount(LSUtable)
SSUdict = hitscount(SSUtable)
SrRNAdict = hitscount(SrRNAtable)
tRNAdict = hitscount(tRNAtable)

names = set(LSUdict.keys() + SSUdict.keys() + SrRNAdict.keys() + tRNAdict.keys())
#maxlen = max([len(n) for n in names])

print '\%s (%d)\tLSU (%d)\tSSU (%d)\tSrRNA (%d)\ttRNA (%d)' % ('read', len(names), len(LSUdict.keys()), len(SSUdict.keys()), len(SrRNAdict.keys()), len(tRNAdict.keys()))
#print "%*s\t%d\t%d\t%d" % (maxlen, 'total # reads', len(LSUdict.keys()), len(SSUdict.keys()), len(SrRNAdict.keys()), len(tRNAdict.keys()))

for n in names:
	print '%s\t%d\t%d\t%d\t%d' % (n, LSUdict.get(n,0), SSUdict.get(n,0), SrRNAdict.get(n,0), tRNAdict.get(n,0))
 

LSUtable.close()
SSUtable.close()
SrRNAtable.close()
tRNAtable.close()
