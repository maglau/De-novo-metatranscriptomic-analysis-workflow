#!/usr/bin python2.6

# ParseSeqsIntoGroups_v3.py
# Created by Maggie Lau (maglau@princeton.edu) and modified by Dennis McRitchie

# Seqs were searched against LSU, SSU, 5S and tRNA databases using USEARCH. Using the default search criteria and id=0.8, it is returned with tables listing the query ID, the subject ID and the matching info (e.g. % id). The query IDs are believed to be related to LSU, SSU, 5S or tRNA. Use UsearchStat.py to generate a summary of all these files. This script is to grasp the corresponding seqeunces and write them into a new file as "rRNA" while the rest as "non-rRNA".

# define usage message
Usage =  """Usage: ParseSeqsIntoGroups_v4.py [USEARCH summary table] [reads.fasta] [outfile prefix]"""

import sys
from Bio import SeqIO

#   Check to make sure that the correct number of arguments are given
if len(sys.argv) != 4:
    print Usage
    exit(1)

usearchfile = open(sys.argv[1], 'rU') # open the file as a text file
fastafile = open(sys.argv[2], 'rU') # open the file as a text file
fileprefix = sys.argv[3]
ncRNAreads = open (fileprefix+"_ncRNA.fasta", 'w')
cRNAreads = open (fileprefix+"_cRNA.fasta", 'w')

idlist = {}

for line in usearchfile:
        col = line.split("\t")
        if len(col) == 0:
            continue
        if col[0] not in idlist:
                idlist[col[0]] = 1

for each in SeqIO.parse(fastafile, "fasta"):
	if each.description in idlist:
               SeqIO.write(each, ncRNAreads, 'fasta')
#                print("'%s' matches" % each.description)
        else:
                SeqIO.write(each, cRNAreads, 'fasta')
#                print("'%s' does not match" % each.description)

usearchfile.close()
fastafile.close()
cRNAreads.close()
ncRNAreads.close()
