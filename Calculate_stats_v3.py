#! /usr/bin/env python2.7

# This script will generate a few files:
# _length.txt that contains the name of contig and its length
# _contig_stat.txt that contains the basic stats of the input contig file
# _converage.txt that contains the name of contig, its length, number of reads mapped, and coverage 


## This version put the "average read length" as an argument in the execution command
# Usage:
# $ python Calculate_stats_v3.py contig-fasta-file n prefix-of-the-"_mappedReads.bed"-file

import sys
from Bio import SeqIO
import collections

contigs = open(sys.argv[1], 'rU')
out = open(sys.argv[3] + "_length.txt", 'w')
stat = open(sys.argv[3] + "_contig_stat.txt", 'w')

seqlen = []
N = 0

for rec in SeqIO.parse(contigs, 'fasta'):
	out.write(rec.id+'\t'+str(len(rec.seq))+'\n')
	seqlen.append(len(rec.seq))

for i in seqlen:
	N += 1  

avseqlen = sum(seqlen)/N

N50 = 0
N50con = 0
testsum = 0 

for i in sorted(seqlen, reverse=True):
	testsum += i
	N50con += 1
	if sum(seqlen)/2.0 < testsum:
		N50 = i
		break # critical 

N90 = 0
N90con = 0
testsum = 0

for i in sorted(seqlen, reverse=True):
	testsum += i
	N90con += 1
	if sum(seqlen)*0.9 < testsum:
		N90 = i
		break # critical

stat.write('input file: %s\n' % (sys.argv[1]))
stat.write('number of contigs: %d\n' % (N))
stat.write('number of bases (bp): %d\n' % (sum(seqlen)))
stat.write('max contig length (bp): %d\n' % (max(seqlen)))
stat.write('min contig length (bp): %d\n' % (min(seqlen)))
stat.write('average contig length (bp): %d\n' % (avseqlen))
stat.write('N50 length (bp): %d \n' % (N50))
stat.write('N50 contig: %d\n' % (N50con))
stat.write('N90 length (bp): %d\n' % (N90))
stat.write('N90 contig: %d\n' % (N90con))



contigs.close()
out.close()
stat.close()

##############################################################

mapresult = open(sys.argv[3] + "mappedReads.bed", 'rU')
contiglen = open(sys.argv[3] + "_length.txt", 'rU')
out2 = open(sys.argv[3] + "_coverage.txt", 'w')

avreadlen = int(sys.argv[2])

contiglenD = {}
mapD = collections.Counter()


for l in contiglen:
        l = l.strip()
        l = l.split("\t")
        contiglenD[l[0]] = int(l[1])
#       print "%s\t%s\n" % (l[0], contiglenD[l[0]])
        
for line in mapresult:
        line = line.strip()
        line = line.split("\t")
        mapD[line[0]] += 1

for k,v in sorted(mapD.items(), key=lambda x: x[1], reverse=True):
#for key in sorted(mapD, key=lambda x: x[1], reverse=True):
#for key in sorted(mapD, key=lambda x: int(x.split('_')[1])):
        C = (mapD[k]*avreadlen)/float(contiglenD[k])
        out2.write("%s\t%d\t%d\t%.2f\n" % (k, contiglenD[k], mapD[k], C))

mapresult.close()
contiglen.close()
out2.close()



