#!/usr/bin/env python

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys


record = SeqIO.read(open(sys.argv[1]), "fasta")
 
print record.id



s = str(record.seq)
lenS = len(s)
count = 1
i=0
while i<lenS:
    # or s[i] == 'n'  now is fine, only block5 and block9 have single a 
    while i<lenS and s[i] == 'N':
        i=i+1;
    start = i;
    while i<lenS and s[i] != 'N':
        i=i+1
    end = i
    if end -start>= 1000000: 
        newRecord = SeqRecord(Seq(s[start:end]),id=str(count)+"_"+str(start+1))
        SeqIO.write(newRecord,'block' + str(count)+"_" + str(start+1) + "_" + str(end) + '.fasta','fasta')
        print "write block", count, start+1, end, end-start
        count = count+1






