#!/usr/bin/env python

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys


record = SeqIO.read(open(sys.argv[1]), "fasta")
 
print record.id

seq2 = re.sub('[NMR]', '', str(record.seq))

#print collections.Counter(seq2)
rec2 = SeqRecord(Seq(seq2),id=record.id)

SeqIO.write(rec2,'N_removed_chr19.fasta','fasta') 




