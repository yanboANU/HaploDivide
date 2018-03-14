#!/usr/bin/env python

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re



record = SeqIO.read(open("sequence.fasta"), "fasta")
 
print record.id

seq2 = re.sub('[NMR]', '', str(record.seq))

print collections.Counter(seq2)
rec2 = SeqRecord(Seq(seq2[0:10000000]),id=record.id)

SeqIO.write(rec2,'N_removed_first_10M.fasta','fasta') 




