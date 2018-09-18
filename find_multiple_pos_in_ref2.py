
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import random

def is_same(a,b):
    return a == b or a == b.upper() or a.upper() == b


def read_positions(filename):
    f = open(filename, "r")
    pos = []
    for line in f:
        words = line.split()
        pos.append( (int(words[0]),  words[1], words[2]) )
    f.close()
    return pos
    #return sorted(pos.items());

def is_same(a,b):
    return a == b or a == b.upper() or a.upper() == b

if __name__ == "__main__":

    #
    if(len(sys.argv) < 2):
        print "Usage: python " + sys.argv[0] + " chr1.fasta snp_position_filename "
        sys.exit()

    record = SeqIO.read(open(sys.argv[1]),"fasta")
    ref = str(record.seq).upper() # only A/T/C/G no a/t/c/g
    
    pos = read_positions(sys.argv[2])
    count = 0
    pos2 = []
    for (p, s1, s2) in pos:
        #print p,s1,s2, ref[p-1:p+2]
        assert ref[p-1] == s1[0]
        #if ref[p]==ref[p+1] or ref[p] == ref[p-1] or ref[p-1] == ref[p+1]: #for delete
        if ref[p] == s2[1] or s2[1]==ref[p-1] or ref[p] == ref[p-1]: #for insert
            count += 1 
        else:
            pos2.append((p, s1, s2))
    print "first homo:", count 
    #sys.exit()
    
    count = 0
    for (p, s1, s2) in pos2:
        seg = ref[p-3:p+4]
        flag = False
        for i in range(4):
            if seg[i:i+2] == seg[i+2:i+4]:
                count += 1
                flag = True
                break
        if flag == False:
            print p, s1, s2 # length 7

    print "second homo:", count    
