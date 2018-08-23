
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
    pos = {}
    for line in f:
        words = line.split()
        pos[ int(words[0]) ] = (words[1], words[2])
    f.close()
    print len(pos)
    return sorted(pos.items());

def is_same(a,b):
    return a == b or a == b.upper() or a.upper() == b

if __name__ == "__main__":

    #
    if(len(sys.argv) < 2):
        print "Usage: python " + sys.argv[0] + " chr1.fasta snp_position_filename "
        sys.exit()

    #pos = read_positions(sys.argv[1])
    record = SeqIO.read(open(sys.argv[1]),"fasta")
    ref = str(record.seq).upper() # only A/T/C/G no a/t/c/g
    refLen = len(ref)
    
    start = 0
    mutipleRange = []
    i=start+1
    # find AAAA/TTTT
    sumLen = 0
    while i<refLen:
        while i<refLen and ref[i] == ref[start]:
            i = i+1;
        l = i - start    
        if l >= 3:
            #fout.write("%d %d\n" % (start,i))
            mutipleRange.append((start,i))
            sumLen = sumLen + l
            #print start,i
            #print ref[start:i] #index, include start, not include i 
        start = i    
        i = start + 1

    print "number of single mutiple:", len(mutipleRange), sumLen

    #fout.write("number of single mutiple: %d\n" % len(mutipleRange))
    #fout.write("single multiple total length: %d\n" % sumLen)
    
    
    #start = 0
    sumLen = 0
    count = 0
    for start in range(0,2):
        i = start + 2
        while i<refLen:
            while i<refLen and ref[start]!=ref[start+1] and ref[start:start+2] == ref[i:i+2]:
                i = i+2
                
            l = i - start    
            if l >= 8:
                #fout.write("%d %d\n" % (start,i))
                count = count + 1
                sumLen = sumLen + l
                mutipleRange.append((start,i))
            start = i    
            i = start + 2

    print "number of double mutiple:", count, sumLen
            
    
    mutiple = sorted(mutipleRange)

    fout = open("chr1_multiple_pos","w")
    mutiplePos= set()
    for (i,j) in mutiple:
        if i>0:
            s = i-2
        else:
            s=i
        e = j+2
        for t in range(s,e+1):
            mutiplePos.add(t+1) #range is index, we want pos
        
    for c in mutiplePos:
        fout.write("%d\n" % (c))  
    fout.close()    
    
    '''
    pos = read_positions(sys.argv[2])
   

    i = 0
    count = 0
     
    fout = open("filter_single3_double8_multiple_delete1s","w")
    #fout = open("filter_single3_double8_multiple_insert1s","w")
    for (p, c) in pos:
        if p not in mutiplePos:        
            fout.write("%d %s %s\n" % (p, c[0], c[1]))
    fout.close()           
    '''
