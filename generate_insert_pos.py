
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import random



singleMultipleThresdhold = 3

DoubleMultipleThresdhold = 8


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
        if l >= singleMultipleThresdhold:
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
            if l >= DoubleMultipleThresdhold:
                #fout.write("%d %d\n" % (start,i))
                count = count + 1
                sumLen = sumLen + l
                mutipleRange.append((start,i))
            start = i    
            i = start + 2

    print "number of double mutiple:", count, sumLen
            
    
    mutiple = sorted(mutipleRange)

    fout = open("chr1_multiple_range","w")
    for (i,j) in mutiple:
        fout.write("%d %d \n" % (i,j))
    fout.close()    

    
    
    pos = read_positions(sys.argv[2])
   
    lenMutiple = len(mutiple) 
    i=0
    while i+1<lenMutiple:
        if mutiple[i][1] > mutiple[i+1][0]:
            print mutiple[i], mutiple[i+1]

        assert mutiple[i][1] <= mutiple[i+1][0]
        i = i+1

    '''
    i = 0
    count = 0
     
    fout = open("filter_single3_double8_multiple_insert1s","w")
    for (p, c) in pos:
        while p > mutiple[i][0]:
            i = i+1
        if p > mutiple[i-1][1]:
            dis = min(p-mutiple[i-1][1], mutiple[i][0]-p)
            if dis >= 3:
                count = count + 1
                print mutiple[i-1], mutiple[i], dis
                print ref[p-5:p+5]
                print p, c[0], c[1]
                fout.write("%d %s %s\n" % (p, c[0], c[1]))
    fout.close()           
    '''
