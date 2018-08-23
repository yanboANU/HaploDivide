#!/usr/bin/python


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import random

snpRate = 0.001 #or Indel 1
chrLen = 243000000 #

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
        print "Usage: python generate_snp_pos chrId_snps snp_position_filename "
        sys.exit()

    pos = read_positions(sys.argv[1])
    record = SeqIO.read(open(sys.argv[2]),"fasta")
    ref = record.seq

    count = 0
    fout = open("good_simulate_insert1s","w")
    for (p, c) in pos:   # v1 is position, v1-1 is index
        index = p-1
        #print p, c, ref[index]
        assert is_same(c[0],ref[index])
        i=index-1

        #    i   ...   index, index+1  ... j
        #  not A ...     A       T         not T
        #                  c[1][1]
        while ref[index] == ref[i]:
            i=i-1
        beforeLen = index-i
        j=index+2
        while ref[index+1] == ref[j]:
            j=j+1
        afterLen = j-index-1
        if is_same(c[1][1], ref[index]) and is_same(c[1][1],ref[index+1]):
            totalLen = beforeLen + 1 + afterLen
            if totalLen >= 3:
                count = count + 1
                continue
                '''
                print p, c, ref[index]
                print i,index,j
                print ref[i:index+1]
                print ref[index+1:j+1]
                print "total len", totalLen 
                '''
        if is_same(c[1][1], ref[index]):
            beforeLen = beforeLen + 1
        if is_same(c[1][1], ref[index+1]):
            afterLen = afterLen + 1
        if beforeLen >=3 or afterLen >= 3:     
            count = count + 1
            '''
            print p, c, ref[index]
            print i,index,j
            print ref[i:index+1]
            print "before len", beforeLen
            print ref[index+1:j+1]
            print "after len", afterLen
            '''
        if beforeLen <3 and afterLen < 3:
            fout.write("%d %s %s\n" % (p,c[0],c[1]))
    fout.close()     
    print count

