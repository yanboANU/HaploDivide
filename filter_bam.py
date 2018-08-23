#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import contig


#input : 0_16.bam
#every read need be covered 0.5*read_length 

if __name__ == "__main__":

    #input: *bam chrId
    samfile = pysam.AlignmentFile(sys.argv[1], "rb")
  

    newBam = pysam.AlignmentFile("chr1.bam", "wb", template=samfile) 
    count = 0
    chrID = sys.argv[2]
    for read in samfile.fetch(chrID):
        #readName = read.query_name
        #readMapScore = read.mapping_quality
        #mapped = float(read.qlen)/len(read.query_sequence)
        #print (read.flag)
        #if mapped >= 0.5 and (read.flag==0 or read.flag==16):
        count += 1
        newBam.write(read)

    print ("after 2 filter, the number of reads:", count)   
    samfile.close()
    newBam.close()

