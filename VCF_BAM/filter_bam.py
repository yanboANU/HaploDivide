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
  

    newBam = pysam.AlignmentFile("filter.bam", "wb", template=samfile) 
    count = 0
    rate = float(sys.argv[2])
    #chrID = sys.argv[2]
    minScore = int(sys.argv[3])
    for read in samfile.fetch():
        #readName = read.query_name
        readMapScore = read.mapping_quality
        mapped = float(read.qlen)/len(read.query_sequence)
        #print (read.flag)
        if mapped >= rate and (read.flag==0 or read.flag==16) and readMapScore >= minScore:
            count += 1
            newBam.write(read)
            #print (readMapScore, read.qlen, len(read.query_sequence), mapped)
            #if count == 100:
                #break

    print ("filter rate:", rate)   
    print ("minScore:", minScore)
    print ("after 2 filter, the number of reads:", count)   
    samfile.close()
    newBam.close()

