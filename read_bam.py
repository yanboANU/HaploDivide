#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import contig


#input : bam 

if __name__ == "__main__":

    samfile = pysam.AlignmentFile(sys.argv[1], "rb")

    count = 0

    print ("read name, read align length, read len, mapping rate")
    for read in samfile.fetch():
        count += 1
        readName = read.query_name
        readMapScore = read.mapping_quality
        mapped = float(read.qlen)/len(read.query_sequence)
        print (readName, read.qlen, len(read.query_sequence), mapped)   
    print("reads number:", count)
    samfile.close()

