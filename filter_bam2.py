
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

    newBam = pysam.AlignmentFile("2_filter.bam", "wb", template=samfile)
    count = 0
    for read in samfile.fetch():
        #readName = read.query_name
        #readMapScore = read.mapping_quality
        mapped = float(read.qlen)/len(read.query_sequence)
        if mapped >= 0.5:
            count += 1
            newBam.write(read)

    print ("after 2 filter, the number of reads:", count)
    newBam.close()
    samfile.close()
