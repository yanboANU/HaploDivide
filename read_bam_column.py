#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import contig


#input : bam 

if __name__ == "__main__":

    samfile = pysam.AlignmentFile(sys.argv[1], "rb")

    print ("position, coverage")

    for pileupcolumn in samfile.pileup(): 
        pp = pileupcolumn.pos
        cov = pileupcolumn.n  
        print (pp, cov)
    samfile.close()

