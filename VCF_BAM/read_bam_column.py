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

    #checkPos = 2142505 + 1
    checkPos = int(sys.argv[2])
    for pileupcolumn in samfile.pileup("5_29553836", checkPos, checkPos+5): 
        pp = pileupcolumn.pos
        cov = pileupcolumn.n
        
        #if pp > checkPos-5 and pp<checkPos+5:    
        if pp == checkPos:
            print (pp, cov)
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # query position is None if is_del or is_refskip is set.
                    print ('\tbase in read %s = %s length=%s aligned_length=%s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position], pileupread.alignment.query_length, pileupread.alignment.query_alignment_length))
            
    samfile.close()

