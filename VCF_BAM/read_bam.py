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
    coverage = {}
    #foutD = open("plot_cov","w")
    for pileupcolumn in samfile.pileup(sys.argv[2]):
        #print ("%s %s" % (pileupcolumn.pos, pileupcolumn.n))
        pos = pileupcolumn.pos
        cov = pileupcolumn.n
        if cov not in coverage:
            coverage[cov] = 1
        else:
            coverage[cov] += 1
            #if pos%1000000 == 0:
                #foutD.write("%s %s\n" % (pos, cov))
        if cov >= 7:
            count += 1
    #foutD.close()
    print ("number of coverage more than 7", count)
    fout = open("cov_count","w")
    #for (cov, count) in coverage.iteritems(): python2
    coverageLen = len(coverage)
    i=0
    while i < coverageLen:
        cov, count = coverage.popitem()
        fout.write("%s %s\n" % (cov, count))
        i += 1
    fout.close()
    
    '''
    print ("read name, read align length, read len, mapping rate")
    
    for read in samfile.fetch():
        
        if count == 10:
            break
            
        count += 1
        readName = read.query_name
        referenceName = read.reference_name
        readMapScore = read.mapping_quality
        mapped = float(read.qlen)/len(read.query_sequence)
        print (readName, referenceName, read.qlen, len(read.query_sequence), mapped)   
        sys.exit()
    print("reads number:", count)
    '''
    samfile.close()

