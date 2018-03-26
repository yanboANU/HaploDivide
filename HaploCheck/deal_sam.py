#!/usr/bin/env python3

import sys
import pysam

def filter_SAM(filename, fileout):
    f = open(filename, "r")
    fout = open(fileout, "w")
    for line in f:
        '''
        if not line.startswith("read"):
            fout.write(line)
        else:
            words = line.split()
            if words[1] == "0" or words[1] == "16":
                fout.write(line)
        ''' 
        if line.startswith("read") or line.startswith('S'):
            words = line.split()
            if words[1] == "0" or words[1] == "16":
                fout.write(line)
        else:
            fout.write(line)
    f.close()
    fout.close()

#get align position from sam
def get_align_pos(samfilename):

    samfile = pysam.AlignmentFile(samfilename, "rb")
    align = {}
    for read in samfile.fetch(): 
        #(x,y) (contigpos ,refpos)
        alignedPairs = read.get_aligned_pairs()
        readName = read.query_name
        print (readName, len(alignedPairs))
        #print (alignedPairs)
        align[readName] = alignedPairs
    return align
 

if __name__ == "__main__":

    #filter_SAM(sys.argv[1], sys.argv[2])

    get_align_pos(sys.argv[1])
