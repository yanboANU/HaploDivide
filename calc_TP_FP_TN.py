#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import HaploDivideBasedonPysam



def read_Actual_Diff(filename):
    f = open(filename, "r")
    actual_diff = {}
    for line in f:
        words = line.strip().split(' ')
        
        actual_diff[int(words[0])] = [] 
        i = 1
        while i < len(words):
            actual_diff[int(words[0])].append(words[i]) 
            i += 1
    return actual_diff


def get_part(actual_diff, l, r):
    actual_diff_part = {}
    for (key, value) in actual_diff:
        if key > l and key <r:
            actual_diff_part[key] = value
    return actual_diff_part 
            

if __name__ == "__main__":

    # sorted bam
    samfile = pysam.AlignmentFile(sys.argv[1], "rb")
   
    ll = 30000
    rr = 50000
    # contig
    contigs = HaploDivideBasedonPysam.read_Contig(sys.argv[2])
    insert = HaploDivideBasedonPysam.get_Insert(samfile)
    stableRange, pileupcolumns, diff = HaploDivideBasedonPysam.pre_Process(samfile,contigs,insert,ll,rr)


    actual_diff = read_Actual_Diff(sys.argv[3]).items()
    actual_diff_part = get_part(actual_diff, ll, rr)
    
     
    print ("length of finding diff: ",len(diff))
    print (sorted(diff.items()))
    
    print ("length of actual diff: ", len(actual_diff_part))  
    print (sorted(actual_diff_part.items()))
    '''
    print ("length of TP: ",len( set(diff) & set(actual_diff_part) ) )
    print (sorted(set(diff) & set(actual_diff_part)))
    
 
    print ("length of FP:", len( set(diff) - set(actual_diff_part) ) )
    print (sorted(set(diff) - set(actual_diff_part)))


    print ("length of TN", len( set(actual_diff_part) - set(diff) ) )
    print (sorted(set(actual_diff_part) - set(diff)))
    '''
    TP = set()
    FP = set()
    TN = set()
    for key in diff:
        if ((key-1) in actual_diff_part) or (key in actual_diff_part) or (key+1 in actual_diff_part):
            TP.add(key)
        else:
            FP.add(key) 
    for key in actual_diff_part:
        if (key not in TP) and (key+1 not in TP) and (key-1 not in TP):
            TN.add(key) 
 
    print ("length of TP: ",len(TP) )
    print (sorted(TP))

 
    print ("length of FP:", len(FP) )
    print (sorted(FP))

    print ("length of TN", len(TN) )

    print (sorted(TN))
