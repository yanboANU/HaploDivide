#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import contig

###################
# map module
###################
def write_Map(fout, nuc):
    for (key,v) in sorted(nuc.items()):
        fout.write("%s %s\n"%(key, v)) 
    fout.write("\n")  

# sorted value
def sorted_Map_Value(m, R=True):
    sortedM = []
    for k, v in [(k, m[k]) for k in sorted(m, key=m.get, reverse=R)]:
        sortedM.append((k,v))  
    return sortedM


# the value of map is a list
# sort the length of the list
def sorted_Map_Value_Len(m, R=True):
    sortedM = sorted(m.items(), key=lambda e:len(e[1]), reverse=True)
    #for k, v in [(k, m[k]) for k in sorted(m, key=m.get, reverse=R)]:
    #    sortedM.append((k,m[k]))  
    return sortedM



# s1 and s2 both string, only include '0'/'1' 
def is_Bool_Reverse(s1, s2):
    assert len(s1) == len(s2)
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            return False
    return True    

def hamming_Distance(s1, s2):
    count = 0
    assert len(s1) == len(s2)

    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count +=1
    return count

# about range
def pos_2_Range(pos):

    stableRange = []
    
    i = 0
    while i < len(stablePos)-1:
        s = stablePos[i]
        e = stablePos[i]
        while (stablePos[i]+1 == stablePos[i+1]):
            i = i+1
            e = stablePos[i]
            if (i >= len(stablePos)-1):
                break 
        i = i+1
        stableRange.append((s,e)) 
    return pos_2_Range 


def is_SubRange(i, j, stableRange):
    for (a,b) in stableRange:
        if i >= a and j <= b:
            return True
        if a >= j:
            return False
    return False
