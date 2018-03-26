#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
##import contig

###################
# map module
###################
def write_Map(fout, nuc):
    for (key,v) in nuc.items():
        fout.write("%s %s %s\n"%(key, len(v), v)) 
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

def reverse_Bool(s1):
    s2 = ""
    for c in s1:
        if c == '1':
            s2 = s2 + '0'
        elif c == '0':
            s2 = s2 + '1'  
        else:
            print ("error")  
    return s2


# s1 and s2 both string, only include '0'/'1' 
def is_Bool_Reverse(s1, s2):
    assert len(s1) == len(s2)
    if s1.find('2') != -1 or s2.find('2') != -1:
        return False
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
            print (i)
    return count

def similar_Distance(s1, s2):
    count = 0
    assert len(s1) == len(s2)
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            count +=1
    return count

def same_Character(s):
    first = s[0]
    for c in s:
        if c != first:
            return False
    return True    

# about range
def pos_2_Range(stablePos):

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
    return stableRange

def get_Range_From_List(small, large):
    

    #print ("small:", small)
    #print ("large:", large)
    assert set(small).intersection(large) == set(small)
    
    #for ele in large:
    s = large.index(small[0])
    e = large.index(small[-1])
    assert large[s:e+1] == small
    return (s,e+1)

def is_SubRange(i, j, stableRange):
    for (a,b) in stableRange:
        if i >= a and j <= b:
            return True
        if a >= j:
            return False
    return False

def get_Cover_Range(rangeList):

    # 3 is not cover
    s = 0
    e = 0
    i = 0
    for c in rangeList:
        if c != 3:
            s = i
            break
        i += 1
    while i < len(rangeList):
        if rangeList[i] == 3:
            e = i
            break    
        i += 1
    if e == 0:
        e = len(rangeList)
    assert (rangeList[s:e].count(3)==0) 
    #print (s,e)
    return (s,e)

def count01(l):
    count = 0
    for c in l:
        if c == 0 or c==1:
            count += 1
    return count 

def only0(a,b,c,d):
    if a != 0 and b == 0 and c == 0 and d == 0:
        return 1

    if a == 0 and b != 0 and c == 0 and d == 0:
        return 2

    if a == 0 and b == 0 and c != 0 and d == 0:
        return 3

    if a == 0 and b == 0 and c == 0 and d != 0:
        return 4

    return 0
