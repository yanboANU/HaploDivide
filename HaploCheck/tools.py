#!/usr/bin/env/ python3

import sys
import os
#import pysam
import string
import copy

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

def convert(lToR, left):
    right = []
    not_find = [] 
    for ll in left:
        if ll in lToR:
           right.append(lToR[ll])
        else:
           not_find.append(ll)
    return right, not_find




# s1 and s2 both string, only include '0'/'1' 
def is_Bool_Reverse(s1, s2):
    assert len(s1) == len(s2)
    if s1.find('2') != -1 or s2.find('2') != -1:
        return False
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            return False
    return True    

def bool_Reverse(s1):
    s2 = ""

    for c in s1:
        if c == "1":
            s2 = s2 + "0"
        elif c == "0":
            s2 = s2 + "1"
        else:
            s2 = s2 + c
    return s2


def hamming_Distance(s1, s2):
    count = 0
    assert len(s1) == len(s2)
    diffPos = []
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count +=1
            diffPos.append(i)
    #print ("diff pos:", diffPos)       
    return count, diffPos

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

# longest_repeating_strings(follow 3)
def partition(suffix_array, start, end):
    if end <= start:
        return
    index1, index2 = start, end
    base = suffix_array[start]
    while index1 < index2 and suffix_array[index2] >= base:
        index2 -= 1
    suffix_array[index1] = suffix_array[index2]
    while index1 < index2 and suffix_array[index1] <= base:
        index1 += 1
    suffix_array[index2] = suffix_array[index1]
    suffix_array[index1] = base
    partition(suffix_array, start, index1 -  1)
    partition(suffix_array, index1 + 1, end)

def find_common_string(str1, str2):
    if not str1 or not str2:
        return 0, ''
    index1, index2 = 0, 0
    length, comm_substr = 0, ''
    while index1 < len(str1) and index2 < len(str2):
        if str1[index1] == str2[index2]:
            length += 1
            comm_substr += str1[index1]
        else:
            break
        index1 += 1
        index2 += 1
    return length, comm_substr

def find_longest_repeating_strings(string):
    if not string:
        return None, None
    suffix_array = []
    # first, get the suffix arrays
    length = len(string)
    for i in range(length):
        suffix_array.append(string[i:])
    # second, sort suffix array
    start, end = 0, len(suffix_array) - 1
    partition(suffix_array, start, end)
    # third, get the longest repeating substring
    max_length,  repeat_substring = 0, ''
    for i in range(len(suffix_array) - 1):
        common_len, common_substring = find_common_string(suffix_array[i], suffix_array[i+1])
        if common_len > max_length:
            max_length, repeat_substring = common_len, common_substring
    return max_length, repeat_substring




def find_lcsubstr(s1, s2):
    p = 0 # end position of s1
    q = 0 # end position of s2
    mmax = 0
    m = [ [ 0 for i in range(len(s2)+1)] for j in range(len(s1)+1) ] 
    for i in range(len(s1)):  
        for j in range(len(s2)):  
            if s1[i]==s2[j]:  
                m[i+1][j+1]=m[i][j]+1  
                if m[i+1][j+1]>mmax:  
                    mmax = m[i+1][j+1]  
                    p = i+1  
                    q = j+1
    return s1[p-mmax:p],mmax,p,q  
  



