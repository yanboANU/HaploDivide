import sys
from itertools import ifilter,imap
import collections
import string



def read_snp(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    for line in f:  
        words = line.strip().split()
        if len(words) >= 5: #for file contig_snp_mutation
            snpPosition.add(int(words[0])) 
            snpContent[int(words[0])] = (words[1], words[3])
        if len(words) == 3: #for file mutation record
            snpPosition.add(int(words[0]))
            snpContent[int(words[0])] = (words[1], words[2])

    return snpPosition, snpContent


def read_snp2(filename, start, end, base):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    #requre snp store in order
    for line in f:  
        words = line.strip().split()
        if len(words) == 5: #for file contig_snp_mutation
            
            pos = int(words[0])
            '''
            if pos < start or pos > end :
                continue
            '''
            if pos < start:
                continue
            if pos > end:
                break
            snpPosition.add(pos+base) 
            snpContent[ pos+base ] = (words[1], words[3])
        if len(words) == 3: #for file mutation record
            pos = int(words[0])
            if pos < start:
                continue
            if pos > end:
                break
            snpPosition.add(pos+base)
            snpContent[pos+base] = (words[1], words[2])

    return snpPosition, snpContent


def read_cov(filename):
    coverage = {}
    f = open(filename, "r")
    for line in f:
        if line.startswith("p"):
            continue
        words = line.strip().split(" ")
        coverage[int(words[0])] = int(words[1])
    return coverage


def read_multiple_pos(filename, start, end, base):
    multiplePos = set() 
    f = open(filename, "r")
    for line in f:
        words = line.strip().split(" ")
        val = int(words[0])  
        if val < start:
            continue
        if val > end:
            break
        multiplePos.add(val+base)
    return multiplePos    
