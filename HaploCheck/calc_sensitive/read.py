import sys
from itertools import ifilter,imap
import collections
import string



def read_pre_snp(filename):
#def read_snp(filename):
    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    snpSupportNum = {}
    for line in f: 
        #print line
        if line.startswith('#'):
            continue
        words = line.strip().split()
        if len(words) == 5: #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            snpContent[pos] = (words[1], words[3])
            snpSupportNum[pos] = (int(words[2]), int(words[4]))
        '''
        if len(words) == 3: #for file mutation record
            pos = int(words[0])
            snpPosition.add(pos)
            snpContent[pos] = (words[1], words[2])
        '''    

    return snpPosition, snpContent, snpSupportNum


def read_delete(filename, start, end, base):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    #requre snp store in order
    for line in f:  
        words = line.strip().split()
        if len(words) == 3: #for file mutation record
            pos = int(words[0])
            if pos < start:
                continue
            if pos > end:
                break
            snpPosition.add(pos+base)
            snpContent[pos+base] = (words[1], words[2])
            deleteLen = len(words[1])-1
            count = 1
            if deleteLen > 1:
                snpPosition.add(pos+base+count)
                count += 1
                deleteLen -= 1
    return snpPosition, snpContent

def read_snp2(filename, start, end, base):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    #requre snp store in order
    for line in f:  
        if line.startswith('#'):
            continue
        words = line.strip().split()
        if len(words) == 5: #for file contig_snp_mutation
            pos = int(words[0])
            if pos < start:
                continue
            if pos > end:
                break
            snpPosition.add(pos+base) 
            snpContent[ pos+base ] = (words[1], words[3])
        if len(words) == 4: #for file mutation record
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
    count = 0
    for line in f:
        if line.startswith("p"):
            continue
        words = line.strip().split(" ")
        #print words
        #sys.exit()
        coverage[int(words[0])] = int(words[1])
        if int(words[1]) >= 8:
            count += 1
    return coverage, count


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
