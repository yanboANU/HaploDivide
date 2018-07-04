import sys
import sequence
import alignment
from itertools import ifilter,imap
import collections
import string

def read_file_list(filename):
    
    f = open(filename, "r")
    fileList = []
    for line in f:
        fileList.append(line.strip())

    return fileList  




def read_align(filename):
    # read align_pos 
    f = open(filename, "r")
    ref_pos = []
    contig_pos = []
    ref_label = 1 
    for line in f:
        words = line.split(',')
        if words[0].isdigit():
            if len(ref_pos) == 0:
                for c in words:
                    ref_pos.append(int(c))   
            else: 
                for c in words:
                    contig_pos.append(int(c))
    
    assert len(ref_pos) == len(contig_pos)
    '''
    fout = open("ref_contig", "w")  
    for i in range(len(ref_pos)):
        fout.write("%s %s\n" % (ref_pos[i], contig_pos[i]))
    '''
    return ref_pos, contig_pos

def read_mutation_list(filename):
    f = open(filename, "r")

    pos = []  
    for line in f:  
        words = line.split()
        if len(words) == 1 or len(words) == 3:
            pos.append(int(words[0]))
    
    return pos

def read_snp(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    for line in f:  
        words = line.strip().split()
        if len(words) == 5: #for file contig_snp_mutation
            snpPosition.add(int(words[0])) 
            snpContent[int(words[0])] = (words[1], words[3])
        if len(words) == 3: #for file mutation record
            snpPosition.add(int(words[0]))
            snpContent[int(words[0])] = (words[1], words[2])

    return snpPosition, snpContent


def read_phasing_result(filename):
    ##########
    #unfinish only read one segment
    ##########
    f = open(filename, "r")
    haplos = [] 
    lineNumber = 0
    for line in f:
        haplotype = {}
        lineNumber += 1
        if lineNumber % 6 == 0: 
            words = line.strip().split(',')
            #print (words)
            #print (len(words))    
            continue
        elif lineNumber % 6 == 1 and lineNumber >1:
            binarySeq = line.strip()
        else:
            continue 
        print ("len binnary",len(binarySeq))
        print ("len position",len(words))
        assert len(binarySeq) == len(words)
        for i in range(len(words)):
            haplotype[int(words[i])] = binarySeq[i] 
        haplos.append(haplotype) 
    return haplos
 

def read_blasr_m5(fileName):

    f = open(fileName,"r")
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        words = line.split()

        if len(words) <= 4:
            continue
        #print words[:4] 
        queryName, queryLen, queryS, queryE = words[:4]
        queryDirection = words[4]
        targetName, targetLen, targetS, targetE = words[5:5+4]
        targetDirection = words[9] 
 
        querySeq = words[-3]   # ATCG-
        align = words[-2]      # |,*
        targetSeq = words[-1]  # ATCG-    
 
    f.close()
    query = sequence.Sequence(queryName, queryLen, queryS, queryE, querySeq, queryDirection) 
    target = sequence.Sequence(targetName, targetLen, targetS, targetE, targetSeq, targetDirection)
    alignObj = alignment.Alignment(query, target, align)  # blasr first reference, second contig

    #alignObj = alignment.Alignment(target, query, align)  # blasr first contig, second reference
    return alignObj

'''
def read_columns(filename):
    f = open(filename, "r")
    columns = {}
    for line in f:
        words = line.strip().split() 
        if line.startswith("reference"):
            #print (words[1])
            pos = int(words[2]) 
            assert pos not in columns
            columns[pos] = column.Column(pos)
            columns[pos]._cov = int(words[4])  
        elif (len(words) > 0 ):
            #columns[pos]._print()
            columns[pos]._diff_cov.append(int(words[1]))
    for pos in columns:
        columns[pos]._diff_cov.sort() 
        columns[pos]._diff_cov.reverse() 
    #columns[0]._print()
    #columns[2925]._print() 
    return columns       
#read_columns(sys.argv[1]) 
#read_align(sys.argv[1])
'''
