#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import contig
import tools

#input : 0_16.bam
#every read need be covered 0.5*read_length 


def get_some_position_pattern(pp, ref):
    pattern = []
    for (c, num) in pp.items():
        s =  ref[c-9:c+10]
        if s.count('N') == 0:
            pattern.append((s, num))
    return pattern        

def calc_same_len(c):
    sameLen = 1
    mid = int(len(c)/2)
    for i in range(mid+1, len(c)):
        if c[i] == c[mid]:
            sameLen += 1
        else:
            break
    for i in range(mid-1, -1, -1):
        if c[i] == c[mid]:
            sameLen += 1
        else:
            break
    return sameLen
    


def stat_pattern(pattern):
    
    xdy = []
    xddy = []
    xdddy = []
    xddddy = []
    xdddddy = []

    xdyNum = 0
    xddyNum = 0
    xdddyNum = 0
    xddddyNum = 0
    xdddddyNum = 0
    patternLen = float(len(pattern))
    mid = int(len(pattern[0][0])/2)
    print (mid)
    totalDelete = 0
    for (c, num) in pattern: # middle is delete position
        #assert not tools.same_aA(c[mid-1], c[mid])       
        totalDelete += num
        sameLen = calc_same_len(c)
        
        if sameLen >= 5: 
            xdddddy.append(c)
            xdddddyNum += num
        elif sameLen == 4:
            xddddy.append(c)
            xddddyNum += num
        elif sameLen == 3:
            xdddy.append(c)
            xdddyNum += num
        elif sameLen == 2:
            xddyNum += num
            xddy.append(c)
        elif sameLen == 1: 
            xdy.append(c)
            xdyNum += num
        else:
            print ("impossible")
           
    print ("len xdy", xdyNum, xdyNum/totalDelete)
    print (xdy[0:10])
    print ("len xddy", xddyNum, xddyNum/totalDelete)
    print (xddy[0:10])
    print ("len xdddy", xdddyNum, xdddyNum/totalDelete)
    print (xdddy[0:10])
    print ("len xddddy", xddddyNum, xddddyNum/totalDelete)
    print (xddddy[0:10])
    print ("len xdddddy", xdddddyNum, xdddddyNum/totalDelete)
    print (xdddddy[0:10])


def calc_delete_rate_in_homo_reigion(delePos, ref):
    pattern = get_some_position_pattern(delePos, ref)
    pattern = get_all_position_pattern(ref)
    stat_pattern(pattern)



def get_ref_pattern(ref):
    if ref[-1] != 'N':
        ref = ref + "N" 
    ref = ref.upper()
    refLen = len(ref)
    lenCnt = {} #key: length value: count
    tempLen = 1
    i=0
    while i<refLen and ref[i] == "N":
        i += 1
    pre = ref[i]
    i += 1
    # must end with N, otherwise count a error
    while i<refLen:
        if ref[i] == pre:
            tempLen += 1
            i += 1
        else:
            if tempLen not in lenCnt:
                lenCnt[tempLen] = 1
            else:
                lenCnt[tempLen] += 1
            tempLen = 1
            while i<refLen and ref[i] == "N":
                i += 1
            if i<refLen:
                pre = ref[i]
                i += 1   
    length = 0
    totalCnt = 0
    for key in lenCnt:
        length += key*lenCnt[key]
        totalCnt += lenCnt[key] 
    for key in lenCnt:    
        print (key, lenCnt[key], float(lenCnt[key])/totalCnt)
    print (length)
    print (ref.count('N'))
    print (refLen)
    assert length + ref.count("N") == refLen 


if __name__ == "__main__":

    #input: *bam contig


    if len(sys.argv) < 2:
        print ("python3 " + sys.argv[0] + " *bam contig")
        sys.exit()
   
    print ("python3 " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2])

    samfile = pysam.AlignmentFile(sys.argv[1], "rb")
    count = 0
    contigs = contig.read_Contig(sys.argv[2])
    contigName, contig = contigs.popitem()
    length = 0 
    inss = 0
    muts =  0
    deles = 0
    insertLen = 0
    delePos = {}
    insePos = []
    seq = ""
    for read in samfile.fetch("1"):

        count += 1
        if count % 100 != 1:
            continue
        if count / 100 >10000:
            break
        readMapScore = read.mapping_quality
        mapped = float(read.qlen)/len(read.query_sequence)
        i = 0
        ins = 0
        mut = 0
        dele = 0 
        alignedPairs = read.get_aligned_pairs()
        lenAlignedPairs = len(alignedPairs)
        #print (alignedPairs)
        #print (lenAlignedPairs)
        #print (read.query_sequence)
        seq = seq + read.query_sequence
        seq = seq + "N"
        '''
        if count >= 10:
            break
        '''    
            #sys.exit()

    get_ref_pattern(seq)
    '''
        length += lenAlignedPairs
        label = 1
        while i+1 < lenAlignedPairs:
            (readPos, referPos) = alignedPairs[i] #check referPos always increse or not
            while label and isinstance(readPos, int) and not type(referPos) == int:
                i += 1
                (readPos, referPos) = alignedPairs[i]  
            label = 0
            # get columns(align info)
            if type(referPos) == int:
                if type(readPos) == int:
                    r = read.query_sequence[readPos]
                    c  = contig.Seq[referPos]
                    if r != c and r.upper() !=c and r != c.upper():
                        mut += 1    
                elif not type(readPos) == int: #delete
                    dele += 1
                    if referPos not in delePos:
                        delePos[referPos] = 0
                    delePos[referPos] += 1
                i += 1    
            inslen = 0
            if isinstance(readPos, int) and not type(referPos) == int:
                ins += 1
                inslen += 1
                while i+1 < lenAlignedPairs: 
                    i += 1
                    (readPos, referPos) = alignedPairs[i]
                    if isinstance(readPos, int) and not type(referPos) == int:
                        inslen += 1 
                    else:
                        insertLen += inslen
                        break
                       
        #print (lenAlignedPairs, mut, ins, dele, mapped, readMapScore, insertLen)
        
        inss += ins
        muts +=  mut
        deles += dele

    print ("reads number:", count)
    print ("total align length:", length)   
    print ("mutate number:", muts, muts/float(length))   
    print ("delete number:", deles, deles/float(length))
    print ("insert number:", inss, inss/float(length))  
    print ("total insert length:", insertLen/float(length))   
    samfile.close()
    print (delePos)
    calc_delete_rate_in_homo_reigion(delePos, contig.Seq.upper()) # delete homo based on ref, not read
    '''
