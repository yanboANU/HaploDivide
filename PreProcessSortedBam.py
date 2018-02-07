#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import contig


#input : sorted_bam contigs(deal with mutiple alignment)
#process: remove reads who align mutiple(only keep one position for one read), move gap to right 
#output: new_bam


def remove_Duplicate(samfile):
   
    reads = {}
    count = 0
    for read in samfile.fetch():
        count += 1
        readName = read.query_name
        readMapScore = read.mapping_quality

        if (readName not in reads) or (reads[readName].mapping_quality < readMapScore):
            reads[readName] = read    
    print ("Original reads number: ",count)
    print ("After remove duplicate read number:", len(reads))
    return reads

def calc_Query_Length(cigar):
    ##########################################
    #according to cigar to calclate query length
    #cigar only have 4 (1,2,3,4) symbols
    ########################################## 
    query_length = 0
    
    for (i,j) in cigar:
        if i != 2:
            query_length += j
    return query_length 

def align2Cigar(readAlignSeq, referAlignSeq, orignalCigar):
    ####################################    
    # AATTC-GG   ==> (0,5) (3, 1) (0,2)
    # AATTCCGG
    # alignment to cigar
    ###################################  
    cigar = []
    i = 0 
    end = len(readAlignSeq)
    if orignalCigar[0][0] == 4:
        cigar.append(orignalCigar[0])
        i = orignalCigar[0][1] 
    if orignalCigar[-1][0] == 4:  
        end = len(readAlignSeq) - orignalCigar[-1][1]
    # maybe dead cycle
    while i < end:
        count = 0
        while i < end and readAlignSeq[i] != '-' and referAlignSeq[i] == '-': # refer gap
            count += 1
            i = i + 1
        if count > 0:
            cigar.append((1, count))
        if i >= end:
            break

        count = 0
        while (i < end and readAlignSeq[i] !='-' and referAlignSeq[i] != '-'):
            count += 1
            i = i+1
        if count > 0:
            cigar.append((0, count)) 
        if i >= end:
            break

        count = 0 
        while i < end and readAlignSeq[i] == '-' and referAlignSeq[i] != '-': # read gap
            count += 1
            i = i + 1
        if count > 0:
            cigar.append((2, count))
        
        #i = i+1
        if i >= end:
            break
    count = 0
    print (i) 
    if i < len(readAlignSeq):
        for c in readAlignSeq[i:]:
            if c != '-':
                count +=1
         
        cigar.append((4,count))
    '''    
    print (count)
    sum1 = 0
    for (i,j) in orignalCigar:
        sum1 += j
    
    sum2 = 0  
    for (i,j) in cigar:
        sum2 += j
    #assert sum1 == sum2 # maybe not equal, differ align way differ insert number, align length(map + insert) not always same
    #assert len(readAlignSeq) == sum2 
    '''
    return cigar     

def alignPos_2_Sequence(readAlignPos, readSeq):  # include gap
# from align pos
    seq = ""
    for c in readAlignPos:
        if type(c) == int:
            seq = seq + readSeq[c] 
        else:
            seq = seq + '-'
    return seq
     
def print_Alignment(readAlignPos,referAlignPos,readSeq,contigSeq):

    
    print (alignPos_2_Sequence(readAlignPos, readSeq)) 
    print (alignPos_2_Sequence(referAlignPos, contigSeq))
    #print (referAlignPos)

def write_No_Repeat(reads, samfile):

    newBam = pysam.AlignmentFile("no_repeat_reads.bam", "wb", template=samfile) 
    for (name, read) in reads.items():
        print ("read len:", read.qlen) # not read length
        print ("read len:", len(read.query_sequence))   
        print ("align len:", len(read.get_aligned_pairs(True,False)))
        print ("mapping quality", read.mapping_quality)
        '''
        newRead = pysam.AlignedSegment()
        newRead.query_name = read.query_name
        newRead.query_sequence=read.query_sequence
        newRead.flag = read.flag
        newRead.reference_id = read.reference_id
        newRead.reference_start = read.reference_start
        newRead.mapping_quality = read.mapping_quality
        newRead.cigar = read.cigar
        newRead.next_reference_id = read.next_reference_id
        newRead.next_reference_start = read.next_reference_start
        newRead.template_length= read.template_length
        newRead.query_qualities = read.query_qualities
        newRead.tags = read.tags
        '''
        newBam.write(read)
        #break  
    newBam.close()

def write_New_Bam(readAndNewAlign, samfile):

    newBam = pysam.AlignmentFile("new.bam", "wb", template=samfile) 
    for (read, newCigar) in readAndNewAlign:
        newRead = pysam.AlignedSegment()
        newRead.query_name = read.query_name
        newRead.query_sequence=read.query_sequence
        newRead.flag = read.flag
        newRead.reference_id = read.reference_id
        newRead.reference_start = read.reference_start
        newRead.mapping_quality = read.mapping_quality
        newRead.cigar = newCigar
        newRead.next_reference_id = read.next_reference_id
        newRead.next_reference_start = read.next_reference_start
        newRead.template_length= read.template_length
        newRead.query_qualities = read.query_qualities
        newRead.tags = read.tags
        newBam.write(newRead)
        #break  
    newBam.close()

def realignment(reads, contigs):
  
    #newAlignedPairs = [] 
    readAndNewAlign = [] 
    contigSeq = contigs[0].Seq
    #readList = ["S1_1592","S1_1926","S1_1181","S1_1778"]
    for (readName, read) in reads.items():
        '''        
        if readName not in set(readList):
            continue
        '''
        alignedPairs = read.get_aligned_pairs()
        readSeq = read.query_sequence
        #print (alignedPairs)
        readSeqLen = len(readSeq)
        print ("query length:",len(readSeq))
        print ("old cigar query length:",calc_Query_Length(read.cigar))
        #print (read.cigarstring)
        print (read.cigar) 
        #print (read.query_alignment_sequence)
        #print (read.get_aligned_pairs(False, True))
        
        alignLen = len(alignedPairs)
         
        readAlignPos = []
        referAlignPos = []
        for (readPos, referPos) in alignedPairs:
            readAlignPos.append(readPos)
            referAlignPos.append(referPos)  
        
        print_Alignment(readAlignPos,referAlignPos,readSeq,contigSeq) 
      
        for i in range(alignLen-1):
            # pushing refer gap
            # i-1    i
            # 5000  none
            # 10     11 
            if (not type(referAlignPos[i]) == int) and (type(referAlignPos[i-1]) == int) and (type(readAlignPos[i]) == int):
                j = i + 1
                while (j  < alignLen):
                    c = referAlignPos[j]
                    if (type(c) == int):
                        if (contigSeq[c] == readSeq[readAlignPos[i]]):
                            referAlignPos[i] = c
                            referAlignPos[j] = None
                            #print ''.join(qstr)     
                            #print ''.join(tstr)
                        break
                    j += 1
            # pushing query gap
            if (not type(readAlignPos[i]) == int) and (type(readAlignPos[i-1]) == int) and (type(referAlignPos[i]) == int):
                j = i + 1
                while (j  < alignLen):
                    c = readAlignPos[j]
                    if (type(c) == int):
                        if not type(referAlignPos[i]) == int:
                            print ("position %s: (%s, %s)"% (i, readAlignPos[i], referAlignPos[i]))

                        assert type(referAlignPos[i]) == int
                        if (readSeq[c] == contigSeq[referAlignPos[i]]):
                            readAlignPos[i] = c
                            readAlignPos[j] = None
                            #print ''.join(qstr)     
                            #print ''.join(tstr)
                        break
                    j += 1
            # 10    None
            # 5000  None
            if ( (not type(referAlignPos[i]) == int) and (type(referAlignPos[i-1]) == int) 
                 and (not type(readAlignPos[i]) == int) and (type(readAlignPos[i-1]) == int) ):
                j = i + 1
                t = i + 1 
                while (j  < alignLen):
                    c = readAlignPos[j]
                    if (type(c) == int):
                        break
                    j += 1
                while (t  < alignLen):
                    c = referAlignPos[t]
                    if (type(c) == int):
                        break
                    t += 1
                if (t < alignLen) and (j < alignLen) and (contigSeq(referAlignPos[t]) == readSeq(readAlignPos[j])):
                        readAlignPos[i] = readAlignPos[j]
                        readAlignPos[j] = None
                        referAlignPos[i] = referAligPos[t]
                        referAlignPos[t] = None
 
        #print (readAlignPos)
        #print (referAlignPos) 

        print_Alignment(readAlignPos,referAlignPos,readSeq,contigSeq) 
        readAlignSeq = alignPos_2_Sequence(readAlignPos, readSeq)
        referAlignSeq = alignPos_2_Sequence(referAlignPos, contigSeq)   
        newCigar = align2Cigar(readAlignSeq, referAlignSeq, read.cigar)
        readAndNewAlign.append((read, newCigar))
        print ("new cigar", newCigar)
        lenByNewCigar =  calc_Query_Length(newCigar) 
        print ("new cigar, query length:",calc_Query_Length(newCigar))
        assert lenByNewCigar == readSeqLen 
        #break
    return readAndNewAlign   
if __name__ == "__main__":

    # sorted bam
    samfile = pysam.AlignmentFile(sys.argv[1], "rb") 
    # contig
    #contigs = contig.read_Contig(sys.argv[2])     

    # step one
    reads = remove_Duplicate(samfile)
    write_No_Repeat(reads,samfile) 
    # step two: realignment, merge gap togther
    #readAndNewAlign = realignment(reads, contigs)
    #write_New_Bam(readAndNewAlign, samfile)
   
    samfile.close() 

 
