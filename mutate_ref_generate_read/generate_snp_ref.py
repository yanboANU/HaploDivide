#!/usr/bin/python


import os
import sys
import random
import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

def read_positions(filename, start, end):
    f = open(filename, "r")
    pos = {}
    for line in f:
        words = line.split()
        p = int(words[0])
        if p > end:
            break
        if p>= start and p <= end:
           pos[p] = (words[1], words[2])
    f.close()
    print "mutate length",len(pos)
    return pos, sorted(pos.keys());


def get_new_record(record, start, end, pos, posKeys, typeSNP):
   
    s = str(record.seq)
    lenS = len(s)
    #assert lenS == end-start+1
    newSeq = ""
    seqS = 0
    #if typeSNP == 0 or typeSNP == 1:   # snp, insert1s   
    insertNum = 0
    snpNum = 0
    deleteNum = 0
    foutI = open("Insert_record","w")
    foutD = open("Delete_record","w")
    foutS = open("SNP_record","w")
    for key in posKeys:
        if len(pos[key][0]) == 1 and len(pos[key][1]) <=2: # means insert or snp
            if key - start <= seqS:
                continue
            newSeq = newSeq + s[seqS : key-start]
            assert s[key-start] == pos[key][0] or s[key-start].upper() == pos[key][0] 
            seqS = (key-start)+1
            newSeq = newSeq + pos[key][1]
            if len(pos[key][1]) == 2:
                insertNum += 1
                foutI.write("%s %s %s\n" % (key, pos[key][0], pos[key][1]))
            else:
                assert(len(pos[key][1]) == 1)
                snpNum += 1
                foutS.write("%s %s %s\n" % (key, pos[key][0], pos[key][1]))
            continue

        if len(pos[key][0]) == 2: # delete
            if key - start + 1 <= seqS:
                continue
            #print key, pos[key]
            foutD.write("%s %s %s\n" % (key, pos[key][0], pos[key][1]))
            newSeq = newSeq + s[seqS : key-start+1]
            assert s[key-start] == pos[key][1] or s[key-start].upper() == pos[key][1] 
            seqS = (key-start)+2
            deleteNum += 1
    newSeq = newSeq + s[seqS: lenS] 
    foutI.close()
    foutD.close()
    foutS.close()
        #print len(newSeq)
    '''
    #if typeSNP == 2:   # delete 1
        #for key in posKeys:
            #newSeq = newSeq + s[seqS : key-start+1]
            #print key, pos[key]
            #print s[key-start]

            if len(pos[key][0]) == 2: # delete 
            assert s[key-start] == pos[key][1] or s[key-start].upper() == pos[key][1] 
            seqS = (key-start)+2
    newSeq = newSeq + s[seqS: lenS]    
    '''
    print "old sequence length:", len(s)
    print "new sequence length:", len(newSeq)
    print "insert num:", insertNum
    print "delete num:", deleteNum
    print "SNP num:", snpNum
    assert len(s) + insertNum -deleteNum == len(newSeq)
    newRecord = SeqRecord(Seq(newSeq),id=record.id+"_1in1kSNP_1in10kIndel")
    return newRecord


if __name__ == "__main__":

    #
    if(len(sys.argv) < 3):
        print "Usage: python generate_snp_ref simulate_snp_position_filename, typeSNP, ref.fasta(block1_*.fasta)"
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3])
    typeSNP = int(sys.argv[2])
    for i in range(len(sys.argv)-3):
        #print "file",i, i+2
        fileName = sys.argv[i+3]
        words = fileName.split('_')
        if len(words) == 1:
            start = 1
            end = sys.maxint
            words[0] = "chr1_baesd_NIST"
        else:
            start = int(words[1])    # position, not index, position = index + 1
            end =   int(words[2].split('.')[0])    # position

      
        pos, posKeys = read_positions(sys.argv[1], start, end)
        record = SeqIO.read(open(fileName), "fasta")
        newRecord = get_new_record(record, start, end, pos, posKeys, typeSNP)
        twoRecord = []
        twoRecord.append(record)
        twoRecord.append(newRecord)
        SeqIO.write(twoRecord, words[0]+ ".fasta",'fasta')
    

