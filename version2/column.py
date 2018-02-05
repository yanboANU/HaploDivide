#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import time
import contig
import tools


class Column:
    # nucleotide: in reference
    # mapContent is a map 'A': support read list
    def __init__(self, mapContent, nucleotide, refPos, cov):

        self._map_content = mapContent
        self._nucleptide = nucleotide
        self._ref_pos = refPos
        self._cov = cov
        self._insert_content = {}
        self._is_insert = -1
        self._is_mutation = -1
        self._is_delete = -1
        self._is_stable = -1    
        #self.best_Content =
        #self.second_Content =
        #self._is_snp_mutation  
       
    def _is_Insert():
        
        insertNum = 0
        insertLength = 0 
        for (content, reads) in _insert_content.items():
            insertNum += len(reads)
            insertLength += len(content)*len(reads)
    
        # waiting statics data 
        if insertNum > _cov*0.2 or insertLength > _cov*0.5:
            return True
        return False

    #mutation
    #stable
    #delete
    def _set_Mutation_or_Delete(self):
        
e       _is_delete = 0
        _is_mutation = 0
        print (self._ref_pos)
        print (self._map_content)
        s = tools.sorted_Map_Value_Len(self._map_content)
        _map_content = s   
        print (s)

        if len(s[0][1]) >= _cov*0.7 and s[0][0] != _nucleotide:
            sys.exit("0.7 cov different from ref")
   
        if len(s[0][1]) < _cov*0.7 and len(s[1][1]) >= _cov*0.7:
            if s[0][0] == '*' or s[1][0] == '*':
                _is_delete = 1
            else:
                _is_mutation = 1 

    
    def _set_Lable(self):
        
        if _is_Insert():
            _is_insert = 1  
        else:
            _is_insert = 0 
       
         _set_Mutation_or_Delete()
        

        if _is_insert ==0 and _is_mutation == 0 and _is_delete == 0:
            is_stable = 1
        else:
            is_stable = 0  


# input .bam file 
def init_Columns(bamfile, contigSeq):
    columns = {}

    for p in bamfile.pileup("contig1001_circular:0_len:172095_cov:39",2000,20100): 

        mapContent = {}
        for pread in p.pileups:
            #print (pread.query_position)
            if pread.is_del:
                cc = '*'
            else:
                cc = pread.alignment.query_sequence[pread.query_position]
            if cc.upper() not in mapContent:
                mapContent[cc.upper()] = []
            mapContent[cc.upper()].append(pread.alignment.query_name)  
        nucleotide = contigSeq[p.pos]  

        column = Column(mapContent, nucleotide, p.pos, p.n)
        columns[p.pos] = column
    return init_Columns_2(bamfile, columns)   
    #return columns   

def init_Columus_2(bamfile, columns):
        
    for read in bamfile.fetch():
        alignedPairs = read.get_aligned_pairs()
        readName = read.query_name  
        
        i = 0
        while i < len(alignedPairs):
            (readPos, referPos) = alignedPairs[i]
            if isinstance(readPos, int) and not type(referPos) == int:
                (readPos_pre, referPos_pre) = alignedPairs[i-1]
                
                if not type(referPos_pre) == int:
                    i = i+1
                    continue  
                if referPos_pre >= read.reference_end-1:
                    break
                
                #columns[referPos_pre]._insert_content = [] 
                insertSeq = read.query_sequence[readPos] 
                while i+1 < len(alignedPairs): 
                    i = i+1
                    (readPos, referPos) = alignedPairs[i]
                    if isinstance(readPos, int) and not type(referPos) == int:
                        insertSeq = insertSeq + read.query_sequence[readPos] 
                    else:
                        break
             
                if insertSeq not in columns[referPos_pre]._insert_content:
                    columns[referPos_pre]._insert_content[insertSeq] = []

                columns[referPos_pre]._insert_content[insertSeq].append(readName)
    return columns
