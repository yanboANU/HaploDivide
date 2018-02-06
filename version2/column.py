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
    
    def __init__(self, refPos):
        self._ref_pos = refPos

        self._map_content = {}
        self._nucleotide = ''
        self._cov = 0
        self._insert_content = {}
        self._is_insert = -1
        self._is_mutation = -1
        self._is_delete = -1
        self._is_stable = -1    
        #self.best_Content =
        #self.second_Content =
    '''
    def __init__(self, mapContent, nucleotide, refPos, cov):

        self._map_content = mapContent
        self._nucleotide = nucleotide
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
    '''   
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
        
        _is_delete = 0
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
'''
def init_Columns_2(bamfile, columns):
    
    #forColumns = {}
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

# input .bam file 
def init_Columns(bamfile, contig):
    columns = {}

    for p in bamfile.pileup(contig._name): 

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
        nucleotide = contig._seq[p.pos]  

        column = Column(mapContent, nucleotide, p.pos, p.n)
        columns[p.pos] = column
    
    write_Columns(contig._name + "_columns1" ,columns)    
    #return init_Columns_2(bamfile, columns)   
    return columns   
'''

def init_Columns_3(bamfile, contig):
    
    columns = {}
    for read in bamfile.fetch():
        alignedPairs = read.get_aligned_pairs()
        readName = read.query_name
        print (readName, len(alignedPairs))
        print (alignedPairs)
        i = 0
        label = 1
        while i < len(alignedPairs):
            (readPos, referPos) = alignedPairs[i]
            '''
            if readName == "read_1730000":
                print (readPos, referPos)
            '''
            # "get reference first not None"
            while label and isinstance(readPos, int) and not type(referPos) == int:
                i += 1
                (readPos, referPos) = alignedPairs[i]
                print (readPos, referPos)
                continue
            label = 0

            if type(referPos) == int:
                if referPos not in columns:
                    columns[referPos] = Column(referPos)
                columns[referPos]._cov += 1
                columns[referPos]._nucleotide = contig._seq[referPos]
                if type(readPos) == int:
                    c = read.query_sequence[readPos] 
                    if c not in columns[referPos]._map_content:
                        columns[referPos]._map_content[c] = []
                    columns[referPos]._map_content[c].append(readName)
                elif not type(readPos) == int:
                    c = '*'
                    if c not in columns[referPos]._map_content:
                        columns[referPos]._map_content[c] = []
                    columns[referPos]._map_content[c].append(readName)
                i += 1
            '''    
            if readName == "read_1730000":
                sys.exit()
            ''' 
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
                if referPos_pre not in columns:
                    columns[referPos_pre] = Column(referPos_pre)
                if insertSeq not in columns[referPos_pre]._insert_content:
                    columns[referPos_pre]._insert_content[insertSeq] = []

                columns[referPos_pre]._insert_content[insertSeq].append(readName)
            #i = i+1    
    write_Columns(contig._name + "_columns3", columns)        
    return columns





def write_Columns(filename, columns):
    fout = open(filename, "w")

    for (refPos, column) in columns.items():
        fout.write("reference position: %s  coverage: %s reference nucleotide: %s \n" % (refPos, column._cov, column._nucleotide))
        tools.write_Map(fout, column._map_content)

