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
    
    def _is_Insert(self):
        
        insertNum = 0
        insertLength = 0 
        for (content, reads) in self._insert_content.items():
            insertNum += len(reads)
            insertLength += len(content)*len(reads)
    
        # waiting statics data 
        if insertNum > self._cov*0.2 or insertLength > self._cov*0.5:
            #print ("insert")
            return True
        return False

    #mutation
    #stable
    #delete
    def _set_Mutation_or_Delete(self):
        
        self._is_delete = 0
        self._is_mutation = 0
        #print (self._ref_pos)
        #print (self._map_content)
        s = tools.sorted_Map_Value_Len(self._map_content)
        self._map_content = s   
        #print (s)

        if self._cov <= 10:
            return       

        if len(s[0][1]) >= self._cov*0.7 and s[0][0] != self._nucleotide:
            print ("reference position %s wrong" % (self._ref_pos))
            #sys.exit("0.7 cov different from ref")
        
        if len(s[0][1]) < self._cov*0.7 and len(s[1][1]) >= self._cov*0.3:
            
            #print ("enter 1")
            if s[0][0] == '*' or s[1][0] == '*':

                #print ("enter 2")
                self._is_delete = 1
            else:
                self._is_mutation = 1 
 
    def _set_Lable(self):
        
        if self._is_Insert():
            self._is_insert = 1  
        else:
            self._is_insert = 0 
       
        self._set_Mutation_or_Delete()
        

        if self._is_insert ==0 and self._is_mutation == 0 and self._is_delete == 0:
            self._is_stable = 1
        else:
            self._is_stable = 0 
        #print (self._is_insert, self._is_mutation, self._is_delete, self._is_stable)    

def init_Columns(bamfile, contig):
    
    columns = {}
    for read in bamfile.fetch(contig._name):
        alignedPairs = read.get_aligned_pairs()
        readName = read.query_name
        #print (readName, len(alignedPairs))
        #print (alignedPairs)
        i = 0
        label = 1
        while i < len(alignedPairs):
            (readPos, referPos) = alignedPairs[i]
            
            # "get reference first not None"
            while label and isinstance(readPos, int) and not type(referPos) == int:
                i += 1
                (readPos, referPos) = alignedPairs[i]
                #print (readPos, referPos)
                #continue
            label = 0

            # get columns(align info)
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
            
            #get insert info
            if isinstance(readPos, int) and not type(referPos) == int:
                #print ("%s get insert" % (i))
                (readPos_pre, referPos_pre) = alignedPairs[i-1]    
                assert type(referPos_pre) == int 
                '''
                if not type(referPos_pre) == int:
                    i = i+1
                    continue  
                '''    
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
                #print ("pos: %s, content: %s" % (referPos_pre, insertSeq))
                if insertSeq not in columns[referPos_pre]._insert_content:
                    columns[referPos_pre]._insert_content[insertSeq] = []

                columns[referPos_pre]._insert_content[insertSeq].append(readName)
        #print ("check", columns[13]._insert_content)        
        #sys.exit()    
        #i = i+1    
    write_Columns(contig._name , columns)        
    return columns





def write_Columns(contigName, columns):
    filename = contigName + "_columns"
    fout1 = open(filename, "w")
    for (refPos, column) in columns.items():
        fout1.write("reference position: %s  coverage: %s reference nucleotide: %s \n" % (refPos, column._cov, column._nucleotide))
        tools.write_Map(fout1, column._map_content)
    fout1.close()

    filename = contigName + "_Ins"
    fout2 = open(filename, "w")
    for (refPos, column) in columns.items():
        if len(column._insert_content) != 0:
            fout2.write("reference position: %s  coverage: %s reference nucleotide: %s \n" % (refPos, column._cov, column._nucleotide))
            tools.write_Map(fout2, column._insert_content)
    fout2.close()


