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
        #self.best_Content =
        #self.second_Content = 
       
  

    #mutation
    #stable
    #delete
    def _is_MDS(self):
        print (self._ref_pos)
        print (self._map_content)
        s = tools.sorted_Map_Value_Len(self._map_content)
        print (s)



# input .bam file 
def init_Columns(bamfile, contigSeq):
    columns = []

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
        columns.append(column)
    return columns   

def init_Columus_2(bamfile):
        

