#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import time
import contig
import tools
import column

class Phasing:

    def __init__(self, columns,  mutation, delete, insert, stable, obLen = 3):
        self_colums = columns
        
        self._stableRange = tools.pos_2_Range(stable)    
        
        self._snp_mutation = _get_SNP(mutation, obLen)
        
        # need filter
        self._snp_delete = delete(delete, obLen)
        self._snp_insert = insert(insert, obLen)

        # temporary 
        self._snp = self._snp_mutation 

    def _get_SNP(self, pos, obLen):
        # unfinish 
        snp = [] 
        for a in pos:
            if tools.is_SubRange(a-obLen, a-1,_stableRange) and tools.is_SubRange(a+1,a+obLen, _stableRange):
                    snp.append(a)

        return snp

    def _label_reads(self):
        
        readsLabel = {} 
	for i in range(len(_snp)):
	    p = _snp[i] 
	    (nuc, supportReadName) = _columns[p]._map_content
		supports = [] #  
		sortedNuc = sorted_Map_Value(nuc)
		a = sortedNuc[0][0]
		b = sortedNuc[1][0]
		for (k, v) in supportReadName.items():
		    #print (k,v)
		    if k != a and k != b:
			supports.extend(v)
		#print (supports)
		#print ("bbb")
		#sys.exit("stop")
		for readId in supportReadName[a]:
		    if readId not in readsFlag:
			readsFlag[readId] = [3]*len(ObSNP)
		    readsFlag[readId][i] = 0
	 
		for readId in supportReadName[b]:
		    if readId not in readsFlag:
			readsFlag[readId] = [3]*len(ObSNP)
		    readsFlag[readId][i] = 1

		for readId in supports:
		    if readId not in readsFlag:
			readsFlag[readId] = [3]*len(ObSNP)
		    readsFlag[readId][i] = 2
	    sortFlag = sorted_Map_Value(readsFlag, False)
	    print ("show")
	    for ele in sortFlag:
		print (ele)
 
    return readsFlag       

def phasing(usefulPhases, phase0, phase1, label0, label1, readsSupport):
     
    print ("in phasing")
    print (label0 , label1)
    print (usefulPhases[0][0], usefulPhases[1][0])
    if len(label0) == 0: 
        label0 = usefulPhases[0][0]
        label1 = usefulPhases[1][0] 
        #phase0.extend(readsSupport[label0])
        #phase1.extend(readsSupport[label1])
    else:
        assert ( label0[-2:] == usefulPhases[0][0][:2] or  
                 label0[-2:] == usefulPhases[1][0][:2] )
        if label0[-2:] == usefulPhases[0][0][:2]: 
            label0 = label0 + usefulPhases[0][0][-1] 

        
   

