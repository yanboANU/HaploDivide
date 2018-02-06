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
 
    def pre_init(self):

        mutation = []
        insert = []
        delete = []
        stable = []  
        for (refPos, c) in _columns.items():
            c._set_Lable()
            #sys.exit()
            if c._is_insert == 1:
                insert.append(refPos)
        
            if c._is_mutation == 1:
                mutation.append(refPos)
        
            if c._is_delete == 1:
                delete.append(refPos)
        
            if c._is_stable == 1:
                stable.append(refPos)
        return stable, mutation, delete, insert

    def __init__(self, columns, obLen = 3):
        self._colums = columns

        stable, mutation, delete, insert = pre_init() 
        
        self._stableRange = tools.pos_2_Range(stable)    
        
        self._snp_mutation = _get_SNP(mutation, obLen)
        
        # need filter
        self._snp_delete = delete(delete, obLen)
        self._snp_insert = insert(insert, obLen)


        # output
        # temporary 
        self._snp = self._snp_mutation
        self._coverage = []
        self._label0 = ""
        self._label1 = ""
        self.phase1 = []
        self.phase2 = [] 
           

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
	    content = _columns[p]._map_content
	    #a = content[0][0]
	    #b = content[1][0]
	    for readId in content[0][1]:
	        if readId not in readsLabel:
		    readsLabel[readId] = [3]*len(_snp)
		readsLabel[readId][i] = 0
	 
	    for readId in content[1][1]:
		if readId not in readsLabel:
	            readsLabel[readId] = [3]*len(_snp)
		readsLabel[readId][i] = 1
            j = 2  
            while j < len(content): 
         	for readId in content[j][1]:
		    if readId not in readsLabel:
		        readsLabel[readId] = [3]*len(_snp)
		    readsLabel[readId][i] = 2
	'''
        sortFlag = sorted_Map_Value(readsLabel, False)
	print ("show")
	for ele in sortFlag:
        print (ele)
        ''' 
    return readsLabel       

    def _phasing(window=3):
            
        readsLabel = _label_reads()
        _coverage = (len(_snp) - window)*[0]   
        for i in range(len(_snp)-window):
            phases = {}
            for (read, label) in readsLabel.items():

                f = ''.join( str(j) for j in label[ i :i+window] )
                    if f not in phases:
                        phases[f] = []
                phases[f].append(read)
        #print (phases)   
        
            print (_snp[i:i+3])
            sortedPhases = tools.sorted_Map_Value_Len(phases)
            cov = 0
            allCoverPhases = []
            for (label, reads) in sortedPhases:
                if label.find('3') == -1:
                    cov += len(reads)
                    allCoverPhases.append((label, len(reads)))
            _coverage.append(cov)  

            if len(allCoverPhases) == 1:
    

            elif ( len(allCoverPhases)>=2 and allCoverPhases[0][1] > cov * 0.2 and allCoverPhases[1][1] > cov * 0.2
               and tools.is_Bool_Reverse(allCoverPhases[0][0], allCoverPhases[1][0]) ):
                _phasing_one_window(allCoverPhases, phases)
                #sys.exit()
            else:
                print ("not statified")
                #reclass_Intersection(phase0, phase1, label0, label1, readsFlag)
                #fout.write("%s %s\n" % (label0, phase0))
                #fout.write("%s %s\n" % (label1, phase1))
                #fout.write("\n")
                #sys.exit()


           
    def _phasing_one_window(allCoverPhases, labelReads):
     
        print ("in phasing")
        print (label0 , label1)
        print (allCoverPhases[0][0], allCoverPhases[1][0])
        if len(label0) == 0: 
            _label0 = allCoverPhases[0][0]
            _label1 = allCoverPhases[1][0] 
        else:
            assert ( _label0[-2:] == allCoverPhases[0][0][:2] or  
                   _label0[-2:] == allCoverPhases[1][0][:2] )
            if _label0[-2:] == allCoverPhases[0][0][:2]: 
                _label0 = _label0 + allCoverPhases[0][0][-1] 
                _label1 = _label1 + allCoverPhases[1][0][-1] 
            elif _label0[-2:] == allCoverPhases[1][0][:2]:     
                _label0 = _label0 + allCoverPhases[1][0][-1]
                _label1 = _label1 + allCoverPhases[0][0][-1]
        else:
            print ("error type1")

   
        for v in allCoverPhases:
            if tools.hamming_Distance(v[0], _label0[-3:]) < tools.hamming_Distance(v[0], _label1[-3:]):
                _phase0.extend(labelReads[v[0]])
 
            elif tools.hamming_Distance(v[0], _label0[-3:]) > tools.hamming_Distance(v[0], _label1[-3:]):
                _phase1.extend(labelReads[v[0]])
            else:
                print ("same distance", v)
            

        print (_label0, len(_phase0), _phase0)      
        print (_label1, len(_phase1), _phase1)      
        
   
        print ("intersection:", set(_phase0).intersection(set(_phase1)) )
