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

    def _write_result(self):
        fout = open(self._contig._name+"_phasing_result","w")
        n =  len(self._label0s)
        fout.write("name: %s len: %s\n" % (self._contig._name, self._contig._len))
        fout.write("snp mutation number: %s\n" % (len(self._snp_mutation)))
        fout.write("snp delete number: %s\n" % (len(self._snp_delete)))
        fout.write("divide in %s sequences\n" % (n)) 
        #assert n == len(self._label1s)
        #assert len(self._phase0s) == len(self._phase1s)
        #assert len(self._positions) == n
        #assert len(self._phase1s) == n
        for i in range(n):
            #print (self._positions[i])
            fout.write(','.join(str(j) for j in self._positions[i]))
            fout.write("\n")
            fout.write(self._label0s[i])

            fout.write("\n")
            fout.write(','.join(j for j in self._phase0s[i]))
            
            fout.write("\n")
            fout.write(self._label1s[i])
            
            fout.write("\n")
            fout.write(','.join(j for j in self._phase1s[i]))

            fout.write("\n")
            fout.write("\n")
        fout.close()

    def _pre_init(self):

        mutation = []
        insert = []
        delete = []
        stable = []  
        for (refPos, c) in self._columns.items():
            c._set_Lable()
            
            #print ("%s finish set label" % (refPos))
            #sys.exit()
            if c._is_insert == 1:
                insert.append(refPos)
        
            if c._is_mutation == 1:
                mutation.append(refPos)
        
            if c._is_delete == 1:
                delete.append(refPos)
        
            if c._is_stable == 1:
                stable.append(refPos)
        print ("stable len:", len(stable))
        print ("mutation len:", len(mutation))
        
        print ( "delete lens:", len(delete) )
        print ( "insert lens:", len(insert) )
        return stable, mutation, delete, insert

    def __init__(self, columns, contig, obLen = 3):
        self._columns = columns
        self._obLen = obLen
        self._contig = contig
        self._stableRange = [] 
        
        self._snp_mutation = []
        
        self._snp_delete = []
        self._snp_insert = []
        self._coverage = []
        # self.label0 = ""
        # self.label1 = ""
        # self.phase0 = []
        # self.phase1 = [] 

        self._snp = []
     
        #output
        self._positions = []
        self._label0s = []
        self._label1s = []
        self._phase0s = []
        self._phase1s = []

        # output
        # temporary 
    def _homopolyer_filter_delete(self):

        tem = copy.deepcopy(self._snp_delete)

        for referPos in tem:
            m = self._columns[referPos]._map_content
            assert m[0][0] == '*' or m[1][0] == '*'
            if m[0][0] != '*':
                c = m[0][0]
            else:
                c = m[1][0]
            
            if c != self._contig._seq[referPos]:
                return

            assert c == self._contig._seq[referPos]
            if ( tools.same_Character(self._contig._seq[referPos:referPos+3]) 
                    or tools.same_Character(self._contig._seq[referPos-1:referPos+2]) 
                    or tools.same_Character(self._contig._seq[referPos-2:referPos+1]) ):
                self._snp_delete.remove(referPos)

    def _pre_Process(self):
        stable, mutation, delete, insert = self._pre_init() 
      
        #print (len)
        self._stableRange = tools.pos_2_Range(stable)    
        #print (self._stableRange) 
        self._snp_mutation = self._get_SNP(mutation)
        
        # need filter
        self._snp_delete = self._get_SNP(delete)
        self._snp_insert = self._get_SNP(insert)
        
        self._homopolyer_filter_delete()

        # output
        # temporary 
        self._snp = self._snp_mutation
                  
        print ("snp insert len:", len(self._snp_insert))
        print ("snp mutation len:", len(self._snp_mutation))
        print ( "snp delete lens:", len(self._snp_delete) )
        #sys.exit()

    def _get_SNP(self, pos):
        # unfinish 
        snp = [] 
        for a in pos:
            if tools.is_SubRange(a-self._obLen, a-1,self._stableRange) and tools.is_SubRange(a+1,a+self._obLen, self._stableRange):
                    snp.append(a)

        return snp

    def _label_reads(self):

        readsLabel = {} 
        for i in range(len(self._snp)):
            p = self._snp[i]
            content = self._columns[p]._map_content
	    #a = content[0][0]
	    #b = content[1][0]
            for readId in content[0][1]:
                if readId not in readsLabel:
                    readsLabel[readId] = [3]*len(self._snp)
                readsLabel[readId][i] = 0 
            for readId in content[1][1]:
                if readId not in readsLabel:
                    readsLabel[readId] = [3]*len(self._snp)
                readsLabel[readId][i] = 1
            j = 2  
            while j < len(content): 
                for readId in content[j][1]:
                    if readId not in readsLabel:
                        readsLabel[readId] = [3]*len(self._snp)
                    readsLabel[readId][i] = 2
                j += 1    
                    
        # sortFlag = tools.sorted_Map_Value(readsLabel, False)
        # print ("show reads label")
        # for ele in sortFlag:
        #     print (ele)
        return readsLabel       

    def _can_link_phases(self, label0, label1, phase0, phase1, position):
       
 
        print ("check can link or not")
        print (phase0.intersection(self._phase0s[-1])) #1
        print (phase0.intersection(self._phase1s[-1])) #0
        print (phase1.intersection(self._phase0s[-1])) #1
        print (phase1.intersection(self._phase1s[-1])) #0
        
        print (self._positions[-1])
        print (position)
  
        pre_position = self._positions[-1]
        if pre_position[-1] >= position[0]:
            return True
        
        if ( len(phase0.intersection(self._phase0s[-1])) > len(phase0.intersection(self._phase1s[-1])) and	
             len(phase1.intersection(self._phase1s[-1])) > len(phase1.intersection(self._phase0s[-1])) ): 
            return True
        elif ( len(phase0.intersection(self._phase0s[-1])) < len(phase0.intersection(self._phase1s[-1])) and	
             len(phase1.intersection(self._phase1s[-1])) < len(phase1.intersection(self._phase0s[-1])) ):
            return True
        else:
            print ("canot link")
            return False

    def _update_labels(self,labelx, labely, s, phasex, phasey, readsLabel):
 
        self._label0s[-1] += labelx[s:]
        self._label1s[-1] += labely[s:]
        phasex.update(self._phase0s[-1].union(phasex))   
        #self._phase0s[-1].update(self._phase0s[-1].union(phase0))
        phasey.update(self._phase1s[-1].union(phasey))

        unphased = phasex.intersection(phasey)
        if len( unphased ) > 0:
            print ("after link have intersection") 
            print (unphased) 
            print (len(self._label0s[-1]) , self._label0s[-1])
            print (len(self._label1s[-1]) , self._label1s[-1])
            print (len(self._positions[-1]) , self._positions[-1])
            self._re_phasing(unphased, self._label0s[-1], self._label1s[-1] , phasex, phasey, readsLabel, self._positions[-1])
           
        self._phase0s[-1].update(phasex)
        self._phase1s[-1].update(phasey)
      
    def _update_in_link_way(self,label0, label1, phase0, phase1, position, readsLabel):

        print ("update in link way ")  
        pre_position = self._positions[-1]
        if pre_position[-1] >= position[0]:
            for i in range(len(position)):
                if position[i] == pre_position[-1]:
                    self._positions[-1].extend(position[i+1:])
                    print ("link type 1:",i, label0[:i+1], label1[:i+1], self._label0s[-1][-(i+1):],self._label1s[-1][-(i+1):])
                    if ( label0[:i+1] == self._label0s[-1][-(i+1):] and
                         label1[:i+1] == self._label1s[-1][-(i+1):] ):
                        self._update_labels(label0, label1, i+1, phase0, phase1, readsLabel)
                        ''' 
                        self._label0s[-1] += label0[i+1:] 
                        self._label1s[-1] += label1[i+1:]
                        self._phase0s[-1].update(self._phase0s[-1].union(phase0))
                        self._phase1s[-1].update(self._phase1s[-1].union(phase1))
                        '''
                    elif ( label0[:i+1] == self._label1s[-1][-(i+1):] and
                           label1[:i+1] == self._label0s[-1][-(i+1):] ):

                        self._update_labels(label1, label0, i+1, phase1, phase0, readsLabel)
                        '''
                        self._label0s[-1] += label1[i+1:] 
                        self._label1s[-1] += label0[i+1:]
                        self._phase0s[-1].update(self._phase0s[-1].union(phase1))
                        self._phase1s[-1].update(self._phase1s[-1].union(phase0))
                        '''
                    else:
                        print ("link error 1")   
                    return  
        print (self._snp)
        print (position)
        print (pre_position[-1])
         
        s = self._snp.index(pre_position[-1]) + 1
        e = self._snp.index(position[0])
        assert (e>=s)
        self._positions[-1].extend(self._snp[s:e])
        self._positions[-1].extend(position)

       
        self._label0s[-1] +=  (e-s)*'*'
        self._label1s[-1] +=  (e-s)*'*'

        print ("link type 2:")
        print (s, e)
        if ( len(phase0.intersection(self._phase0s[-1])) > len(phase0.intersection(self._phase1s[-1])) and	
             len(phase1.intersection(self._phase1s[-1])) > len(phase1.intersection(self._phase0s[-1])) ):
 
            self._update_labels(label0, label1, 0, phase0, phase1, readsLabel)
            '''
            self._label0s[-1] += label0
            self._label1s[-1] += label1
            #self._phase0s[-1].extend(phase0)
            #self._phase1s[-1].extend(phase1)

            self._phase0s[-1].update(self._phase0s[-1].union(phase0))
            self._phase1s[-1].update(self._phase1s[-1].union(phase1))
            '''
        elif ( len(phase0.intersection(self._phase0s[-1])) < len(phase0.intersection(self._phase1s[-1])) and	
             len(phase1.intersection(self._phase1s[-1])) < len(phase1.intersection(self._phase0s[-1])) ): 
            
            self._update_labels(label1, label0, 0, phase1, phase0, readsLabel)
            '''
            self._label0s[-1] += label1
            self._label1s[-1] += label0
            #self._phase0s[-1].extend(phase1)
            #self._phase1s[-1].extend(phase0)

            self._phase0s[-1].update(self._phase0s[-1].union(phase1))
            self._phase1s[-1].update(self._phase1s[-1].union(phase0))
            '''  
        else: 
            print ("link error 2")

        return
         
 
    def _update(self, label0, label1, phase0, phase1, readsLabel, position):

        if len(position) > 0:
            unphased = phase0.intersection(phase1)
            if len( unphased ) > 0:
                self._re_phasing(unphased, label0, label1, phase0, phase1, readsLabel, position)
           
            if len(self._phase0s) > 0: 
                print ("phasing more longer change 1")
                print (phase0)
                print (phase1)
                print (self._phase0s[-1])
                print (self._phase1s[-1])
                if self._can_link_phases(label0, label1, phase0, phase1, position):
                    self._update_in_link_way(label0, label1, phase0, phase1, position, readsLabel) 
                    return  
            self._label0s.append(label0)
            self._label1s.append(label1)
            self._phase0s.append(phase0)
            self._phase1s.append(phase1)
            self._positions.append(position) 

    def _phasing(self, window=3):
            
        label0 = ""
        label1 = ""
        phase0 = set()
        phase1 = set() 
        readsLabel = self._label_reads()
        self._coverage = (len(self._snp) - window)*[0]
        position = []
        for i in range(len(self._snp)-window):
            phases = {}
            for (read, label) in readsLabel.items():
                f = ''.join( str(j) for j in label[ i :i+window] )
                if f not in phases:
                    phases[f] = []
                phases[f].append(read)
            #print (phases)   
        
            print (self._snp[i:i+3])
            sortedPhases = tools.sorted_Map_Value_Len(phases)
            cov = 0
            allCoverPhases = []
            for (label, reads) in sortedPhases:
                if label.find('3') == -1:
                    #cov += len(reads)
                    self._coverage[i] += len(reads)  
                    allCoverPhases.append((label, len(reads)))
            '''
            if len(allCoverPhases) == 1:
                print ("waiting dealing")
                if len(position) > 0:
                    self._label0s.append(label0)
                    self._label1s.append(label1)
                    self._phase0s.append(phase0) 
                    self._phase1s.append(phase1) 
                    self._positions.append(position) 
                label0 = ""
                label1 = ""
                phase0 = []
                phase1 = [] 
                position = []
            '''
            if ( len(allCoverPhases)>=2 and allCoverPhases[0][1] > cov * 0.2 and allCoverPhases[1][1] > cov * 0.2
               and tools.is_Bool_Reverse(allCoverPhases[0][0], allCoverPhases[1][0]) ):
                label0, label1 = self._phasing_one_window(allCoverPhases, phases, label0, label1, phase0, phase1)
                if len(position) == 0:
                    position.extend(self._snp[i:i+3])
                else:
                    position.extend(self._snp[i+2:i+3])
                #sys.exit()
            else:
 
                self._update(label0, label1, phase0, phase1, readsLabel, position)                
                '''
                if len(position) > 0:

                    unphased = phase0.intersection(phase1)     
                    #print ("intersection:", unphased)
                    if len( unphased ) > 0:
                        self._re_phasing(unphased, label0, label1, phase0, phase1, readsLabel, position)
                    self._label0s.append(label0)
                    self._label1s.append(label1)
                    self._phase0s.append(phase0) 
                    self._phase1s.append(phase1)
                    self._positions.append(position) 
                '''
                #unphased = phase0.intersection(phase1)     
                #print ("after re_phasing intersection:", unphased)
                
                print ("not statified")
                label0 = ""
                label1 = ""
                phase0 = set()
                phase1 = set()
                position = []

                #reclass_Intersection(phase0, phase1, label0, label1, readsFlag)
                #fout.write("%s %s\n" % (label0, phase0))
                #fout.write("%s %s\n" % (label1, phase1))
                #fout.write("\n")
                #sys.exit()

        self._update(label0, label1, phase0, phase1, readsLabel, position)                
        '''
        if len(position) > 0:

            unphased = phase0.intersection(phase1)     
                #print ("intersection:", unphased)
            if len( unphased ) > 0:
                self._re_phasing(unphased, label0, label1, phase0, phase1, readsLabel, position)
            self._label0s.append(label0)
            self._label1s.append(label1)
            self._phase0s.append(phase0) 
            self._phase1s.append(phase1)
            self._positions.append(position) 
        '''
    def _re_phasing(self, unphased, label0, label1, phase0, phase1, readLabel, position):
        (s,e) = tools.get_Range_From_List(position, self._snp)
        for read in unphased:
            readL = ''.join(str(j) for j in readLabel[read][s:e])
            if ( tools.hamming_Distance(label0, readL) < 
                    tools.hamming_Distance(label1, readL) ):
                phase1.remove(read)
            elif ( tools.hamming_Distance(label0, readL) >
                    tools.hamming_Distance(label1, readL) ):
                phase0.remove(read)
            else:
                phase0.remove(read)
                phase1.remove(read)
           
    def _phasing_one_window(self, allCoverPhases, labelReads, label0, label1, phase0, phase1):
     
        #print ("in phasing")
        print (label0 , label1)
        print (allCoverPhases[0][0], allCoverPhases[1][0])
        if len(label0) == 0: 
            label0 = allCoverPhases[0][0]
            label1 = allCoverPhases[1][0] 
        else:
            assert ( label0[-2:] == allCoverPhases[0][0][:2] or  
                   label0[-2:] == allCoverPhases[1][0][:2] )
            if label0[-2:] == allCoverPhases[0][0][:2]: 
                label0 = label0 + allCoverPhases[0][0][-1] 
                label1 = label1 + allCoverPhases[1][0][-1] 
            elif label0[-2:] == allCoverPhases[1][0][:2]:     
                label0 = label0 + allCoverPhases[1][0][-1]
                label1 = label1 + allCoverPhases[0][0][-1]
            else:
                print ("error type1")

   
        for v in allCoverPhases:
            if ( tools.hamming_Distance(v[0], label0[-3:]) < tools.hamming_Distance(v[0], label1[-3:]) ):
                phase0.update(phase0.union(set(labelReads[v[0]])))
            elif ( tools.hamming_Distance(v[0], label0[-3:]) > tools.hamming_Distance(v[0], label1[-3:]) ):
                phase1.update(phase1.union(set(labelReads[v[0]])))
            else:
                print ("same distance", v)
            

        print (label0, len(phase0), phase0)      
        print (label1, len(phase1), phase1)      
        
   
        print ("intersection:", phase0.intersection(phase1) )
        return label0, label1
