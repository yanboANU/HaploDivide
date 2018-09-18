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


    def _write_result2(self, start, end):
        fout = open(self._contig._name+"_phasing_result_" + str(start) + "_" + str(end),"w")
        n =  len(self._label0s)
        fout.write("name: %s len: %s\n" % (self._contig._name, end - start))
        fout.write("snp mutation number: %s\n" % (len(self._snp_mutation)))
        fout.write("snp mutation rate: %s\n" % (len(self._snp_mutation)/float(end-start)))
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




    def _write_result(self):
        fout = open(self._contig._name+"_phasing_result","w")
        n =  len(self._label0s)
        fout.write("name: %s len: %s\n" % (self._contig._name, self._contig._len))
        fout.write("snp mutation number: %s\n" % (len(self._snp_mutation)))
        fout.write("snp mutation rate: %s\n" % (len(self._snp_mutation)/float(self._contig._len)))
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

    def _pre_init(self, start, end, contigLen, averageCov):

        mutation = []
        insert = []
        delete = []
        stable = []  
        it = self._columns.items()

        #remove head and tail 10k,  
        if start == 0:
            start = 10000;
        
        if end == contigLen:
            assert end > 10000
            end = contigLen - 10000 
        #for (refPos, c) in self._columns.items():
        for (refPos, c) in it:
            c._is_insert = 0
            c._is_mutation = 0
            self._is_delete = 0
            self._is_stable = 0    
            if refPos < start or refPos > end or c._cov<7 or c._cov>=2*averageCov:
                continue
            #print ("%s finish set label" % (refPos))
            #sys.exit() 
            c._set_Lable()
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

    def _pre_Process(self, start, end, contigLen, averageCov):

        time1 = time.clock()
        stable, mutation, delete, insert = self._pre_init(start, end, contigLen, averageCov) 
        time2 = time.clock()

        print ("snp insert len:", len(insert))
        print ("snp len:", len(mutation))
        print ( "snp delete lens:", len(delete) )
        print ( "pre_init running %s Seconds" % (time2 - time1) )
        #print (len)
        if self._obLen == 0:  #in this case, don't need to ensure neighours in stable range
            self._snp_mutation, self._snp_insert, self._snp_delete = mutation, insert, delete

        else:    
            #self._stableRange = tools.pos_2_Range(stable)    
            #print (self._stableRange)
            #sys.exit()
            #print (mutation)
            #print (insert)
            #print (delete)
            self._snp_mutation, self._snp_insert, self._snp_delete = self._get_SNP(mutation, insert, delete)
            #print (self._snp_mutation)
            #print (self._snp_insert)
            #print (self._snp_delete)
        time3 = time.clock()
        
        print ( "get SNP part running %s Seconds" % (time3 - time2) )
        # need filter, now version, don't use _snp_delete and _snp_insert
        # self._snp_delete = self._get_SNP(delete)
        # self._snp_insert = self._get_SNP(insert)
        #self._homopolyer_filter_delete() 
        # output
        # temporary 
        #change July, eight
        self._snp = self._snp_mutation
        
        print ("snp insert len:", len(self._snp_insert))
        print ("snp len:", len(self._snp_mutation))
        print ( "snp delete lens:", len(self._snp_delete) )

        self._write_SNP(start, end)
        self._write_delete(start, end)
        self._write_insert(start, end)
    

    # in the following three reading, reference start with 1
    def _write_insert(self, start, end):
        fout = open(self._contig._name+"_insert_" + str(start) + "_" + str(end) ,"w")
        fout.write("name: %s len: %s\n" % (self._contig._name, self._contig._len))
        fout.write("insert number: %s \n" % (len(self._snp_insert)))
        for sm in self._snp_insert:
            assert sm in self._columns
            fout.write("%s %s %s" % (sm+1, self._columns[sm]._insert_content, self._columns[sm]._insert_number))
            #fout.write("%s %s" % (self._columns[sm]._insert_content[1][0], len(self._columns[sm]._insert_content[1][1])))
            #fout.write(" %s" % (self._columns[sm]._cov))
            fout.write("\n")
        fout.close()
   
    
    def _write_delete(self, start, end):
        fout = open(self._contig._name+"_delete_" + str(start) + "_" + str(end) ,"w")
        fout.write("name: %s len: %s\n" % (self._contig._name, self._contig._len))
        fout.write("snp mutation number: %s \n" % (len(self._snp_delete)))
        for sm in self._snp_delete:
            # reference start with 1, and dbSNP150, dbSNP158 both record deletion at position 11
            # A   T  C      
            # 11  12 13
            # A   T  C    reads1
            # A   *  C    reads2
            fout.write("%s %s %s " % (sm, self._columns[sm]._map_content[0][0], len(self._columns[sm]._map_content[0][1])))
            if len(self._columns[sm]._map_content) >= 2:
                fout.write("%s %s" % (self._columns[sm]._map_content[1][0], len(self._columns[sm]._map_content[1][1])))
                #fout.write(" %s" % (self._columns[sm]._cov))
                fout.write("\n")
            else:
                fout.write("\n")
        fout.close()
    
    def _write_SNP(self, start, end):
        fout = open(self._contig._name+"_snp_mutation_" + str(start) + "_" + str(end) ,"w")
        fout.write("name: %s len: %s\n" % (self._contig._name, self._contig._len))
        fout.write("snp mutation number: %s \n" % (len(self._snp_mutation)))
        for sm in self._snp_mutation:
            fout.write("%s %s %s " % (sm+1, self._columns[sm]._map_content[0][0], len(self._columns[sm]._map_content[0][1])))
            fout.write("%s %s" % (self._columns[sm]._map_content[1][0], len(self._columns[sm]._map_content[1][1])))
            #fout.write(" %s" % (self._columns[sm]._cov))
            fout.write("\n")
        fout.close()
    

    def _get_SNP(self, mutate, insert, delete):
        # unfinish 
        if self._obLen == 0:
            return pos
        snp = []
        for pos in mutate:
            snp.append((pos, 'M'))
        for pos in insert:
            snp.append((pos, 'I'))
        for pos in delete:
            snp.append((pos, 'D'))
        snp = sorted(snp)
        snpLen = len(snp)
        print ("snp+indel length:", snpLen)
        i = 0
        mutate = []
        insert = []
        delete = []
        while i+1 < snpLen:
            #print (i)
            if i+1 < snpLen and snp[i][0] < snp[i+1][0] - self._obLen:
                if snp[i][1] == 'M': 
                    mutate.append(snp[i][0])    
                if snp[i][1] == 'I': 
                    insert.append(snp[i][0])
                if snp[i][1] == 'D': 
                    delete.append(snp[i][0])
                i = i+1
                continue
            while i+1<snpLen and snp[i][0] >= snp[i+1][0] - self._obLen:
                i = i+1
            i += 1    
        # deal last element
        i=snpLen-1 
        if snp[i][0] > snp[i-1][0] + self._obLen:
            if snp[i][1] == 'M': 
                mutate.append(snp[i][0])    
            if snp[i][1] == 'I': 
                insert.append(snp[i][0])
            if snp[i][1] == 'D': 
                delete.append(snp[i][0])
        return mutate, insert, delete

       
    '''   
    def _get_SNP(self, pos):
        # unfinish 
        if self._obLen == 0:
            return pos

        snp = []
        for a in pos:
            if tools.is_SubRange(a-self._obLen, a-1,self._stableRange) and tools.is_SubRange(a+1,a+self._obLen, self._stableRange):
                    snp.append(a)
        return snp
    '''
    def _label_reads(self):

        readsLabel = {}
        lenSNP = len(self._snp)

        # self._obLen is different with window
        # self._obLen used in looking SNP
        for i in range(lenSNP):
            p = self._snp[i]
            content = self._columns[p]._map_content
	    #a = content[0][0]
	    #b = content[1][0]
            for readId in content[0][1]:
                if readId not in readsLabel:
                    readsLabel[readId] = [3]*lenSNP
                readsLabel[readId][i] = 0 
            for readId in content[1][1]:
                if readId not in readsLabel:
                    readsLabel[readId] = [3]*lenSNP
                readsLabel[readId][i] = 1    
            j = 2
            lenContent = len(content) 
            #while j < len(content):
            while j < lenContent: 
                for readId in content[j][1]:
                    if readId not in readsLabel:
                        readsLabel[readId] = [3]*lenSNP
                    readsLabel[readId][i] = 2
                j += 1    
         
        #fout = open("readsLabel", "w")
        fout = open(self._contig._name + "readsLabel", "w") 
        for readId in readsLabel:
            coverRange = tools.get_Cover_Range(readsLabel[readId])
            coverLength = coverRange[1] - coverRange[0] 
             
            coverLabel = readsLabel[readId][coverRange[0]:coverRange[1]]
            fout.write("readId: %s coverLength %s 0/1 length %s\n" % (readId, coverLength, tools.count01(coverLabel)))
            fout.write( ''.join(str(c) for c in coverLabel ))
            fout.write("\n") 

        # sortFlag = tools.sorted_Map_Value(readsLabel, False)
        # print ("show reads label")
        # for ele in sortFlag:
        #     print (ele)
        return readsLabel 

    '''
    def _can_link_phases_loose_way(self, label0, label1, phase0, phase1, position):
       
 
        print ("check can link or not")
        print ("0 intersection new 0", len(phase0.intersection(self._phase0s[-1]))) #1
        print ("0 intersection new 1", len(phase0.intersection(self._phase1s[-1]))) #0
        print ("1 intersection new 0", len(phase1.intersection(self._phase0s[-1]))) #1
        print ("1 intersection new 1", len(phase1.intersection(self._phase1s[-1]))) #0
        
        print ("pre pos:", self._positions[-1])
        print ("new pos:", position)
        
        _0id = tools.only0( len(phase0.intersection(self._phase0s[-1])), len(phase0.intersection(self._phase1s[-1])) ,	
             len(phase1.intersection(self._phase1s[-1])) , len(phase1.intersection(self._phase0s[-1])) )
        # three are 0, only one not 0
        if _0id > 0:
            return True
        #
        else:
            print ("canot loose link")
            return False
    '''



    def _can_link_phases(self, label0, label1, phase0, phase1, position):
       

        lenIntersection00 = len(phase0.intersection(self._phase0s[-1]))
        lenIntersection01 = len(phase0.intersection(self._phase1s[-1]))
        lenIntersection10 = len(phase1.intersection(self._phase0s[-1]))
        lenIntersection11 = len(phase1.intersection(self._phase1s[-1]))
        print ("check can link or not")
        print ("0 intersection new 0", lenIntersection00) #1
        print ("0 intersection new 1", lenIntersection01) #0
        print ("1 intersection new 0", lenIntersection10) #1
        print ("1 intersection new 1", lenIntersection11) #0
        
        print ("pre pos:", self._positions[-1])
        print ("new pos:", position)
  
        pre_position = self._positions[-1]
        if pre_position[-1] >= position[0]:
            return True
       
        if lenIntersection00 == 0 and lenIntersection01 == 0 and lenIntersection10 == 0 and lenIntersection11 == 0:
            print ("canot link both 0")
            return False

        if ( lenIntersection00 >= lenIntersection01 and	
             lenIntersection11 >= lenIntersection10 ): 
            return True
        elif ( lenIntersection00 <= lenIntersection01 and	
               lenIntersection11 <= lenIntersection10 ):
            return True
        else:
            print ("canot link")
            return False

    def _update_labels(self,labelx, labely, s, phasex, phasey, readsLabel):
 
        self._label0s[-1] += labelx[s:]
        self._label1s[-1] += labely[s:]
        #phasex = self._phase0s[-1].union(phasex)   
        
        phasex.update(self._phase0s[-1])
        #self._phase0s[-1].update(self._phase0s[-1].union(phase0))
        #phasey = self._phase1s[-1].union(phasey)
        phasey.update(self._phase1s[-1])

        unphased = phasex.intersection(phasey)
        if len( unphased ) > 0:
            print ("after link have intersection") 
            print (unphased) 
            print (len(self._label0s[-1]) , self._label0s[-1])
            print (len(self._label1s[-1]) , self._label1s[-1])
            print ("before intersection check after link", len(phasex) , len(phasey))
            self._re_phasing(unphased, self._label0s[-1], self._label1s[-1] , phasex, phasey, readsLabel, self._positions[-1])
           
            print ("after intersection check after link1 ", len(phasex) , len(phasey))
        assert len(phasex.intersection(phasey)) == 0

        print ("after intersection check after link2", len(phasex) , len(phasey))
        # a.update(b) # b small than a
        # a.update(b) is add element
        self._phase0s[-1] = phasex
        self._phase1s[-1] = phasey
        print ("after intersection check after link3", len(self._phase0s[-1]) , len(self._phase1s[-1]))

    def _update_in_link_way(self,label0, label1, phase0, phase1, position, readsLabel):

        print ("update in link way ")  
        pre_position = self._positions[-1]
        if pre_position[-1] >= position[0]:
            lenP = len(position)
            #for i in range(len(position)):
            for i in range(lenP):
                if position[i] == pre_position[-1]:
                    self._positions[-1].extend(position[i+1:])
                    print ("link type 1:",i, label0[:i+1], label1[:i+1], self._label0s[-1][-(i+1):],self._label1s[-1][-(i+1):])
                    if ( label0[:i+1] == self._label0s[-1][-(i+1):] and
                         label1[:i+1] == self._label1s[-1][-(i+1):] ):
                        self._update_labels(label0, label1, i+1, phase0, phase1, readsLabel)
                    elif ( label0[:i+1] == self._label1s[-1][-(i+1):] and
                           label1[:i+1] == self._label0s[-1][-(i+1):] ):

                        self._update_labels(label1, label0, i+1, phase1, phase0, readsLabel)
                    else:
                        print ("link error 1")   
                    return  
        #print (self._snp)
        #print (position)
        #print (pre_position[-1])
         
        s = self._snp.index(pre_position[-1]) + 1
        e = self._snp.index(position[0])
        assert (e>=s)
        self._positions[-1].extend(self._snp[s:e])
        self._positions[-1].extend(position)

       
        self._label0s[-1] +=  (e-s)*'*'
        self._label1s[-1] +=  (e-s)*'*'

        print ("link type 2:")
        print (s, e)

        lenIntersection00 = len(phase0.intersection(self._phase0s[-1]))
        lenIntersection01 = len(phase0.intersection(self._phase1s[-1]))
        lenIntersection10 = len(phase1.intersection(self._phase0s[-1]))
        lenIntersection11 = len(phase1.intersection(self._phase1s[-1]))

        #if ( len(phase0.intersection(self._phase0s[-1])) >= len(phase0.intersection(self._phase1s[-1])) and	
        #     len(phase1.intersection(self._phase1s[-1])) >= len(phase1.intersection(self._phase0s[-1])) ):

        if ( lenIntersection00 >= lenIntersection01  and 	
             lenIntersection11 >= lenIntersection10 ):
            self._update_labels(label0, label1, 0, phase0, phase1, readsLabel)
        #elif ( len(phase0.intersection(self._phase0s[-1])) <= len(phase0.intersection(self._phase1s[-1])) and	
        #     len(phase1.intersection(self._phase1s[-1])) <= len(phase1.intersection(self._phase0s[-1])) ):     
        elif ( lenIntersection00 <= lenIntersection01 and	
               lenIntersection11 <= lenIntersection10 ): 
            self._update_labels(label1, label0, 0, phase1, phase0, readsLabel)
        else: 
            print ("link error 2")

        return
         
 
    def _update(self, label0, label1, phase0, phase1, readsLabel, position):

        if len(position) > 0:
            print ("before ratio check ",len(phase0), len(phase1) )

            self._cover_01_check(label0, label1, phase0, phase1, position, readsLabel)             
            #self._cover_ratio_check(label0, label1, phase0, phase1, position, readsLabel)             
            
            print ("after ratio check ",len(phase0), len(phase1) )
            unphased = phase0.intersection(phase1)
            if len( unphased ) > 0:
                self._re_phasing(unphased, label0, label1, phase0, phase1, readsLabel, position)
          
                print ("after intersection check ",len(phase0), len(phase1) )
                
            if len(self._phase0s) > 0: 
                print ("phasing more longer change 1")
                if self._can_link_phases(label0, label1, phase0, phase1, position):
                    print ("link")
                    self._update_in_link_way(label0, label1, phase0, phase1, position, readsLabel) 
                    return  
                '''  
                if self._can_link_phases_loose_way(label0, label1, phase0, phase1, position):
                    print ("loose link")
                    self._update_in_loose_link_way(label0, label1, phase0, phase1, position, readsLabel) 
                    return  
                '''
            self._label0s.append(label0)
            self._label1s.append(label1)
            assert len(phase0.intersection(phase1)) == 0
            self._phase0s.append(phase0)
            self._phase1s.append(phase1)
            self._positions.append(position) 

    def _phasing(self, window=3):
            
        label0 = ""
        label1 = ""
        phase0 = set()
        phase1 = set() 
        readsLabel = self._label_reads()
        
        lenSNP = len(self._snp)
        self._coverage = (lenSNP - window + 1)*[0]
        position = []
        #for i in range(len(self._snp) - window + 1):
        print ("window distrubution")
        Wout = open( str(window)+"mer_distrubution","w" )
        FSout= open("1st_2nd_"+str(window)+"mer_distrubution", "w")
        for i in range(lenSNP - window + 1):
            phases = {}
            # every position go through all reads, not good 
            '''
            for (read, label) in readsLabel.items():
                f = ''.join( str(j) for j in label[ i :i+window] )
                if f not in phases:
                    phases[f] = []
                phases[f].append(read)
            ''' 

            # only consider these position correspond reads
            p = self._snp[i]
            content = self._columns[p]._map_content
            lenContent = len(content) 
            j = 0
            while j < lenContent: 
                for readId in content[j][1]:
                    label = readsLabel[ readId ]
                    f = ''.join( str(c) for c in label[ i :i+window] )
                    if f not in phases:
                        phases[f] = []
                    phases[f].append(readId)
                j = j + 1    


            sortedPhases = tools.sorted_Map_Value_Len(phases)
            cov = 0
            allCoverPhases = []
            Wout.write( ("%d %d") % (i,p) )
            for (label, reads) in sortedPhases:
                if label.find('3') == -1:
                    #cov += len(reads)
                    lenReads= len(reads)
                    self._coverage[i] += lenReads  
                    allCoverPhases.append((label, lenReads))
                    Wout.write(" %s %d" % (label, lenReads))
            Wout.write("\n")

            if len(allCoverPhases)>=2:
                FSout.write("%d %d %.2f %.2f\n" % (i, p, float(allCoverPhases[0][1]) /self._coverage[i] , float(allCoverPhases[1][1]/self._coverage[i])))
            elif len(allCoverPhases) == 1:
                FSout.write("%d %d %.2f 0\n" % (i, p, float(allCoverPhases[0][1]) /self._coverage[i]))

            if ( len(allCoverPhases)>=2 and allCoverPhases[0][1] > cov * 0.2 and allCoverPhases[1][1] > cov * 0.2
               and tools.is_Bool_Reverse(allCoverPhases[0][0], allCoverPhases[1][0]) and 
                 (len(label0) == 0 or label0[-2:] == allCoverPhases[0][0][:2] or label0[-2:] == allCoverPhases[1][0][:2] ) ):
                label0, label1 = self._phasing_one_window(allCoverPhases, phases, label0, label1, phase0, phase1)
                if len(position) == 0:
                    position.extend(self._snp[i:i+3])
                else:
                    position.extend(self._snp[i+2:i+3])
                #sys.exit()
            else:
 
                print ("before update ",len(phase0), len(phase1) )
                self._update(label0, label1, phase0, phase1, readsLabel, position)                
                #unphased = phase0.intersection(phase1)     
                #print ("after re_phasing intersection:", unphased)
                
                #print ("not statified")
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
        Wout.close() 

    def _cover_ratio_check(self, label0, label1, phase0, phase1, position, readsLabel):
        (s,e) = tools.get_Range_From_List(position, self._snp)
        remove = set()
        for read in phase0: 
            # read cover range
            coverRange = tools.get_Cover_Range(readsLabel[read])
            
            if read == "S1_2737" or read == "S1_6599" or read == "S1_3419":
            #if read == "S1_298" or read == "S1_5935" or read == "S1_4331":
                print (read, s, e, "read cover:", coverRange)
                print ("readsLabel: for check coverRange rigth or not")
                print (readsLabel[read])
                readL = ''.join(str(j) for j in readsLabel[read][s:e])
                print ("read ", readL)
                print ("label", label0)
            
            # (6, 12)  phasing (10, 15), remove read only cover head part
            coverLength = coverRange[1] - coverRange[0]
            coverLabel = readsLabel[read][coverRange[0]:coverRange[1]]
            if tools.count01(coverLabel) < 3:
                remove.add(read)
                print ("remove ", read, "becase cover label:", coverLabel) 
                continue 
    
            if (coverRange[0] < s and coverRange[1] > s  and  coverRange[1] < e and 
                  (coverRange[1] - s) < coverLength * 0.8 ):           
                remove.add(read)
                continue
           
            # (12, 20)  phasing (10, 15), remove read only cover tail part
            if (coverRange[0] < e and coverRange[1] > e and 
                  (e - coverRange[0]) < coverLength * 0.8 ):           
                remove.add(read)
        for r in remove:
            phase0.remove(r)   
        remove = set()
        for read in phase1: 
            # (6, 12)  phasing (10, 15)
            coverRange = tools.get_Cover_Range(readsLabel[read])
            coverLength = coverRange[1] - coverRange[0] 
            
            if read == "S1_2737" or read == "S1_6599" or read == "S1_3419":
                print (read, s, e, "read cover", coverRange)
                print ("readsLabel: for check coverRange rigth or not")
                print (readsLabel[read])
                readL = ''.join(str(j) for j in readsLabel[read][s:e])
                print ("read ", readL)
                print ("label", label1)
           
             
            coverLabel = readsLabel[read][coverRange[0]:coverRange[1]]
            if tools.count01(coverLabel) < 3:
                remove.add(read)
                print ("remove ", read, "becase cover label:", coverLabel) 
                continue 
 
            if (coverRange[0] < s and coverRange[1] > s  and  coverRange[1] < e and 
                  (coverRange[1] - s ) < coverLength * 0.8 ):            
                remove.add(read)
                continue
                
            if (coverRange[0] < e and coverRange[1] > e and 
                  (e - coverRange[0] ) < coverLength * 0.8 ):            
                remove.add(read)
        for r in remove:
            phase1.remove(r)   

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
        #print (label0 , label1)
        #print (allCoverPhases[0][0], allCoverPhases[1][0])
        if len(label0) == 0: 
            label0 = allCoverPhases[0][0]
            label1 = allCoverPhases[1][0] 
        else:
            # check sufix == prefix 
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
                #phase0 = phase0.union(set(labelReads[v[0]]))
                # must update, using = , error
                phase0.update(set(labelReads[v[0]])) 
            elif ( tools.hamming_Distance(v[0], label0[-3:]) > tools.hamming_Distance(v[0], label1[-3:]) ):
                #phase1 = phase1.union(set(labelReads[v[0]]))
                phase1.update(set(labelReads[v[0]])) 
            #else:
                #print ("same distance", v)

        #print (label0, len(phase0), phase0)      
        #print (label1, len(phase1), phase1)      
        #print ("intersection:", phase0.intersection(phase1) )
        return label0, label1

    def _cover_01_check(self, label0, label1, phase0, phase1, position, readsLabel):
        (s,e) = tools.get_Range_From_List(position, self._snp)
        remove = set()
        for read in phase0: 
            # read cover range
            coverRange = tools.get_Cover_Range(readsLabel[read])
            coverLength = coverRange[1] - coverRange[0]
            coverLabel = readsLabel[read][coverRange[0]:coverRange[1]]
            if tools.count01(coverLabel) < 3:
                remove.add(read)
                print ("remove ", read, "becase cover label:", coverLabel) 
    
        for r in remove:
            phase0.remove(r)   
        remove = set()
        for read in phase1: 
            # (6, 12)  phasing (10, 15)
            coverRange = tools.get_Cover_Range(readsLabel[read])
            coverLength = coverRange[1] - coverRange[0] 
             
            coverLabel = readsLabel[read][coverRange[0]:coverRange[1]]

            if tools.count01(coverLabel) < 3:
                remove.add(read)
                print ("remove ", read, "becase cover label:", coverLabel) 
    
        for r in remove:
            phase1.remove(r)   
