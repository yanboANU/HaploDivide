#!/usr/bin/env python

import os
import sys
import string
from itertools import ifilter,imap
import collections

class Sequence(object):
    def __init__(self, name, length, s, e, seq, direction):
        self.name = name
        self.length = int(length)
        self.s = int(s)
        self.e = int(e)
        self.seq = seq   # ATCG- 
        self.seqPos = [] 
        self.dir = direction 
        # seqPos[i] = j , the position in self sequence(have no -) j is reposible to i in alignSeq(self.seq) 
    
    def printSeq(self):
        print self.seq 
    
    def generateSeqPos(self):
        #align seqence to position sequence in original seq 
        if self.dir == '+':
            ss = self.s
            for ii in self.seq:
                self.seqPos.append(ss)
                if ii !=  '-':
                    ss += 1 
        
        if self.dir == '-':
            ss = self.s
            for ii in self.seq:
                self.seqPos.append(self.length-1 - ss)
                if ii !=  '-':
                    ss += 1 
    
    def printSeqPart(self,start,ll):
        print self.seq[start:start+ll] 

    def writeName(self, f):
        f.write("%s %s %s %s"  % (self.name, self.length, self.s, self.e))

    def writeSeq(self, f):
        f.write(self.seq) 

class Alignment(object):
    def __init__(self, left, right, align, leftDir, rightDir):
        self.left = left    # Sequence, align seq (include '-')
        self.right = right  # Sequence
        self.align = align  # ||| ***
        self.score = [0, 0, 0, 0, 0, 0] 
        self.left_dir = leftDir
        self.right_dir = rightDir 
        # socre1, matchNum, insertNum, delNum, gapNum, score2
 
    def writeScore(self, f):
        for ss in self.score:
            f.write("%s " % (ss))     

    def write_m5(self, filename):
        f = open(filename, "w")
        self.left.writeName(f) # name length start end
        f.write(" + ")
        self.right.writeName(f) 
        f.write(" + ")
        self.writeScore(f)
        self.left.writeSeq(f)
        f.write(" %s " % (self.align))
        self.right.writeSeq(f)
        f.close()  
    
    def printAlign(self):
        self.printAlignPart(0, len(self.align))

    def printAlignPart(self,s,l):
        print "\n"
        print self.left.name, self.right.name
        print "in align postion: ", s, s+l
        l = min(l, len(self.align)-s) 
        print self.left.seqPos[s:s+l]  # left seq pos  
        self.left.printSeqPart(s,l)
        print self.align[s:s+l]
        self.right.printSeqPart(s,l)
        print self.right.seqPos[s:s+l]  # right seq pos
        return
 
    def getAlignDiffStat(self, mergeLen):

        diff = {}
        
        checkRange = self.calc_diffPosition(mergeLen)
        for (s,e) in checkRange:
            #for i in range(s,e+1):
            i = s
            while i < e+1: 
	        if self.align[i] == '*':
		    if (self.left.seq[i] != '-' and self.right.seq[i] != '-'):
                        j = self.left.seqPos[i]
		        print (i, j, "mutation")
		        if j not in diff:
			    diff[j] = []
		        diff[j].append('M'+self.left.seq[i]+self.right.seq[i])
                        #i = i+1
             
                
                    if (self.left.seq[i] == '-' and self.right.seq[i] != '-'):
                        j = self.left.seqPos[i]
                        if j not in diff:
                            diff[j] = []
                        start = i
                        while (i+1<e+1 and self.left.seq[i+1] == '-' and self.right.seq[i+1] != '-'):
                            i = i+1
                        end = i+1
                        print (start, j, "insert")
                        diff[j].append('I'+self.right.seq[start:end])
                        

                    if (self.left.seq[i] != '-' and self.right.seq[i] == '-'):
                        j = self.left.seqPos[i]
                        if j not in diff:
                            diff[j] = []
                        start = i
                        while (i+1<e+1 and self.left.seq[i+1] != '-' and self.right.seq[i+1] == '-'):
                            i = i+1
                        end = i+1 
                        print (start, j, "delete")
                        diff[j].append('D'+self.left.seq[start:end])
                i = i+1
        fout = open("actual.diff.output","w") 
        for (key, value) in sorted(diff.items()):
            print (key, value)
            #assert len(value) == 1 
            fout.write("%s" % (key)) 
            for v in value:
                fout.write(" %s" % (v))
            fout.write("\n") 

        return diff     

    def calc_diffPosition(self, mergeLen):  
 
        diffPos = []
        for i in range(len(self.align)):
            if self.align[i] == '*':
                diffPos.append(i)
        checkRange = []
        i = 0
        print "Number dismatch:", len(diffPos)
        while i < len(diffPos)-1: 
            s = diffPos[i]
            assert ( (i +1) < len(diffPos) )
            while (diffPos[i+1] <= (diffPos[i] + mergeLen)):
                i = i+1
                if i >= len(diffPos)-1:
                    break
            #if i > len(diffPos):
            i = min(i, len(diffPos)-1)  
            e = diffPos[i]
            checkRange.append((s,e))
            i = i+1
        print "position of dismatch:", diffPos
        print "dismatch range:", checkRange
        print "number of dismatch range:", len(checkRange)
        return checkRange  
       


    def printAlignDiff(self, mergeLen):
        check = self.calc_diffPosition(mergeLen)
        for (s,e) in check:
            self.printAlignPart(max(s-mergeLen,0),e-s+2*mergeLen)
        return 

    #######################
    # move gaps togther, push gap to right
    #######################
    def PreProcess(self):
        assert len(self.left.seq) == len(self.right.seq)
        alignLen = len(self.left.seq)
        assert len(self.align) == alignLen
        qstr = list(self.left.seq)
        tstr = list(self.right.seq)
        #print collections.Counter(qstr)
        #print collections.Counter(tstr) 

        # push gap to the right, but not past the end

        assert tstr[0] != '-'
        assert qstr[0] != '-'  
        for i in range(alignLen-1):
            # pushing target gap
            if tstr[i] == '-' and tstr[i-1] != '-':
                j = i + 1
                while (j  < alignLen):
                    c = tstr[j]
                    #print "tstr[j]", j, tstr[j]
                    if (c != '-'):
                        if (c == qstr[i]):
                            tstr[i] = c
                            tstr[j] = '-'
                            #print ''.join(qstr)     
                            #print ''.join(tstr)
                        break
                    j += 1
            # pushing query gap
            if qstr[i] == '-' and qstr[i-1] != '-':
                j = i + 1
                while (j < alignLen):
                    c = qstr[j]
                    #print "qstr[j]", j, qstr[j]
                    if (c != '-'):
                        if (c == tstr[i]):
                            qstr[i] = c
                            qstr[j] = '-'
                            #print ''.join(qstr)     
                            #print ''.join(tstr)
                        break
                    j += 1
        '''
        print ''.join(qstr)     
        print ''.join(tstr)
        '''
        self.left.seq = ''.join(qstr)     
        self.right.seq = ''.join(tstr)
        newAlign = ""
        for i in range(len(self.left.seq)):
            if self.left.seq[i] == self.right.seq[i]:
                newAlign = newAlign + '|'
            else:   
                newAlign = newAlign + '*'
        self.align = newAlign
        '''
        print self.left.seq
        print self.align
        print self.right.seq   
        '''

    def calcMatchRate(self, ll, rr):

        match = 0
        for c in self.align[ll:rr]:
	    if c == '|':
	        match += 1    
        matchRate =  match/float(rr-ll)
        print ("(%s, %s): %s/%s = %s" % (ll, rr, match, rr-ll, matchRate)) 
        return matchRate

##########################
# same function calc_diffPosition
##########################

def giveCheckPos(align, mergeLen):
    diffPos = []
    for i in range(len(align)):
        if align[i] == '*':
            diffPos.append(i)
    checkRange = []
    i = 0
    print "Number dismatch:", len(diffPos)
    while i < len(diffPos)-1: 
        s = diffPos[i]
        assert ( (i +1) < len(diffPos) )
        while (diffPos[i+1] <= (diffPos[i] + mergeLen)):
            i = i+1
            if i >= len(diffPos)-1:
                break
        #if i > len(diffPos):
        i = min(i, len(diffPos)-1)  
        e = diffPos[i]
        checkRange.append((s,e))
        i = i+1
    print "position of dismatch:", diffPos
    print "dismatch range:", checkRange
       
    return checkRange  
       


def read_blasr_m5(fileName):

    f = open(fileName,"r")
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        words = line.split()

        #print words[:4] 
        queryName, queryLen, queryS, queryE = words[:4]
        queryDirection = words[4]
        targetName, targetLen, targetS, targetE = words[5:5+4]
        targetDirection = words[9] 
 
        querySeq = words[-3]   # ATCG-
        align = words[-2]      # |,*
        targetSeq = words[-1]  # ATCG-    
 
    f.close()
    query = Sequence(queryName, queryLen, queryS, queryE, querySeq, queryDirection) 
    target = Sequence(targetName, targetLen, targetS, targetE, targetSeq, targetDirection)
    alignObj = Alignment(query, target, align, queryDirection, targetDirection)  
    return alignObj
 


if __name__ == "__main__":
    ''' 
    if(len(sys.argv) <= 4):
        print "Usage: python calc_swith.py mutation_ref_m5.blasr contig_ref_m5.blasr contig_mutation_m5.blasr mergeLen"   
    '''
    # blasr ref mutation output(m5) 
    # blasr ref contig
    # blasr mutation contig
    # get swith number    

    alignRM = read_blasr_m5(sys.argv[1])
    #alignRM.printAlign()
    alignRM.PreProcess() 
    #alignRM.write_m5(sys.argv[1] + "_shift")
    alignRM.left.generateSeqPos()
    alignRM.right.generateSeqPos()    

    alignRM.printAlign() 
    
    check = []   
    ''' #naive way give check
    for i in range(len(alignRM.align)):
        if alignRM.align[i] == '*':
            check.append(i)
    '''
   
    '''  
    mergeLen = int(sys.argv[4])
    check = giveCheckPos(alignRM.align, mergeLen)
        
    alignRC = read_blasr_m5(sys.argv[2])    
    #alignRC.printAlign()
    alignRC.PreProcess() 
   
    alignRC.left.generateSeqPos()
    alignRC.right.generateSeqPos()
    
    alignRC.printAlign()

    #alignRC.write_m5(sys.argv[2] + "_shift")
    alignMC = read_blasr_m5(sys.argv[3])

    #alignMC.printAlign()
    alignMC.PreProcess()

    alignMC.left.generateSeqPos()
    alignMC.right.generateSeqPos()

    alignMC.printAlign()
    #alignMC.write_m5(sys.argv[3] + "_shift")

    
    RorM = {}   
    for i in range(len(check)):
        #print check[i]

        RseqPosStart = alignRM.left.seqPos[check[i][0]]
        RseqPosEnd = alignRM.left.seqPos[check[i][1]]

        MseqPosStart = alignRM.right.seqPos[check[i][0]] 
        MseqPosEnd = alignRM.right.seqPos[check[i][1]] 
        
        if RseqPosStart > alignRC.left.s and MseqPosStart > alignMC.left.s and RseqPosEnd < alignRC.left.e and MseqPosEnd < alignMC.left.e:
	    l = max(check[i][0]-mergeLen, 0)
            r = min(check[i][1]+mergeLen, len(alignRM.align)) 
            #r = check[i]+5  

            RseqPosStart = alignRM.left.seqPos[l]
            RseqPosEnd = alignRM.left.seqPos[r]

            MseqPosStart = alignRM.right.seqPos[l] 
            MseqPosEnd = alignRM.right.seqPos[r] 
	    
            print ("(%s , %s) , (%s, %s), (%s, %s)"% (l , r, RseqPosStart, RseqPosEnd, MseqPosStart, MseqPosEnd))
            #print alignRM.left.seqPos[l:l+20]
            alignRM.printAlignPart(l, r-l)
	    #print alignRM.right.seqPos[l:l+20]
           
            # l is align sequence position, RseqPos is ref seq position, MseqPos is mutation seq position
            RCl = alignRC.left.seqPos.index(RseqPosStart) 
            RCr = alignRC.left.seqPos.index(RseqPosEnd) 

            MCl = alignMC.left.seqPos.index(MseqPosStart)
            MCr = alignMC.left.seqPos.index(MseqPosEnd)
            
	    #print alignRC.left.seqPos[ll:ll+20]
            alignRC.printAlignPart(RCl, RCr-RCl)
            #print alignRC.right.seqPos[ll:ll+20]
            rateRC =  alignRC.calcMatchRate(RCl, RCr)

	    #print alignMC.left.seqPos[rr:rr+20]
            alignMC.printAlignPart(MCl, MCr-MCl)            
	    #print alignMC.right.seqPos[rr:rr+20]
            rateMC =  alignMC.calcMatchRate(MCl, MCr)
            
            if rateRC > rateMC:
                print "choose Ref\n"
                RorM[(l,r)] = "R"
            elif rateRC == rateMC: 
                print "hard to say\n"
                RorM[(l,r)] = "N"
            else:
                print "choose Mutation\n"
                RorM[(l,r)] = "M"
        
    res = ""
    for key in RorM:
        res = res + RorM[key] 
    print res 
    print collections.Counter(res)
    switchNum = 0
    for i in range(len(res)-1):
        if res[i+1] != res[i] and res[i+1] != 'N' and res[i] != 'N':
            switchNum += 1
    print "switchNum number: ", switchNum
    print "total number: ", len(res)
    ''' 


