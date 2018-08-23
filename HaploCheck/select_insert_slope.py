import os
import sys
import contig
def read_1st(filename):
    f = open(filename, "r")
    val = {}
    poss = {}
    count = 1
    for line in f:
        words = line.split()
        wordsLen = len(words) 
        if count%2 == 1:
            assert wordsLen == 2
            rate = float(words[0])
            val[rate] = int(words[1])
            poss[rate] = []
        else:
            for i in range(wordsLen):
                poss[rate].append(int(words[i]))
        count += 1        
    return val, poss;
def is_same(s):
    c = s[0]
    sLen = len(s)
    for i in range(1,sLen):
        if c!=s[i]:
            return False
    return True    


def find_neighbor(a, b, ref):
    aa = []
    for c in a:
        aa.append((c,"FN"))
    for c in b:
        aa.append((c,"FP"))

    bb = sorted(aa)    
    bbLen = len(bb)
    ans = 0
    for i in range(bbLen-1):
        if bb[i][1] != bb[i+1][1] and (bb[i+1][0]-bb[i][0] <= 2) and is_same(ref._seq[bb[i][0]:bb[i+1][0]+1]):
            if bb[i+1][0]-bb[i][0] == 2:
                print (bb[i], bb[i+1])
                print (ref._seq[bb[i][0]:bb[i+1][0]+1])
            #if ans == 100:
            #    sys.exit()
            ans += 1
    return ans




if __name__ == "__main__":
    
    snp1st, snpPoss = read_1st(sys.argv[1])
    nonSNP1st, nonSNPPoss = read_1st(sys.argv[2])

    '''
    contigs = contig.read_Contig(sys.argv[3])
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
    '''
    slope = 0.1

    print ("slope, FP, FN, FN+FP, totalNode, FP/float(TP+FP), float(FN)/(realNode)")
    
    minFPFN = 50000
    bestCaseFPFN = []
    while slope <= 1:
        FN = 0   # Snp, think it is not
        FP = 0   # nonSnp, think it is
        TP = 0
        FPPos = []
        FNPos = []
        realNode = 0
        slope = round(slope,2)
        for x in snp1st:
            realNode = realNode + snp1st[x] 
            if x < slope:
                #print (x, slope, snp1st[x])
                FN = FN + snp1st[x]
                FNPos.extend(snpPoss[x])
            else:
                TP = TP + snp1st[x]
        for x in nonSNP1st:
            if x >= slope:
                FP = FP + nonSNP1st[x]
                FPPos.extend(nonSNPPoss[x])
                if nonSNP1st[x] != len(nonSNPPoss[x]):
                    print (x)
                    print (nonSNP1st[x])
                    print (len(nonSNPPoss[x]))
                    print ("wrong")

                assert nonSNP1st[x] == len(nonSNPPoss[x])
        if (TP+FP) == 0:
            break    
        print (FN, FP, len(FNPos), len(FPPos))
        assert FN == len(FNPos)
        assert FP == len(FPPos)
        FNPos = sorted(FNPos)
        FPPos = sorted(FPPos)
        #print FNPos
        #print FPPos
        #remove neighbor
        ''' 
        neiNum = find_neighbor(FNPos,FPPos, contig)
        print ("wrong FN/FP number ",neiNum)
     
          
        #sys.exit()
        FP = FP-neiNum 
        FN = FN-neiNum
        TP = TP+neiNum
        '''
        rate1 = FP/float(TP+FP) 
        rate2 = FN/float(realNode)
        if FP+FN < minFPFN:
            bestCaseFPFN = []
            minFPFN = FP+FN
            bestCaseFPFN.append(slope)
            bestCaseFPFN.append(FP)
            bestCaseFPFN.append(FN)
            bestCaseFPFN.append(rate1*100)
            bestCaseFPFN.append(rate2*100)
            #print bestCaseFPFN
        print (slope, FP, FN, FN+FP, rate1, rate2)
        slope = slope + 0.01
    print ("slope FP FN error_rate 1-sensitive")
    for v in bestCaseFPFN:
        print ("%.2f\t" % v, end='')
