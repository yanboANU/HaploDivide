import os
import sys

def read_1st(filename):
    f = open(filename, "r")
    val = {}
    poss = {}
    for line in f:
        words = line.split()
        wordsLen = len(words) 
        if wordsLen == 2:
            rate = float(words[0])
            val[rate] = int(words[1])
            poss[rate] = []
        else:
            for i in range(wordsLen):
                poss[rate].append(int(words[i]))
                
    return val, poss;

def find_neighbor(a, b):
    aa = []
    for c in a:
        aa.append((c,"FN"))
    for c in b:
        aa.append((c,"FP"))

    bb = sorted(aa)    
    bbLen = len(bb)
    ans = 0
    for i in range(bbLen-1):
        if bb[i][1] != bb[i+1][1] and (bb[i+1][0]-bb[i][0] <= 2):
            print bb[i], bb[i+1]
            ans += 1
    return ans




if __name__ == "__main__":
    
    snp1st, snpPoss = read_1st(sys.argv[1])
    nonSNP1st, nonSNPPoss = read_1st(sys.argv[2])

    slope = 0.1

    print "slope, FP, FN, FN+FP, totalNode, FP/float(TP+FP), float(FN)/(realNode)"
    
    minFPFN = 50000
    bestCaseFPFN = []
    while slope <= 1:
        FN = 0   # Snp, think it is not
        FP = 0   # nonSnp, think it is
        TP = 0
        FPPos = []
        FNPos = []
        realNode = 0
        for x in snp1st:
            realNode = realNode + snp1st[x] 
            if x < slope:
                FN = FN + snp1st[x]
                FNPos.extend(snpPoss[x])
            else:
                TP = TP + snp1st[x]
        for x in nonSNP1st:
            if x >= slope:
                FP = FP + nonSNP1st[x]
                FPPos.extend(nonSNPPoss[x])
        if (TP+FP) == 0:
            break    
        print FN, FP, len(FNPos), len(FPPos)
        assert FN == len(FNPos)
        assert FP == len(FPPos)
        FNPos = sorted(FNPos)
        FPPos = sorted(FPPos)
        #print FNPos
        #print FPPos
        #remove neighbor
        
        #neiNum = find_neighbor(FNPos,FPPos)
        #print neiNum


        #sys.exit()
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
        print slope, FP, FN, FN+FP, rate1, rate2
        slope = slope + 0.05

    for v in bestCaseFPFN:
        print ("%.2f\t" % v),
