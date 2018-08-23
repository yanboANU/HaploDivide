import os
import sys
import contig
def read_1st(filename):
    f = open(filename, "r")
    val = {}
    for line in f:
        words = line.split()
        wordsLen = len(words) 
        assert wordsLen == 2
        rate = float(words[0])
        val[rate] = int(words[1])
    return val;




if __name__ == "__main__":
    
    snp1st = read_1st(sys.argv[1])
    nonSNP1st = read_1st(sys.argv[2])

    slope = 0.1
    print ("slope, FP, FN, FN+FP, totalNode, FP/float(TP+FP), float(FN)/(realNode)")
    minFPFN = 50000
    bestCaseFPFN = []
    while slope <= 1:
        FN = 0   # Snp, think it is not
        FP = 0   # nonSnp, think it is
        TP = 0
        realNode = 0
        slope = round(slope,2)
        for x in snp1st:
            realNode = realNode + snp1st[x] 
            if x < slope:
                #print (x, slope, snp1st[x])
                FN = FN + snp1st[x]
            else:
                TP = TP + snp1st[x]
        for x in nonSNP1st:
            if x >= slope:
                FP = FP + nonSNP1st[x]
        if (TP+FP) == 0:
            break    
        #print (FN, FP)
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
        print ("%s %s %s %s %.2f %.2f"%(slope, FP, FN, FN+FP, rate1, rate2))
        slope = slope + 0.01
    print ("slope FP FN error_rate 1-sensitive")
    for v in bestCaseFPFN:
        print ("%.2f\t" % v, end='')
    print ("\n")    
