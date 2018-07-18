import os
import sys
import read

def read_1st2nd(filename):
    f = open(filename, "r")
    val = {}
    for line in f:
        words = line.split()
        val[(float(words[0]), float(words[1]))] = int(words[2])
    return val;

def read_2file(filename1, filename2):

    f = open(filename1, "r")
    val = {}
    for line in f:
        words = line.split()
        val[(float(words[0]), float(words[1]))] = int(words[2])
    
    f = open(filename2, "r")
    for line in f:
        words = line.split()
        key = (float(words[0]), float(words[1]))  
        if key in val:
            val[key] = val[key] + int(words[2])
        else:
            val[key] = int(words[2])
    
    return val;

if __name__ == "__main__":
     
    snp1st2nd = read_1st2nd(sys.argv[1])
    nonSNP1st2nd = read_1st2nd(sys.argv[2])
    
    #realSNP = read.read

    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[3], 29553836, 121757928 , -29553836) # mutation_record
    #print "real snp number", len(real_snpPosition)
    #snp1st2nd = read_2file(sys.argv[1], sys.argv[3])
    #nonSNP1st2nd = read_2file(sys.argv[2], sys.argv[4])

    slope = 1.0

    #print "slope, FN, FP, FN+FP, totalNode, float(FN+FP)/(totalNode)"
    minFPFN = 50000
    minRate = 2
    bestCaseFPFN = []
    bestCaseRate = []
    while slope <= 3.5:
        FN = 0   # Snp, think it is not
        FP = 0   # nonSnp, think it is
        totalNode = 0
        TP = 0
        realNode = 0
        for (x,y) in snp1st2nd:
            #print x,y
            totalNode = totalNode + snp1st2nd[(x,y)] 
            realNode  = realNode + snp1st2nd[(x,y)] 
            if x > slope*y:
                #print "FN", x, y
                FN = FN + snp1st2nd[(x,y)]
            else:
                TP = TP + snp1st2nd[(x,y)]


        for (x,y) in nonSNP1st2nd:

            totalNode = totalNode + nonSNP1st2nd[(x,y)] 
            if x <= slope*y:
                #print "FP", x, y
                FP = FP + nonSNP1st2nd[(x,y)]
                
        rate1 = FP/float(TP+FP) 
        rate2 = FN/float(realNode)
        print "real Node number", realNode
        print slope,  FP+FN, rate1, rate2, rate1+rate2
        if FP+FN < minFPFN:
            bestCaseFPFN = []
            minFPFN = FP+FN

            bestCaseFPFN.append(slope)
            bestCaseFPFN.append(FP)
            bestCaseFPFN.append(FN)
            bestCaseFPFN.append(rate1*100)
            bestCaseFPFN.append(rate2*100)
        if rate1+rate2 < minRate:
            bestCaseRate = []
            minRate = rate1+rate2
            bestCaseRate.append(slope)
            bestCaseRate.append(FP)
            bestCaseRate.append(FN)
            bestCaseRate.append(rate1*100)
            bestCaseRate.append(rate2*100)
 
        slope = slope + 0.1

    for v in bestCaseFPFN:
        print ("%.1f\t" % v),
    
    for v in bestCaseRate:
        print ("%.1f\t" % v),
    print ("\n")
