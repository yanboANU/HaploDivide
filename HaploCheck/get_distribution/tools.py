import sys


#nucleotide in string s are same, like AAAA/TTT return True
def is_same(s):
    c = s[0]
    sLen = len(s)
    for i in range(1,sLen):
        if c!=s[i]:
            return False
    return True    

def same_len(s):
    c = s[0]
    sLen = len(s)
    same_len = 1
    for i in range(1,sLen):
        if c!=s[i]:
            return same_len
        else:
            same_len += 1
    return same_len        

#A A   => A
#1 2
#delte in position 1 equal to delete in position 2
def find_neighbor(a, b, ref):
    aa = []
    for c in a:
        aa.append((c,"TN"))
    for c in b:
        aa.append((c,"FP"))

    bb = sorted(aa)    
    bbLen = len(bb)
    ans = 0
    for i in range(bbLen-1):
        if bb[i][1] != bb[i+1][1] and (bb[i+1][0]-bb[i][0] <= 2) and is_same(ref._seq[bb[i][0]:bb[i+1][0]+1]):
            '''
            if bb[i+1][0]-bb[i][0] <= 2:
                print (bb[i], bb[i+1])
                print (ref._seq[bb[i][0]:bb[i+1][0]+1])
            '''    
            if bb[i][1] == "TN":
                if (bb[i][0] in a) and (bb[i+1][0] in b):
                    a.remove(bb[i][0])
                    b.remove(bb[i+1][0])
                    ans += 1
            else:
                if (bb[i][0] in a) and (bb[i+1][0] in b):
                    b.remove(bb[i][0])
                    a.remove(bb[i+1][0])
                    ans += 1
    return ans, a, b



def print_TP(snpP1, snpP2, contig):

    realLen = len(snpP1)
    preLen = len(snpP2) 
    TP =  sorted(snpP1.intersection(snpP2)) 
    FP = sorted(snpP2 - snpP1)
    TN = sorted(snpP1 - snpP2)

    print "FP length and TN length", len(FP), len(TN)
    neiNum, TN, FP = find_neighbor(TN, FP, contig)
    print "wrong TN near FP number ",neiNum
    print "remove TN near FP,FP length and TN length", len(FP), len(TN)

    print "first 10 FP ", FP[0:10]
    print "first 10 TN", TN[0:10]

    if preLen == 0:
        return FP, TN
    TPLen = len(TP)
    FPLen = len(FP)
    TNLen = len(TN)
    TPLen = TPLen + neiNum

    TPRate  = float(TPLen)/preLen*100 
    rate1 = round(FPLen/float(TPLen+FPLen) * 100, 1)
    rate2 = round(TNLen/float(realLen) * 100, 1)

    print "FPLen  TNLen FPLen+TNLen TPLen realLen preLen ErrorRate 1-sensitive"
    print ("%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f"% (FPLen, TNLen, FPLen + TNLen, TPLen, realLen, preLen, rate1, rate2))

    return FP, TN, TP


def calc_TP_coverage(snpP1, snpP2, columns, minCov, contig): 
    snpP1F = set()
    snpP2F = set()
    for p in snpP1:
        if p not in columns:
            continue
        if columns[p] >= minCov:
            snpP1F.add(p)

    for p in snpP2:
        if columns[p] >= minCov:
            snpP2F.add(p)
    print (">=%s" % minCov) 
    print_TP(snpP1F, snpP2F, contig)
    return snpP1F, snpP2F

def calc_TP_coverage_equal(snpP1, snpP2, columns, cov, contig): 
    snpP1F = set()
    snpP2F = set()
    for p in snpP1:
        if p not in columns:
            continue
        if columns[p] == cov:
            snpP1F.add(p)
    for p in snpP2:
        if columns[p] == cov:
            snpP2F.add(p)
    print ("=%s" % cov) 
    print_TP(snpP1F, snpP2F, contig)

def calc_TP_coverage_between(snpP1, snpP2, columns, minCov, maxCov, contig): 
    snpP1F = set()
    snpP2F = set()
    for p in snpP1:
        if p not in columns:
            continue
        if columns[p] >= minCov and columns[p] <= maxCov:
            snpP1F.add(p)

    for p in snpP2:
        if columns[p] >= minCov and columns[p] <= maxCov:
            snpP2F.add(p)
    print (">=%s <=%s" % (minCov, maxCov)) 
    print_TP(snpP1F, snpP2F, contig)
    return snpP1F, snpP2F


#pre_snpPosition - multiplePos
def filter_in_mutiple_range(pre_snpPosition, multiplePos):
   pre = set()
   for p in pre_snpPosition:
       if p not in multiplePos:
          pre.add(p)

   return pre       
