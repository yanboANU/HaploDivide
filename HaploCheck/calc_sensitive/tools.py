import sys


#nucleotide in string s are same, like AAAA/TTT return True
def is_same(s):
    c = s[0]
    sLen = len(s)
    for i in range(1,sLen):
        if c!=s[i]:
            return False
    return True    



def same_aA(a,b):
    return a==b or a==b.upper() or a.upper() == b

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
    #print "enter find neighbor"
    aa = []
    for c in a:
        aa.append((c,"TN")) 
    for c in b:
        aa.append((c,"FP"))

    bb = sorted(aa)    
    bbLen = len(bb)
    ans = 0
    for i in range(bbLen-1):
        #print bb[i]
        '''
        if bb[i][1] != bb[i+1][1] and (bb[i+1][0]-bb[i][0] <= 2):
            print "neigbhor"
            print (bb[i], bb[i+1])
            print (ref[bb[i][0]:bb[i+1][0]+1])
        '''    
        if bb[i][1] != bb[i+1][1] and (bb[i+1][0]-bb[i][0] <= 2) and is_same(ref[bb[i][0]:bb[i+1][0]+1]):
            '''
            if bb[i+1][0]-bb[i][0] == 2:
                print (bb[i], bb[i+1])
                print (ref[bb[i][0]:bb[i+1][0]+1])
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

# mutil delete means delete a sequence "ATGAC", "A"
def remove_multi_delete(snpPos):
    ans = set()
    snpLen = len(snpPos)
    i = 0
    # remove all multiple position
    '''
    while i < snpLen-1:
        if snpPos[i] < snpPos[i+1]-1:
            if (i==0) or (i>=1 and snpPos[i-1] != snpPos[i]-1):
                #print i, snpPos[i]
                ans.add(snpPos[i])
        i += 1
    '''
    # remove all multiple position expect first one
    while i < snpLen-1:
        #if (i==0) or (i>=1 and snpPos[i-1] != snpPos[i]-1):

        if (i==0) or (i>=1 and snpPos[i-1] < snpPos[i]-50):
            #print i, snpPos[i]
            ans.add(snpPos[i])
        i += 1
        
    return ans    


def merge_multi_delete(snpPos):
    snpStart = set()
    snpRange = []
    snpLen = len(snpPos)
    i = 0
    start = 0
    end = 0
    while i < snpLen-1:
        if (i==0) or (i>=1 and snpPos[i-1] < snpPos[i]-100):
            snpStart.add(snpPos[i])
            if start != 0 and end != 0:
                snpRange.append((start, end, end-start+1, count))
            start = snpPos[i]
            end = snpPos[i]
            count = 0
        else:    
            end = snpPos[i]    
        i += 1
        count += 1
        
    snpRange.append((start, end, end-start+1, count))
    #print "number of snpRange:", len(snpRange)    
    #print "numberStart: ", len(snpStart)
    ans = set()
    for (a,b,c, d) in snpRange:
        if c <= 1:
            ans.add(a)
        #if c > 10:
        #print a, b, c, d
    return ans    


def same_aA(a,b):
    return a==b or a==b.upper() or a.upper() == b


def print_TP(snpP1, snpP2, contig):

    realLen = len(snpP1)
    preLen = len(snpP2) 
    TP =  sorted(snpP1.intersection(snpP2)) 
    FP = sorted(snpP2 - snpP1)
    TN = sorted(snpP1 - snpP2)

    print ("FP length and TN length", len(FP), len(TN))
    if len(contig) > 0:
        neiNum, TN, FP = find_neighbor(TN, FP, contig)
    else:
        neiNum = 0
    print ("wrong TN near FP number ",neiNum)
    print ("remove TN near FP,FP length and TN length", len(FP), len(TN))

    print ("first 10 FP ", FP[0:10])
    print ("first 10 TN", TN[0:10])

    if preLen == 0:
        return FP, TN
    TPLen = len(TP)
    FPLen = len(FP)
    TNLen = len(TN)
    TPLen = TPLen + neiNum

    TPRate  = float(TPLen)/preLen*100 
    rate1 = round(FPLen/float(TPLen+FPLen) * 100, 1)
    rate2 = round(TNLen/float(realLen) * 100, 1)

    #print "FPLen  TNLen FPLen+TNLen TPLen realLen preLen ErrorRate 1-sensitive"
    #print ("%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f"% (FPLen, TNLen, FPLen + TNLen, TPLen, realLen, preLen, rate1, rate2))
    '''
    print "FPLen  TNLen FPLen+TNLen TPLen realLen preLen ErrorRate 1-sensitive"
    print ("%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f"% (FPLen, TNLen, FPLen + TNLen, TPLen, realLen, preLen, rate1, rate2))
    '''

    print ("preLen TPLen FPLen TNLen ErrorRate 1-sensitive realLen FPLen+TNLen ")
    print ("%d\t%d\t%d\t%d\t %.1f\t%.1f\t %d\t%d"% (preLen, TPLen, FPLen, TNLen,  rate1, rate2, realLen, FPLen+TNLen))
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
        if p not in columns:
            continue

        if columns[p] >= minCov:
            snpP2F.add(p)
    print (">=%s" % minCov) 
    print_TP(snpP1F, snpP2F, contig)
    return snpP1F, snpP2F


def calc_TP_coverage_ratio(realSnp, preSnp, preSupportNum, threshold, contig): 
    
    preSnp2 = set()
    #threshold = 2
    for p in preSnp:
        h1Cov, h2Cov = preSupportNum[p]
        if max(h1Cov, h2Cov) < threshold*min(h1Cov, h2Cov):
            preSnp2.add(p)
    print_TP(realSnp, preSnp2, contig)
    return realSnp, preSnp2



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
