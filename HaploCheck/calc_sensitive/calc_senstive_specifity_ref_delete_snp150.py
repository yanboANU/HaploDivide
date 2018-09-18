import sys
import read
import tools  
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import contig


def read_delete(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    starNumber = {}
    for line in f:  
        words = line.strip().split()
        if len(words) >= 5: #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            snpContent[pos] = (words[1], words[3])
            if words[1] == "*":
                starNumber[pos] = int(words[2])    
            if words[3] == "*":
                starNumber[pos] = int(words[4])
    
        if len(words) == 3: #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            snpContent[pos] = (words[1], words[2])
    return snpPosition, snpContent, starNumber

def write_delete_rate_distribution(FP, startNumber, cov, filename):
    
    FPRateCount = {}
    #fout = open(filename, "w")
    for c in FP:
        if c not in cov:
            print c, "not in cov"
            continue
        if c not in starNumber:
            print c, "not in starNumber"
            continue
        rate = round(starNumber[c]/float(cov[c]),2)
        if rate > 0.8:
            print(c, starNumber[c], cov[c], rate)
        if rate not in FPRateCount:
            FPRateCount[rate] = 0 
        FPRateCount[rate] += 1
    RC = sorted(FPRateCount.items())
    ''' 
    for (r,c) in RC:
        fout.write("%.2f %d\n" %(r, c))
    fout.close() 
    '''



def calc_TP_coverage(snpP1, snpP2, columns, minCov, contig):
    #print "coverage: ", minCov        
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


def remove_multi_delete(snpPos):
    ans = set()
    snpLen = len(snpPos)
    i = 0
    while i < snpLen-1:
        if snpPos[i] < snpPos[i+1]-1:
            if (i==0) or (i>=1 and snpPos[i-1] != snpPos[i]-1):
                #print i, snpPos[i]
                ans.add(snpPos[i])
        i += 1
    return ans    



if __name__ == "__main__":

    if len(sys.argv) < 5:
        print ("python " + sys.argv[0] + " real_snp pre_snp which_block  contig(file) *cov(file)")
        sys.exit()

    print ("python " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4] + " " + sys.argv[5])
    
    if sys.argv[3] == "0":
        start, end, base = 0, sys.maxint, 0
    
    ignoreLen = 10000
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base) #snp 138 record real delete position
    pre_snpPosition, pre_snpContent, starNumber = read_delete(sys.argv[2])
    print "delete number: ",len(pre_snpPosition)
    #print sorted(pre_snpPosition)[0:10]
    pre_snpPosition = remove_multi_delete(sorted(pre_snpPosition))
    print "after remove multiple delete, delete number: ", len(pre_snpPosition)
    #print sorted(pre_snpPosition)[0:10]
    #sys.exit()
    contigs = contig.read_Contig(sys.argv[4])
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
    
    # check and correction position for snp138
    #count = 0
    '''
    for key in real_snpContent:
        #print key
        #print real_snpContent[key]
        #print contig._seq[key-4:key+5]
        if real_snpContent[key][1] != contig._seq[key-1] and real_snpContent[key][1] != contig._seq[key-1].upper():
            print key
            print real_snpContent[key]
            print contig._seq[key-4:key+5]
        assert real_snpContent[key][1] == contig._seq[key-1] or real_snpContent[key][1] == contig._seq[key-1].upper()  
    print "snp 150 real snp Content accord with this reference"
    count = 0 
    for key in pre_snpContent:
        if pre_snpContent[key][0] != contig._seq[key] and pre_snpContent[key][0] != contig._seq[key].upper():
            #print key
            #print pre_snpContent[key]
            #print contig._seq[key-4:key+5] 
            count += 1
        #assert pre_snpContent[key][0] == contig._seq[key] or real_snpContent[key][1] == contig._seq[key].upper()  
    
    print "this individual has ", count, "position different with reference"
    '''
    FP, TN, TP = tools.print_TP(real_snpPosition, pre_snpPosition, contig._seq)

    print "FP length and TN length", len(FP), len(TN)
    
    print "first 10 FP ", FP[0:10]
    print "first 10 TN", TN[0:10]

    cov, count = read.read_cov(sys.argv[5])
      
    #print "write TP distribution"
    #write_delete_rate_distribution(TP, starNumber, cov, "TP_delete_rate_distribution")

    #print "write FP distribution"
    #write_delete_rate_distribution(FP, starNumber, cov, "FP_delete_rate_distribution")
    
    
    #print "write TN distribution"
    #write_delete_rate_distribution(TN, starNumber, cov, "TN_delete_rate_distribution")

  
    '''
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,0, contig)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,7, contig)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,10, contig)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,15, contig)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,20, contig)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,25, contig)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, cov,30, contig)
    '''
