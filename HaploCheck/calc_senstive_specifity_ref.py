import sys
import read
import tools  

def read_columns(filename):
    columns= {}
    f = open(filename, "r")
    for line in f:
        words = line.strip().split(" ")
        if words[0].startswith("reference"):
            pos = int(words[2])
            coverage = int(words[5])
            columns[pos] = coverage
    '''        
    for key in columns:
        columns[key] = sorted(columns[key], reverse=True)
    '''
    return columns
def calc_TP_coverage(snpP1, snpP2, columns, minCov):
    
    #print "coverage: ", minCov        
    snpP1F = set()
    snpP2F = set()
    for p in snpP1:
        if p not in columns:
            #print p
            continue
        if columns[p] >= minCov:
            snpP1F.add(p)

    for p in snpP2:
        if columns[p] >= minCov:
            snpP2F.add(p)
    print ">=", minCov, 
    print_TP(snpP1F, snpP2F)

def print_TP(snpP1, snpP2):

    #print ("number of real snp:", len(snpP1))
    #print ("number of pre snp:", len(snpP2)) 
    TP =  sorted(snpP1.intersection(snpP2))
    
    FP = sorted(snpP2 - snpP1)
    #print ("FP:", len(FP))

    TN = sorted(snpP1 - snpP2)
    #print ("TN:", len(TN))
    
    print ("%.2f\t%.2f" % (100*(1-len(TP)/float(len(snpP1))),100*(1- len(TP)/float(len(snpP2))) ) )
    return FP, TN

if __name__ == "__main__":

    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record
    #*snp_mutation
    #*columns
    if sys.argv[3] == "0":
        start, end, base = 0, sys.maxint, 0
    if sys.argv[3] == "1":
        start, end, base = 585989, 2702781, -585989
    if sys.argv[3] == "2":
        start, end, base = 2746291, 12954384, -2746291 
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553836
    
    ignoreLen = 10000
    # function: calc TP 
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base)
    pre_snpPosition, pre_snpContent = read.read_snp(sys.argv[2])
 
    print "1-sensitive(%) error_rate(%)"
    # function: calc switch number
    FP,TN = print_TP(real_snpPosition, pre_snpPosition)

    #columns = read_columns(sys.argv[3])
    '''
    print "FP"
    for c in FP:
        print c, columns[c]
        
    print "TN"
    for c in TN:
        print c, columns[c]
    '''
    '''
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,0)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,5)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,10)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,15)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,20)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,25)
    calc_TP_coverage(real_snpPosition,pre_snpPosition, columns,30)
    '''
