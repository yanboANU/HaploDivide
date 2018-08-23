import sys
import read
import tools  
import contig

def read_insert(filename):

    f = open(filename, "r")
    snpPosition = set()
    #snpContent = {}
    starNumber = {}
    for line in f:  
        words = line.strip().split()
        if len(words) == 3 and (not words[0].startswith("in")) : #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            starNumber[pos] = int(words[2])    
    return snpPosition, starNumber
   


def write_insert_rate_distribution(FP, startNumber, cov, filename, averageCov):
    
    FPRateCount = {}
    fout = open(filename, "w")
    for c in FP:
        if c not in cov:
            print c, "not in cov"
            continue
        if c not in starNumber:
            print c, "not in starNumber"
            continue
        if cov[c] != averageCov:
            continue
        rate = round(starNumber[c]/float(cov[c]),2)
        if rate not in FPRateCount:
            FPRateCount[rate] = 0 
        FPRateCount[rate] += 1
    RC = sorted(FPRateCount.items())
    for (r,c) in RC:
        fout.write("%.2f %d\n" %(r, c))
    fout.close()    



if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print ("python " + sys.argv[0] + " insert_record *insert* which_block contig *_cov average_coverage")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3]+ " " + sys.argv[4] + " " + sys.argv[5] + " " +sys.argv[6])

    
    if sys.argv[3] == "10":
        start, end, base = 0, 10000000, 0 # first 10m
    if sys.argv[3] == "0":
        start, end, base = 0, sys.maxint, 0 # real and pre accord to each other now
    if sys.argv[3] == "1":
        start, end, base = 585989, 2702781, -585988
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
   
    ignoreLen = 10000
    real_insertPosition, real_insertContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base) # mutation_record
    insertPos, starNumber = read_insert(sys.argv[2]) # 

    contigs = contig.read_Contig(sys.argv[4])
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
    FP, TN, TP = tools.print_TP(real_insertPosition, insertPos, contig)

    cov = read.read_cov(sys.argv[5])
    averageCov = int(sys.argv[6])
    print "write TP distribution"
    write_insert_rate_distribution(TP, starNumber, cov, "TP_insert_rate_distribution_" +str(averageCov), averageCov)

    print "write FP distribution"
    write_insert_rate_distribution(FP, starNumber, cov, "FP_insert_rate_distribution_"+str(averageCov), averageCov)

