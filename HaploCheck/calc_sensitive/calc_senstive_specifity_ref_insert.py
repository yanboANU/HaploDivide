import sys
import read
import tools  
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import contig


def read_insert(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    for line in f:  
        words = line.strip().split()

        if len(words) == 2: #for file mutation record
            snpPosition.add(int(words[0]))
            snpContent[int(words[0])] = (words[1])
            continue

        if len(words) == 3 and (not words[0].startswith("in")) : #for file mutation record
            snpPosition.add(int(words[0]))
            snpContent[int(words[0])] = (words[1])
    return snpPosition, snpContent


if __name__ == "__main__":


    if len(sys.argv) < 6:
        print ("python " + sys.argv[0] + " real_snp pre_snp which_block  contig(file) *cov(file) average_coverage")
        sys.exit()

    print ("python " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4] + " " + sys.argv[5] + " " +sys.argv[6])
  
    if sys.argv[3] == "0":
        start, end, base = 0, sys.maxint, 0
    if sys.argv[3] == "1":
        start, end, base = 585989, 2702781, -585988
    if sys.argv[3] == "2":
        start, end, base = 2746291, 12954384, -2746290 
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
    
    ignoreLen = 10000
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base) # real start with 1
    pre_snpPosition, pre_snpContent = read_insert(sys.argv[2])

    #print "real snp number:", len(real_snpPosition)
    #print "pre snp number:", len(pre_snpPosition)
    #sys.exit()
    print ("error_rate(%) 1-sensitive(%)")
    contigs = contig.read_Contig(sys.argv[4])
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
    
    # check and correction position
    '''
    count = 0
    for key in real_snpContent:
        print real_snpContent[key]
        print contig._seq[key-4:key+5]
        count += 1
        if count == 10:
            sys.exit()
    
    count = 0
    for key in pre_snpContent:
        print pre_snpContent[key]
        print contig._seq[key-4:key+5]
        count += 1
        if count == 10:
            sys.exit()
    ''' 
    
    
    
    FP,TN, TP = tools.print_TP(real_snpPosition, pre_snpPosition, contig._seq)
    '''    
    cov = read.read_cov(sys.argv[5]) 

    averageCov = int(sys.argv[6])
    snpP1F, snpP2F = tools.calc_TP_coverage(real_snpPosition, pre_snpPosition, cov, averageCov/2, contig._seq)
    #snpP1FF, snpP2FF = calc_TP_coverage(snpP1F, snpP2F, cov,30, contig)
     
    snpP1FF, snpP2FF = tools.calc_TP_coverage_between(snpP1F, snpP2F, cov, averageCov/2, averageCov*2, contig._seq)
    #snpP1F, snpP2F = calc_TP_coverage_between(real_snpPosition,pre_snpPosition, cov, 30, 90,  contig)
    '''
