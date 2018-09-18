import sys
import read
import tools  
import contig


if __name__ == "__main__":


    if len(sys.argv) < 5:
        print ("python " + sys.argv[0] + " real_snp pre_snp which_block *cov(file)")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4])

    if sys.argv[3] == "0":
        #start, end, base = 0, sys.maxint, -1

        start, end, base = 0, sys.maxint, 0 # pre and real accord to each other
    if sys.argv[3] == "1":
        #start, end, base = 585989, 2702781, -585989
        start, end, base = 585989, 2702781, -585988 # after 16 Aug, for new version main
    if sys.argv[3] == "2":
        start, end, base = 2746291, 12954384, -2746290 
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
    
    ignoreLen = 10000
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base)
    pre_snpPosition, pre_snpContent = read.read_snp(sys.argv[2])
 
    print "error_rate(%) 1-sensitive(%)"
    # function: calc switch number
    FP,TN,TP = tools.print_TP(real_snpPosition, pre_snpPosition, "")
   
    contigs = contig.read_Contig(sys.argv[4]) # contig seq start with 0
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()

    cov, count = read.read_cov(sys.argv[5])
    print "length larger than 8 ", count 
    #averageCov = int(sys.argv[5])
   

    snpP1F, snpP2F = tools.calc_TP_coverage(real_snpPosition, pre_snpPosition, cov, 8, contig._seq)
