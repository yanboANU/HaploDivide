import sys
import read
import tools  
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import contig


def extend_real_snpPosition(pos, content):

    for key in content:
        assert key in pos
        deleteLen = len(content[key][0]) - len(content[key][1])
        assert deleteLen >= 1
        tempKey = key
        while deleteLen > 1:
            tempKey = tempKey + 1
            pos.add(tempKey)
            deleteLen -= 1
    #return     


def check_consistence(real_snpContent, pre_snpContent, contig):
    
    for key in real_snpContent:
        assert real_snpContent[key][1] == contig._seq[key-1] or real_snpContent[key][1] == contig._seq[key-1].upper() 
    
    print "snp 150 real snp Content accord with this reference"
    count = 0
    for key in pre_snpContent:
        if pre_snpContent[key][0] != "*":
            delete_content = pre_snpContent[key][0]  
        else:
            delete_content = pre_snpContent[key][1]
        if delete_content != contig._seq[key] and delete_content != contig._seq[key].upper():   
            count += 1
            print key, delete_content, contig._seq[key:key+2]  
    print count, "not accord with this reference"
    print "total snp: ", len(pre_snpContent)


if __name__ == "__main__":

    if len(sys.argv) < 6:
        print ("python " + sys.argv[0] + " real_snp pre_snp which_block  contig(file) *cov(file) average_coverage")
        sys.exit()

    print ("python " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4] + " " + sys.argv[5] + " " + sys.argv[6])
    
    if sys.argv[3] == "0":
        start, end, base = 0, sys.maxint, 0 # real and pre accord to each other now
    if sys.argv[3] == "1":
        #start, end, base = 585989, 2702781, -585989
        start, end, base = 585989, 2702781, -585988 # after 16 Aug, for new version main
    if sys.argv[3] == "2":
        start, end, base = 2746291, 12954384, -2746291 
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
    if sys.argv[3] == "138":
        start, end, base = 0, sys.maxint, -2
    
    ignoreLen = 10000
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base)
    #snp 138 record real delete position 
    #reference start with 1, nomal record real delete position - 1
    pre_snpPosition, pre_snpContent, pre_snpSupportNum = read.read_pre_snp(sys.argv[2])  
    #reference start with 1, nomal record real delete position - 1, nucleotide is delete one
    #pre_snpPosition = tools.remove_multi_delete(sorted(pre_snpPosition))
    #pre_snpPosition = tools.merge_multi_delete(sorted(pre_snpPosition))
    #print "after remove multiple delete, delete number: ", len(pre_snpPosition)
    #sys.exit()
    contigs = contig.read_Contig(sys.argv[4]) # contig seq start with 0
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
  

    check_consistence(real_snpContent, pre_snpContent, contig)
    print ("real length", len(real_snpPosition))
    extend_real_snpPosition(real_snpPosition, real_snpContent)
    print ("extend real length", len(real_snpPosition))
     
    FP,TN, TP = tools.print_TP(real_snpPosition, pre_snpPosition, contig._seq)
    print "\n" 
     
    cov, count = read.read_cov(sys.argv[5])
    averageCov = int(sys.argv[6])
    
    snpP1F, snpP2F = tools.calc_TP_coverage(real_snpPosition, pre_snpPosition, cov, 8, contig._seq)
    
    print "\n" 
    #snpP1FF, snpP2FF = tools.calc_TP_coverage(snpP1F, snpP2F, cov, averageCov/2+10, contig._seq)
    #snpP1FF, snpP2FF = tools.calc_TP_coverage(snpP1F, snpP2F, cov, averageCov/2+20, contig._seq)
    snpP1FF, snpP2FF = tools.calc_TP_coverage_between(snpP1F, snpP2F, cov, averageCov/2, averageCov*2, contig._seq)
    
    print "\n" 
    tools.calc_TP_coverage_ratio(snpP1FF, snpP2FF, pre_snpSupportNum, 2, contig._seq)
    
    print "\n" 
    tools.calc_TP_coverage_ratio(snpP1FF, snpP2FF, pre_snpSupportNum, 1.5, contig._seq)
    
    print "\n" 
    tools.calc_TP_coverage_ratio(snpP1FF, snpP2FF, pre_snpSupportNum, 2.5, contig._seq)
    #snpP1FF, snpP2FF = tools.calc_TP_coverage_between(snpP1F, snpP2F, cov, averageCov/2+10, averageCov*2-10, contig._seq)
    
