import sys
import read
import tools  
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import contig


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
    pre_snpPosition, pre_snpContent = read.read_snp(sys.argv[2])  
    #reference start with 1, nomal record real delete position - 1, nucleotide is delete one
    
    contigs = contig.read_Contig(sys.argv[4]) # contig seq start with 0
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
   
    count = 0
    #sameLenDis = {}
    for key in real_snpContent:
        #if real_snpContent[key][1] != contig._seq[key-1] and real_snpContent[key][1] != contig._seq[key-1].upper():
        #print key
        #print real_snpContent[key]
        #print contig._seq[key-1:key+5]
        '''
        l = tools.same_len(contig._seq[key:key+10])
        if l not in sameLenDis:
            sameLenDis[l] = 0
        sameLenDis[l] += 1 
        '''
        assert real_snpContent[key][1] == contig._seq[key-1] or real_snpContent[key][1] == contig._seq[key-1].upper() 
    
    #print sameLenDis
    print "snp 150 real snp Content accord with this reference"
    count = 0

    #sameLenDis = {}
    for key in pre_snpContent:
        if pre_snpContent[key][0] != "*":
            delete_content = pre_snpContent[key][0]  
        else:
            delete_content = pre_snpContent[key][1]
        ''' 
        l = tools.same_len(contig._seq[key:key+10])
        if l not in sameLenDis:
            sameLenDis[l] = 0
        sameLenDis[l] += 1 
        '''
        if delete_content != contig._seq[key] and delete_content != contig._seq[key].upper():   
        #assert delete_content == contig._seq[key] or delete_content == contig._seq[key].upper()  
            count += 1
    #print sameLenDis        
    print count, "not accord with this reference"
    print "total snp: ", len(pre_snpContent)

    FP,TN, TP = tools.print_TP(real_snpPosition, pre_snpPosition, contig._seq)
    
    '''
    cov = read.read_cov(sys.argv[5])
    averageCov = int(sys.argv[6])
    snpP1F, snpP2F = tools.calc_TP_coverage(real_snpPosition, pre_snpPosition, cov, averageCov/2, contig._seq)
    snpP1FF, snpP2FF = tools.calc_TP_coverage(snpP1F, snpP2F, cov, averageCov/2+10, contig._seq)
    snpP1FF, snpP2FF = tools.calc_TP_coverage(snpP1F, snpP2F, cov, averageCov/2+20, contig._seq)
    snpP1FF, snpP2FF = tools.calc_TP_coverage_between(snpP1F, snpP2F, cov, averageCov/2, averageCov*2, contig._seq)
    '''
    
