import sys
import read
import tools  
import contig



def read_snp(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    for line in f:  
        words = line.strip().split()
        if len(words) == 2: #for file mutation record
            snpPosition.add(int(words[0]))
            snpContent[int(words[0])] = (words[1])

        if len(words) == 3: #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            snpContent[pos] = (words[1], words[2])
    return snpPosition, snpContent

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print ("python " + sys.argv[0] + " pre_snp contig(file)")
        sys.exit()

    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2])

    pre_snpPosition, pre_snpContent = read_snp(sys.argv[1])

    contigs = contig.read_Contig(sys.argv[2])
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
    seq = contig._seq 


    positionS = sorted(pre_snpPosition) 
    filename = sys.argv[1] + "_upper"
    fout = open(filename, "w")
    for p in positionS:
        if seq[p].isupper():
            fout.write("%s %s\n" % (p, pre_snpContent[p]))
    fout.close()


    filename = sys.argv[1] + "_lower"
    fout = open(filename, "w")
    for p in positionS:
        if seq[p].islower():
            fout.write("%s %s\n" % (p, pre_snpContent[p]))
    fout.close()        
