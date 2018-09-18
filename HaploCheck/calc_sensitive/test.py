import sys
import read
import tools  
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import contig
import calc_senstive_specifity_ref_delete


if __name__ == "__main__":

    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 0, sys.maxint, 0)
    print sorted(real_snpPosition)
   
    calc_senstive_specifity_ref_delete.extend_real_snpPosition(real_snpPosition, real_snpContent)
    print sorted(real_snpPosition)
