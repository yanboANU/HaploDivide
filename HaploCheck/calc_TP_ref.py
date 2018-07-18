import sys
import read
import tools  

if __name__ == "__main__":

    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record
    #*snp_mutation
    #*phasing_result

    # function: calc TP
    start = 585989
    end = 2702781
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start + 10000, end-10000, 0)
    pre_snpPosition, pre_snpContent = read.read_snp2(sys.argv[2], 10000, (end-start+1)-10000, start)

    print ("number of real snp:", len(real_snpPosition))
    print ("number of pre snp:", len(pre_snpPosition))
    
    TP =  sorted(pre_snpPosition.intersection(real_snpPosition))

    print ("TP number:", len(TP))
    print ("FP number:", len(pre_snpPosition - real_snpPosition))
    print ("Sensitive: %.2f" % (float(len(TP))/len(real_snpPosition))) 
    print ("Error rate: %.2f" % (1- float(len(TP))/len(pre_snpPosition)))
    
    FP = sorted(pre_snpPosition - real_snpPosition)
    print ("FP:", len(FP))
    print ("TN number:", len(real_snpPosition - pre_snpPosition))
    print ("TN:", sorted(real_snpPosition - pre_snpPosition))
  
