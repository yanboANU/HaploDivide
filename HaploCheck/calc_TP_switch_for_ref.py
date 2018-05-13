import sys
import read
import tools  

if __name__ == "__main__":

    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record
    #*snp_mutation
    #*phasing_result

    # function: calc TP 
    real_snpPosition, real_snpContent = read.read_snp(sys.argv[1])
    pre_snpPosition, pre_snpContent = read.read_snp(sys.argv[2])

    print ("number of real snp:", len(real_snpPosition))
    print ("number of pre snp:", len(pre_snpPosition))
    
    TP =  sorted(pre_snpPosition.intersection(real_snpPosition))

    print ("TP number:", len(TP))
    print ("FP number:", len(pre_snpPosition - real_snpPosition))
    print ("FP:", sorted(pre_snpPosition - real_snpPosition))
    print ("TN number:", len(real_snpPosition - pre_snpPosition))
    print ("TN:", sorted(real_snpPosition - pre_snpPosition))
  
    # function: calc switch number
    phasingHaplo = read.read_phasing_result(sys.argv[3])
 
    phasingHaploTP = "" 
    for cc in TP:
        if cc not in phasingHaplo:
            #print cc
            continue
        phasingHaploTP = phasingHaploTP + phasingHaplo[cc]

    refHaploTP = ""
    for cc in TP:
        if cc not in phasingHaplo:
            #print cc
            continue
        if (real_snpContent[cc][0] == pre_snpContent[cc][0]): 
        #or real_snpContent[cc][0] == reverseNucleotide(pre_snpContent[cc][1])):
            refHaploTP = refHaploTP + "0"
        elif (real_snpContent[cc][0] == pre_snpContent[cc][1]):
        #or real_snpContent[cc][0] == reverseNucleotide(pre_snpContent[cc][0])):
            refHaploTP = refHaploTP + "1"
        else:
            refHaploTP = refHaploTP + "2"
    print (len(TP),len(phasingHaploTP),len(refHaploTP))
    print (phasingHaploTP)
    print (refHaploTP)
 
    print (tools.hamming_Distance(refHaploTP, phasingHaploTP))
