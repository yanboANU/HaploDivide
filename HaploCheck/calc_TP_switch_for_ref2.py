import sys
import read
import tools  

if __name__ == "__main__":

    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record
    #*snp_mutation
    #*phasing_result
    #reference length

    # function: calc TP 
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], int(sys.argv[4]))
    pre_snpPosition, pre_snpContent = read.read_snp2(sys.argv[2], int(sys.argv[4]))

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
  
    # function: calc switch number
    phasingHaplos = read.read_phasing_result(sys.argv[3])
    print ("phasingHaplos size", len(phasingHaplos))

    FPf =open("FP_content", "w")
    FPf.write("FP content\n")

    FPf.write("#FP is %d\n" % len(FP))
    count = 0

    for cc in FP:
        FPf.write("FP: %s" % cc)
        for phasingHaplo in phasingHaplos: 
            if cc in phasingHaplo:
                if phasingHaplo[cc] == '*':
                    count += 1
                FPf.write(" FP's label: %s" % phasingHaplo[cc])
        FPf.write("\n")    

    FPf.write("#* is %d" % count)    
    FPf.close()

    unPhasingNum = 0
    for phasingHaplo in phasingHaplos: 
        phasingHaploTP = ""
        print len(TP), len(phasingHaplo)
        for cc in TP:
            if cc not in phasingHaplo:
                #print cc, "not in"
                continue
            phasingHaploTP = phasingHaploTP + phasingHaplo[cc]

        refHaploTP = ""
        for cc in TP:
            if cc not in phasingHaplo:
                #print cc, "not in"
                continue
            #print "aa"
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
     
        count1, diffPos1 = tools.hamming_Distance(refHaploTP, phasingHaploTP)

        count2, diffPos2 = tools.hamming_Distance(tools.bool_Reverse(refHaploTP), phasingHaploTP)
        if count1 < count2:
            print "diff length:", count1
            for c in diffPos1:
                if phasingHaploTP[c] == '*':
                    unPhasingNum += 1
                else:    
                    print "diff", c, phasingHaploTP[c]
        else:
            print "diff length:", count2
            for c in diffPos2:
                if phasingHaploTP[c] == '*':
                    unPhasingNum += 1
                else:    
                    print "diff", c, phasingHaploTP[c]
    print "total unphasing number:", unPhasingNum          
