#import deal_sam
import sys
import tools
import read

def reverseNucleotide(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'  



if __name__ == "__main__":

    

    #before run this script
    #step1: blasr: scaffold.fasta ref.fasta 
    
    if len(sys.argv) < 4:
        print ("python " + sys.argv[0] + "(1)/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record (2)/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/haplo/after_filter/contig_1_snp_mutation (3)scaffold_ref.blasr.m5 (4)*phasing_result")
        sys.exit()

    # pos in ref
    real_snpPosition, real_snpContent = read.read_snp(sys.argv[1])
    
    #pos in contig
    pre_snpPosition, pre_snpContent = read.read_snp(sys.argv[2])

    alignScaffoldRef = read.read_blasr_m5(sys.argv[3])

    print ("number of real snp:", len(real_snpPosition))
    print ("number of pre snp:", len(pre_snpPosition))
    
    alignScaffoldRef._left._generate_seq_pos() # real
    alignScaffoldRef._right._generate_seq_pos()# pre    

    alignScaffoldRef._print_align_part(1530113,90)

    leftToRight, rightToLeft = alignScaffoldRef._generate_pos_to_pos()

    #convert contie pos to ref pos
    #print rightToLeft
    #print pre_snpPosition

    left_snpPosition, not_find = tools.convert(rightToLeft, pre_snpPosition)
    TPInRef = set(left_snpPosition).intersection(real_snpPosition)
    print len(TPInRef)
    print len(not_find)


    right_snpPosition, not_find = tools.convert(leftToRight, real_snpPosition) 
    TPInContig = sorted(set(right_snpPosition).intersection(pre_snpPosition))
    print len(TPInContig)
    print len(not_find)




    print ("FP number:", len(set(left_snpPosition) - real_snpPosition))
    print ("TN number:", len(real_snpPosition - set(left_snpPosition)))
  
 
    #sys.exit() 


    Haplos = read.read_phasing_result(sys.argv[4])
    for phasingHaplo in Haplos: 
        phasingHaploTP = ""  
        refHaploTP = ""
        refHaploTP2 = "" 
        for cc in TPInContig:
            '''
            print cc
            print phasingHaplo
            print rightToLeft
            sys.exit()    
            '''
            if (cc not in phasingHaplo) or (cc not in rightToLeft):
                #print cc
                continue
            phasingHaploTP = phasingHaploTP + phasingHaplo[cc]
            
            print ("posision in contig %s position in reference %s" % (cc, rightToLeft[cc]))
            print ("content in contig " , pre_snpContent[cc])
            print ("content in reference " , real_snpContent[ rightToLeft[cc] ])
            ''' 
            if pre_snpContent[cc] != real_snpContent[ rightToLeft[cc] ]:
                print ("contig position %s content %s" % ( cc, pre_snpContent[cc]))
                print ("reference position %s content %s" % (rightToLeft[cc], real_snpContent[ rightToLeft[cc] ]))
            else:
                print ("equal")
                print ("contig position %s content %s" % ( cc, pre_snpContent[cc]))
                print ("reference position %s content %s" % (rightToLeft[cc], real_snpContent[ rightToLeft[cc] ]))
            '''
            if (real_snpContent[ rightToLeft[cc] ][0] == pre_snpContent[cc][0]):  
                #or real_snpContent[ rightToLeft[cc] ][0] == reverseNucleotide(pre_snpContent[cc][1])):
                refHaploTP = refHaploTP + "0"
            elif (real_snpContent[ rightToLeft[cc] ][0] == pre_snpContent[cc][1]):
                #or real_snpContent[ rightToLeft[cc] ][0] == reverseNucleotide(pre_snpContent[cc][0])):
                refHaploTP = refHaploTP + "1"
            else:
                refHaploTP = refHaploTP + "2"


            if (real_snpContent[ rightToLeft[cc] ][0] == reverseNucleotide(pre_snpContent[cc][1])):
                refHaploTP2 = refHaploTP2 + "0"
            elif (real_snpContent[ rightToLeft[cc] ][0] == reverseNucleotide(pre_snpContent[cc][0])):
                refHaploTP2 = refHaploTP2 + "1"
            else:
                refHaploTP2 = refHaploTP2 + "2"

        print (len(TPInContig),len(phasingHaploTP),len(refHaploTP))
        print (phasingHaploTP)
        print (refHaploTP)
        print (refHaploTP2)

        print (tools.hamming_Distance(refHaploTP, phasingHaploTP))

        print (tools.hamming_Distance(refHaploTP2, phasingHaploTP))



