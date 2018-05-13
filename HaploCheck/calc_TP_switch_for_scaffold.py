#import deal_sam
import sys
import tools
import read

def read_align(filename):
    f = open(filename, "r")

    align = {}
    ref_pos = []
    contig_pos = []
    name = ""
    ref_label = 1 
    for line in f:
        if line.startswith("NC"):
            name = line.strip().split()[1]
            continue
        '''
        if not line.startswith("1") or not line.startswith("0"):
            continue
        '''
        words = line.strip().split(',')
        #print (words)

        if len(words) >1 and ref_label == 1:
            #print (words) 
            ref_label = 0   
            for c in words:
                ref_pos.append(int(c))
            continue  

        if len(words) >1 and ref_label == 0:
            for c in words:
                contig_pos.append(int(c))
        
    #print ("ref Pos", ref_pos)

    #print ("contig Pos", contig_pos)
    fout = open("contig_pos_ref_pos", "w")
    align[name] =[]   
    for i in range(len(contig_pos)):
        align[name].append((contig_pos[i], ref_pos[i]))
        fout.write("%s %s\n" %  (contig_pos[i], ref_pos[i]))
    fout.close()
    #print ("align", align)
    return align

def read_mutation_record(filename):
    f = open(filename, "r")

    pos = []
    #mutation = {}
    for line in f:  
        words = line.split()
        pos.append(int(words[0]))
        #mutation[pos] = (words[1], words[2])
    #print ("real snp pos:")
    #print (pos)
    return pos#, mutation

def read_snp_mutation(filename):
    snp_pos = []
    f = open(filename)
    name = ""
    for line in f:
        words = line.strip().split()
        if words[0] == "name:":
            name = words[1]
        elif len(words) == 5:
            snp_pos.append(int(words[0]))
    
    #print ("contig snp pos:")
    #print (snp_pos)
    return (name, snp_pos)

def convert(align, name, pos):
    assert name in align

    align_pairs = align[name]
    align_map = {} # key is contig_pos, value is ref_pos
    for (i, j) in align_pairs:
        if isinstance(i, int) and isinstance(j, int): 
            align_map[i] = j

    #print ("align pairs number:", len(align_pairs))
    #print ("both int number", len(align_map))
    result = []
    not_find = [] #alignment not report 
    for p in pos:
        if p in align_map:
            #print (p, align_map[p])
            result.append(align_map[p])
        else:
            not_find.append(p)

  
    return result, not_find 



def convert_reverse(align, name, pos):
    assert name in align
    align_pairs = align[name]
    align_map = {} # key is contig_pos, value is ref_pos
    for (i, j) in align_pairs:
        if isinstance(i, int) and isinstance(j, int): 
            align_map[j] = i

    #print ("align pairs number:", len(align_pairs))
    #print ("both int number", len(align_map))
    result = []
    not_find = [] #alignment not report 
    for p in pos:
        if p in align_map:
            #print (p, align_map[p])
            result.append(align_map[p])
        else:
            not_find.append(p)

    return result, not_find 




if __name__ == "__main__":

    

    #before run this script
    #step1: blasr: scaffold.fasta ref.fasta 
    
    if len(sys.argv) < 4:
        print ("python " + sys.argv[0] + "(1)/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record (2)/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/haplo/after_filter/contig_1_snp_mutation (3)scaffold_ref.blasr.m5")
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

    alignScaffoldRef._print_align_part(80,90)

    leftToRight, rightToLeft = alignScaffoldRef._generate_pos_to_pos()

    #convert contie pos to ref pos
    #print rightToLeft
    #print pre_snpPosition

    left_snpPosition, not_find = tools.convert(rightToLeft, pre_snpPosition) 
    print len(set(left_snpPosition).intersection(real_snpPosition))
    print len(not_find)


    right_snpPosition, not_find = tools.convert(leftToRight, real_snpPosition) 
    print len(set(right_snpPosition).intersection(pre_snpPosition))
    print len(not_find)


    print ("FP number:", len(set(left_snpPosition) - real_snpPosition))
    print ("TN number:", len(real_snpPosition - set(left_snpPosition)))
  
 
    #sys.exit() 

    '''
    TP_in_contig_from_real, nf = convert_reverse(align, name, real_snp_pos)

    TP_in_contig = set(contig_snp_pos).intersection(set(TP_in_contig_from_real))


    fout2 = open("TP_in_"+name, "w")
    print ("TP in contig number:", len(TP_in_contig))
    for cc in sorted(TP_in_contig):
        fout2.write("%s " % (cc))
    fout2.close()


    TP, nf = convert(align, name, TP_in_contig)

    print ("TP number:", len(TP))
    print ("not found number:", len(not_find))

    fout = open("TP_in_reference", "w")
    for cc in sorted(TP):
        fout.write("%s " % (cc))
    fout.close() 
    ''' 




