import deal_sam
import sys
#import ../calc_swith

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

    #filter_SAM(sys.argv[1], sys.argv[2])
    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/scaffols.ref.sam
    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record
    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/haplo/after_filter/contig_1_snp_mutation
    #align = deal_sam.get_align_pos(sys.argv[1]) # bwa mem result
    

    #before run this script
    #step1: blasr:  scaffold.fasta ref.fasta 
    #step2: input: blasr result, output:  align_pos(calc_switch)
    if len(sys.argv) < 4:
        print ("python " + sys.argv[0] + " /media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/align_pos  /media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record /media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/haplo/after_filter/contig_1_snp_mutation")
        sys.exit()

    align = read_align(sys.argv[1]) # blasr result
    # pos in ref
    real_snp_pos = read_mutation_record(sys.argv[2])

    #pos in contig
    (name, contig_snp_pos) = read_snp_mutation(sys.argv[3])  
 
    #convert contie pos to ref pos
    ref_snp_pos, not_find = convert(align, name, contig_snp_pos)
    
    #print (set(ref_snp_pos).intersection(set(real_snp_position)))  
 
    print ("number of real snp:", len(real_snp_pos))
    print ("number of pre snp:", len(ref_snp_pos))

    '''
    TP  =  set(ref_snp_pos).intersection(set(real_snp_pos))

    print ("TP number:", len(TP))
    print ("not found number:", len(not_find))

    fout = open("TP_in_reference", "w")
    for cc in sorted(TP):
        fout.write("%s " % (cc))
    fout.close()

    #TP_in_contig, nf = convert_reverse(align, name, TP)
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


    print ("FP number:", len(set(ref_snp_pos) - set(real_snp_pos)))

    #print ("FP:", sorted(set(ref_snp_pos) - set(real_snp_pos)))

    print ("TN number:", len(set(real_snp_pos)- set(ref_snp_pos)))

    #print ("TN:", sorted(set(real_snp_pos)- set(ref_snp_pos)))


