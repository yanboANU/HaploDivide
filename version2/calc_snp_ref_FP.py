import sys

def read_mutation_list(filename):
    f = open(filename, "r")

    pos = []  
    for line in f:  
        words = line.split()
        if len(words) == 1 or len(words) == 3:
            pos.append(int(words[0]))
    
    return pos


  

if __name__ == "__main__":

    #filter_SAM(sys.argv[1], sys.argv[2])
    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/scaffols.ref.sam
    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/mutation_record
    #/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/mutation1/longest_yeast_mutation1/flye/haplo/after_filter/contig_1_snp_mutation
    #align = deal_sam.get_align_pos(sys.argv[1]) # bwa mem result
    # position in ref
    real_snp_pos = read_mutation_list(sys.argv[1])
    pre_snp_pos = read_mutation_list(sys.argv[2])


    print ("number of real snp:", len(real_snp_pos))
    print ("number of pre snp:", len(pre_snp_pos))
    '''
    #position in contig
    (name, contig_snp_position) = read_snp_mutation(sys.argv[3])  
 
    #convert contig pos to ref pos
    ref_snp_position, not_find = convert(align, name, contig_snp_position)
    '''
    #print (set(ref_snp_position).intersection(set(real_snp_position)))  
 
    print ("TP number:", len(set(pre_snp_pos).intersection(set(real_snp_pos))))

    print ("FP number:", len(set(pre_snp_pos) - set(real_snp_pos)))

    print ("FP:", sorted(set(pre_snp_pos) - set(real_snp_pos)))

    print ("TN number:", len(set(real_snp_pos)- set(pre_snp_pos)))

    print ("TN:", sorted(set(real_snp_pos)- set(pre_snp_pos)))



