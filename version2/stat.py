import sys
import read


def stat_pre_snp_cov(preMfile, columnsfile):

    pre_mutation = read.read_mutation_list(preMfile)
    #record in reference pos
    columns = read.read_columns(columnsfile)


    fout = open("pre_snp_first_cov", "w")    
 
    fout2 = open("pre_snp_second_cov", "w")     
    for pos in pre_mutation:
        fout.write("%s\n" % (columns[pos]._diff_cov[0]))
        #fout.write("%s\n" % (columns[pos]._diff_cov[0]/columns[pos]._cov))
        columns[pos]._print()
        if pos not in columns:
            break    
        if len(columns[pos]._diff_cov) == 1:
            fout2.write("0\n")
        else:
            fout2.write("%s\n" % (columns[pos]._diff_cov[1]))
            #fout2.write("%s\n" % (columns[pos]._diff_cov[1]/columns[pos]._cov))

    fout.close()
    fout2.close() 

def stat_real_snp_cov(realMfile, columnsfile, alignfile):

    real_mutation = read.read_mutation_list(realMfile)
    #record in reference pos
    columns = read.read_columns(columnsfile)

    ref_pos, contig_pos = read.read_align(alignfile) # blasr result
    real_mutation_in_contig = []
    #print ("aaa")
    align_map = {}
    for i in range(len(ref_pos)):
        if isinstance(ref_pos[i], int) and isinstance(contig_pos[i], int): 
            align_map[ref_pos[i]]= contig_pos[i]


    for pos in real_mutation:
        # some position of reference donot have alignment
        if pos not in align_map:
            continue     
        real_mutation_in_contig.append(align_map[pos])
        #print (pos, align_map[pos])
    #cover reference pos to contig pos
    fout = open("snp_first_cov_ratio", "w")    
 
    fout2 = open("snp_second_cov_ratio", "w")     
    for pos in real_mutation_in_contig:
        #fout.write("%s\n" % (columns[pos]._diff_cov[0]))
        fout.write("%s\n" % (columns[pos]._diff_cov[0]/columns[pos]._cov))
        columns[pos]._print()
        if pos not in columns:
            break    
        if len(columns[pos]._diff_cov) == 1:
            #fout2.write("0\n")
            fout2.write("0\n")

        else:
            #fout2.write("%s\n" % (columns[pos]._diff_cov[1]))
            fout2.write("%s\n" % (columns[pos]._diff_cov[1]/columns[pos]._cov))
    fout.close()
    fout2.close() 


def stat_distance_snp(snpMfile,snpDisfile):

    snp_pos = read.read_mutation_list(snpMfile)
    #record in reference pos

    fout = open(snpDisfile, "w")    
    
    i=1 
    while i < len(snp_pos):
        fout.write("%s\n" % (snp_pos[i] - snp_pos[i-1]))
        i += 1
    fout.close()

if __name__ == '__main__':

    #stat_real_snp_cov(sys.argv[1], sys.argv[2], sys.argv[3])
    #stat_pre_snp_cov(sys.argv[1], sys.argv[2])
    stat_distance_snp(sys.argv[1], sys.argv[2])
