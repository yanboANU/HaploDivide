import sys
import tools


#input: TP_in_reference TP_in_contig
#       mutation_record contig_snp_mutation
#       contig_phasing_result
     
def read_TP(filename):

    f = open(filename, "r")
    words = []
    for line in f:
        words = line.strip().split(' ')

    return words



def read_mutation_record(filename):
    
    f = open(filename, "r")
    mutation = {}
    for line in f:  
        words = line.strip().split()
        mutation[words[0]] = (words[1], words[2])
    #print ("real snp pos:")
    #print (pos)
    return mutation 

def read_snp_mutation(filename):

    f = open(filename, "r")
    mutation = {}
    for line in f:  
        words = line.strip().split()
        mutation[words[0]] = (words[1], words[3])
    return mutation

def read_phasing_result(filename):

    f = open(filename, "r")
    haplotype = {} 
    lineNumber = 0
    for line in f:
        lineNumber += 1
        if lineNumber == 6:  
            words = line.strip().split(',')
        if lineNumber == 7:
            binarySeq = line.strip()
    print (len(binarySeq))
    print (len(words))
    assert len(binarySeq) == len(words)
    for i in range(len(words)):
        haplotype[words[i]] = binarySeq[i] 
    return haplotype
 
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

   
    if len(sys.argv) < 6:
        print ("python " + sys.argv[0] + "TP_in_ref, TP_in_contig, mutation_record, contig_snp_mutation, contig_phasing_result")
        sys.exit()

    TP_in_ref = read_TP(sys.argv[1])
    print (TP_in_ref) 
 
    TP_in_con = read_TP(sys.argv[2])
    print (TP_in_con) 
    # in reference
    mutation_record = read_mutation_record(sys.argv[3])
    print (sorted(mutation_record))

    # in contig
    snp_mutation = read_snp_mutation(sys.argv[4])
    print (sorted(snp_mutation))

    haplotype = read_phasing_result(sys.argv[5])
    #print (start, end)
    print ((sorted(haplotype)))       
    
    #fout = open("haplotype_compare","w")
    hap_both_con = "" 
    for cc in TP_in_con:
        if cc not in haplotype:
            continue
        hap_both_con = hap_both_con + haplotype[cc]

    print (hap_both_con)
    hap_both_ref = ""
    for i in range(len(TP_in_con)):
        pos_ref = TP_in_ref[-(i+1)]
        pos_con = TP_in_con[i]
        #print (pos_con, pos_ref)
        #if pos_ref not in mutation_record or pos_con not in snp_mutation:
        #    break 
        if pos_con not in haplotype:
            continue
        #print (pos_ref, pos_con) 
        #print (mutation_record[pos_ref], snp_mutation[pos_con])
        if (mutation_record[pos_ref][0] == snp_mutation[pos_con][0] or 
            mutation_record[pos_ref][0] == reverseNucleotide(snp_mutation[pos_con][1])):
            hap_both_ref = hap_both_ref + "0"
        elif (mutation_record[pos_ref][0] == snp_mutation[pos_con][1] or
            mutation_record[pos_ref][0] == reverseNucleotide(snp_mutation[pos_con][0])):
            hap_both_ref = hap_both_ref + "1"
        else:
            hap_both_ref = hap_both_ref + "2"

    #print (hap_both_ref)
    print (tools.reverse_Bool(hap_both_ref))
    #print (tools.hamming_Distance(hap_both_ref, hap_both_con))
 
    print (tools.hamming_Distance(hap_both_ref, tools.reverse_Bool(hap_both_con)))
