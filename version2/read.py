import sys
import column



def read_align(filename):
    # read align_pos 
    f = open(filename, "r")
    ref_pos = []
    contig_pos = []
    ref_label = 1 
    for line in f:
        words = line.split(',')
        if words[0].isdigit():
            if len(ref_pos) == 0:
                for c in words:
                    ref_pos.append(int(c))   
            else: 
                for c in words:
                    contig_pos.append(int(c))
    
    assert len(ref_pos) == len(contig_pos)
    '''
    fout = open("ref_contig", "w")  
    for i in range(len(ref_pos)):
        fout.write("%s %s\n" % (ref_pos[i], contig_pos[i]))
    '''
    return ref_pos, contig_pos

def read_mutation_list(filename):
    f = open(filename, "r")

    pos = []  
    for line in f:  
        words = line.split()
        if len(words) == 1 or len(words) == 3:
            pos.append(int(words[0]))
    
    return pos


def read_columns(filename):
    f = open(filename, "r")
    columns = {}
    for line in f:
        words = line.strip().split() 
        if line.startswith("reference"):
            #print (words[1])
            pos = int(words[2]) 
            assert pos not in columns
            columns[pos] = column.Column(pos)
            columns[pos]._cov = int(words[4])  
        elif (len(words) > 0 ):
            #columns[pos]._print()
            columns[pos]._diff_cov.append(int(words[1]))
    for pos in columns:
        columns[pos]._diff_cov.sort() 
        columns[pos]._diff_cov.reverse() 
    #columns[0]._print()
    #columns[2925]._print() 
    return columns       
#read_columns(sys.argv[1]) 
#read_align(sys.argv[1])
