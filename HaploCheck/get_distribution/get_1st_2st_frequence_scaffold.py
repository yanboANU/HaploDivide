import sys
import read
import tools  

def read_columns(filename):
    columns= {}
    f = open(filename, "r")
    for line in f:
        #print line
        words = line.strip().split(" ")
        #print len(words)
        if len(words) <= 1:
            continue
        if words[0].startswith("reference"):
            pos = int(words[2])
            coverage = int(words[5]) 
        else:
            if coverage <= 10:
                continue
            if pos not in columns:
                columns[pos] = []
            columns[pos].append(float(words[1])/coverage)
            #print pos, columns[pos]
    for key in columns:
        columns[key] = sorted(columns[key], reverse=True)
    return columns    
    

    

if __name__ == "__main__":
    
    
    real_snpPosition, real_snpContent = read.read_snp(sys.argv[1]) # mutation_record
    columns = read_columns(sys.argv[2]) # *columns
    alignScaffoldRef = read.read_blasr_m5(sys.argv[3]) #ref_scaffols.blasr.m5
    alignScaffoldRef._left._generate_seq_pos() # real
    alignScaffoldRef._right._generate_seq_pos()# pre    
    leftToRight, rightToLeft = alignScaffoldRef._generate_pos_to_pos()

    # real_snpPosition => position correspond to reference
    right_snpPosition, not_find = tools.convert(leftToRight, real_snpPosition) 
    print "not find number", len(not_find)
    
    SNPCount = {}

    print "SNPs strange case" 
    for key in right_snpPosition:
        #if key not in columns:
            #continue 
        if key not in columns:
            continue
        if len(columns[key]) <= 1:
            #print columns[key][0], 0
            print key, columns[key]
            if (columns[key][0], 0) not in SNPCount:
                SNPCount[ (columns[key][0], 0)] = 0
            SNPCount[ (columns[key][0], 0)] += 1
        else:
            if columns[key][0] > 0.78:
                print key, columns[key]
            #print columns[key][0], columns[key][1]
            if (columns[key][0], columns[key][1]) not in SNPCount:
                SNPCount[ (columns[key][0], columns[key][1]) ] = 0
            SNPCount[ (columns[key][0], columns[key][1]) ] += 1


    print "non SNPs strange case"
    nonSNPCount = {}
    for key in columns:
        if key not in right_snpPosition:
            if len(columns[key]) <= 1:
                if (columns[key][0], 0) not in nonSNPCount:
                    nonSNPCount[ (columns[key][0], 0) ] = 0
                nonSNPCount[ (columns[key][0], 0) ] += 1
            else:
                if columns[key][0] < 0.7 and columns[key][1] > 0.25:
                    print key, columns[key]
                if (columns[key][0], columns[key][1]) not in nonSNPCount:
                    nonSNPCount[ (columns[key][0], columns[key][1]) ] = 0
                nonSNPCount[ (columns[key][0], columns[key][1]) ] += 1

     
    ''' 
    fout = open("SNP_1st_2nd_frequence.txt", "w")
    print "real SNP Position"
    for (v1,v2) in sorted(SNPCount.items()):
        #print ('%.2f %.2f %d ' % (v1[0], v1[1], v2))
        fout.write('%.2f %.2f %d\n' % (v1[0], v1[1], v2))
    fout.close()


    fout2 = open("nonSNP_1st_2nd_frequence.txt", "w")
    print "non SNP Position"
    for (v1,v2) in sorted(nonSNPCount.items()):
        #print v1[0], v1[1], v2
        #print ('%.2f %.2f %d ' % (v1[0], v1[1], v2))
        fout2.write('%.2f %.2f %d\n' % (v1[0], v1[1], v2))
    fout2.close()    
    '''
