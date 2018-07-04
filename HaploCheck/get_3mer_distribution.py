import sys
import read
import tools  

def read_1st_2nd_3mer(filename):
    columns= {}
    f = open(filename, "r")
    for line in f:
        #print line
        words = line.strip().split(" ")
        columns[int(words[1])] = (words[2], words[3])    
    return columns    
    

    

if __name__ == "__main__":
    
    
    real_snpPosition, real_snpContent = read.read_snp(sys.argv[1]) # mutation_record
    columns = read_1st_2nd_3mer(sys.argv[2]) # 
   
    ALLSNPCount = {}
    IncludeNONSNPCount = {}
    keys = sorted(columns)

    #for key in sorted(columns):
    content = sorted(columns.items())
    lenf = len(keys) - 2 #window +1
    for i in range(lenf):
        if (keys[i] not in real_snpPosition) or (keys[i+1] not in real_snpPosition) or (keys[i+2] not in real_snpPosition):
            print keys[i], columns[keys[i]]
            if (columns[ keys[i] ][0], columns[ keys[i] ][1]) not in IncludeNONSNPCount:
                IncludeNONSNPCount[ (columns[ keys[i] ][0], columns[ keys[i] ][1]) ] = 0
            IncludeNONSNPCount[ (columns[ keys[i] ][0], columns[ keys[i] ][1]) ] += 1
        else:    
            if (columns[ keys[i] ][0], columns[ keys[i] ][1]) not in ALLSNPCount:
                ALLSNPCount[ (columns[ keys[i] ][0], columns[ keys[i] ][1]) ] = 0
            ALLSNPCount[ (columns[ keys[i] ][0], columns[ keys[i] ][1]) ] += 1
    
    fout = open("ALLSNP_1st_2nd_3mer_frequence.txt","w")
    #print "real SNP Position"
    for (v1,v2) in sorted(ALLSNPCount.items()):
        #print ('%.2f %.2f %d ' % (v1[0], v1[1], v2))
        fout.write('%s %s %d\n' % (v1[0], v1[1], v2))
    fout.close()
    
    #print "non SNP Position"
    fout2 = open("IncludeNONSNP_1st_2nd_3mer_frequence.txt","w")
    for (v1,v2) in sorted(IncludeNONSNPCount.items()):
        #print v1[0], v1[1], v2
        #print ('%.2f %.2f %d ' % (v1[0], v1[1], v2))
        fout2.write('%s %s %d\n' % (v1[0], v1[1], v2))
    fout2.close()    
    
