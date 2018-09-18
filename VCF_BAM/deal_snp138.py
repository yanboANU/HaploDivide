import sys


def read_txt(filename):
    f = open(filename, "r")
    insert1s = {}
    delete1s = {}
    snps = {}
    for line in f:
        words = line.split()
        chrID = words[1]
        if chrID not in snps:
            insert1s[chrID] = []
            delete1s[chrID] = []
            snps[chrID] = []
        s = words[9].split('/')
        s1 = s[0]
        if len(s) < 2:
            print line
            continue
        s2 = s[1]
        #print s, s[0], s[1]
        #sys.exit()
        s1Len = len(s1)
        s2Len = len(s2)
        
        if s1Len == 1 and s2Len == 1 and s1 != '-' and s2 != '-':
            snps[chrID].append((words[3], s1, s2))
            continue
        
        '''
        if s1Len == 1 and s2Len >= 3:
            sub = s2.split(',')
            subLen = len(sub)
            if subLen >= 2:
                #print s2
                i=0
                flag = True
                for i in range(subLen):
                    if len(sub[i]) != 1:
                        flag = False
                if flag == True:
                    #print s2
                    snps[chrID].append((words[1], s1, s2))
            #sys.exit()
            #continue
        '''
        if s1 == '-' and s2Len == 1:
            #print words[10], words[11]
            #sys.exit()
            if words[11] == "deletion":
                delete1s[chrID].append((words[3], s1, s2))
                continue

            elif words[11] == "insertion":
                insert1s[chrID].append((words[3], s1, s2))
                continue
            
    return insert1s, delete1s, snps    
 



if __name__ == "__main__":
    if(len(sys.argv) <= 1):
        print "Usage: python fasta_to_fastq.py base_dir [name.fasta...]"
        sys.exit()

    insert1s, delete1s, snps = read_txt(sys.argv[1])    
    for key in snps:
         
        foutI = open(key+"_insert1s", "w")
        for c in insert1s[key]:
            foutI.write("%s %s %s\n" % (c[0],c[1],c[2]))
        foutI.close()

        foutD = open(key+"_delete1s", "w")
        for c in delete1s[key]:
            foutD.write("%s %s %s\n" % (c[0],c[1],c[2]))
        foutD.close()
        
        foutS = open(key+"_snps", "w")
        for c in snps[key]:
            foutS.write("%s %s %s\n" % (c[0],c[1],c[2]))
        foutS.close()
