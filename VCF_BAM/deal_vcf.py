import sys


def read_vcf(filename):
    f = open(filename, "r")
    insert1s = {}
    delete1s = {}
    snps = {}
    for line in f:
        if line.startswith('#'):
            continue
        words = line.split()
        chrID = words[0]
        if chrID not in snps:
            insert1s[chrID] = []
            delete1s[chrID] = []
            snps[chrID] = []
        s1 = words[3].strip()
        s2 = words[4].strip()
        homo = words[9].split(':')[0]
        #print homo
        #sys.exit()
        s1Len = len(s1)
        s2Len = len(s2)
        if s1 == '.' or s2 == '.':
            continue
        
        #print s1Len, s2Len
        #if chrID != "1":
            #break
           
        
        if s1Len == 1 and s2Len == 1:
            snps[chrID].append((words[1], s1, s2, homo))
            continue
                
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
                    snps[chrID].append((words[1], s1, s2, homo))
                    continue

        if s1Len == 1 and s2Len >= 2 and s2[0] == s1:
            if s2Len == 3 and s2.count(',')>=1:
                print s2
                continue
            insert1s[chrID].append((words[1], s1, s2, homo))
            continue        
           
        
        if s1Len >= 2 and s2Len == 1 and s1[0] == s2:
            delete1s[chrID].append((words[1], s1, s2, homo))
            continue
         
    return insert1s, delete1s, snps    
 



if __name__ == "__main__":
    if(len(sys.argv) <= 1):
        print "Usage: python *vcf"
        sys.exit()

    insert1s, delete1s, snps = read_vcf(sys.argv[1])    
    
    for key in snps:
        print key, len(insert1s[key]) 
        if len(insert1s[key]) > 0:
            foutI = open(key+"_inserts_all", "w")
            for c in insert1s[key]:
                foutI.write("%s %s %s %s\n" % (c[0],c[1],c[2],c[3]))
            foutI.close()
       
        if len(delete1s[key]) > 0: 
            foutD = open(key+"_deletes_all", "w")
            for c in delete1s[key]:
                #foutD.write("%s %s %s\n" % (c[0],c[1],c[2]))
                foutD.write("%s %s %s %s\n" % (c[0],c[1],c[2],c[3]))
            foutD.close()
       
        if len(snps[key]) > 0: 
            foutS = open(key+"_snps_all", "w")
            for c in snps[key]:
                foutS.write("%s %s %s %s\n" % (c[0],c[1],c[2],c[3]))
            foutS.close()
           
