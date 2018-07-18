import sys
import read
import tools  

def read_columns(filename):
    f = open(filename, "r")
    covColumns ={}

    contentColumns = {}
    for line in f:
        #print line
        words = line.strip().split(" ")
        #print len(words)
        if len(words) <= 1:
            continue
        if words[0].startswith("reference"):
            pos = int(words[2])
            coverage = int(words[5])
            content = words[8] 
            contentColumns[pos] = content
        else:
            
            if coverage not in covColumns:
                covColumns[coverage] = {}
                
            if pos not in covColumns[coverage]:
                covColumns[coverage][pos] = []
            if words[0] != '*':    
                covColumns[coverage][pos].append(float(words[1])/coverage)
            #print pos, columns[pos]
    f.close()
    for cov in covColumns:
        for pos in covColumns[cov]:
            covColumns[cov][pos] = sorted(covColumns[cov][pos], reverse=True)
    return covColumns, contentColumns    
    

    

if __name__ == "__main__":
    # block1: 585989, 2702781,
    # block2: 2746291, 12954384,
    # block5: 29553836, 121757928
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 585989, 2702781, -585989) # mutation_record
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 2746291, 12954384, -2746291) # mutation_record

    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 29553836, 121757928 , -29553836) # mutation_record
    covColumns,contentColumns = read_columns(sys.argv[2]) # 
   
    SNPCount = {}


    for key in real_snpContent:
        assert key in contentColumns

    for cov in covColumns:
        snp = {}
        nonSnp= {}
        for key in covColumns[cov]:
            if key in real_snpPosition:
                if real_snpContent[key][0] != contentColumns[key] and real_snpContent[key][0] != contentColumns[key].upper():
                    print "maybe something wrong"
                    print real_snpContent[key][0], contentColumns[key] 
                assert real_snpContent[key][0] == contentColumns[key] or real_snpContent[key][0] == contentColumns[key].upper() 
                if cov == 8 and covColumns[cov][key][0] >= 0.9:
                    print "snp",key, covColumns[cov][key], real_snpContent[key] 
                
                if len(covColumns[cov][key]) == 0:
                    print "nothing", key, cov
                    assert cov <= 5
                    continue
                if len(covColumns[cov][key]) <= 1:
                    F1st_2nd = (covColumns[cov][key][0], 0 ) 
                else:    
                    F1st_2nd = (covColumns[cov][key][0], covColumns[cov][key][1] ) 
                if F1st_2nd not in snp:
                    snp[ F1st_2nd ] = 0
                snp[ F1st_2nd ] = snp[ F1st_2nd ] + 1
            else:
                if len(covColumns[cov][key]) <= 1:
                    if len(covColumns[cov][key]) == 0:
                        print "nothing", key, cov
                        assert cov <= 5 
                        continue
                    F1st_2nd = (covColumns[cov][key][0], 0 ) 
                else:   

                    if cov == 8 and covColumns[cov][key][1] >= 0.3:
                        print "non",key, covColumns[cov][key]
                    F1st_2nd = (covColumns[cov][key][0], covColumns[cov][key][1] )  
                if F1st_2nd not in nonSnp:
                    nonSnp[ F1st_2nd ] = 0
                nonSnp[ F1st_2nd ] = nonSnp[ F1st_2nd ] + 1

        #print cov, len(snp), len(nonSnp)
        if len(covColumns[cov]) > 100000 and cov>5:
            foutSnp = open(str(cov)+"SNP_F1st_2nd_frequence.txt","w")
            for (v1,v2) in snp:
                foutSnp.write('%.2f %.2f %d\n' % (v1, v2, snp[(v1, v2)]))
            foutSnp.close()
            
            foutNon = open(str(cov)+"nonSNP_F1st_2nd_frequence.txt","w")
            for (v1,v2) in nonSnp:
                foutNon.write('%.2f %.2f %d\n' % (v1, v2, nonSnp[(v1, v2)]))
            foutNon.close()    
    
