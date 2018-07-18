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
            if words[0] == '*':    
                covColumns[coverage][pos] = float(words[1])/coverage
    f.close()

    print sorted(contentColumns.items())
    print sorted(covColumns.items())
    return covColumns, contentColumns    

    
def read_multiple_pos(filename, start, end, base):
    multiplePos = set() 
    f = open(filename, "r")
    for line in f:
        words = line.strip().split(" ")
        val = int(words[0])  
        if val > end or val <start:
            continue
        multiplePos.add(val+base)

    return multiplePos    


if __name__ == "__main__":
    # block1: 585989, 2702781,
    # block2: 2746291, 12954384,
    # block5: 29553836, 121757928
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 585989, 2702781, -585989) # mutation_record
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 2746291, 12954384, -2746291) # mutation_record

    if sys.argv[4] == "1":
        start, end, base = 585989, 2702781, -585989
    if sys.argv[4] == "5":
        start, end, base = 29553836, 121757928 , -29553836
    
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start, end, base) # mutation_record
    covColumns,contentColumns = read_columns(sys.argv[2]) # 
    sys.exit() 
    SNPCount = {}

    ''' 
    for key in real_snpContent:
        if key not in contentColumns:
            print "noReport", key, "have delete, but report nothing"
        #assert key in contentColumns
    '''
    multiplePos = read_multiple_pos(sys.argv[3], start, end, base) 
    for cov in range(7,10):
        snp = {}
        nonSnp= {}
        for key in covColumns[cov]:
            if key in real_snpPosition:

                assert key not in multiplePos
                if real_snpContent[key][1] != contentColumns[key] and real_snpContent[key][1] != contentColumns[key].upper():
                    print "maybe something wrong"
                    print real_snpContent[key], contentColumns[key] 
                assert real_snpContent[key][1] == contentColumns[key] or real_snpContent[key][1] == contentColumns[key].upper() 
                if cov == 8 and covColumns[cov][key] >= 0.9:
                    print "snp",key, covColumns[cov][key], real_snpContent[key] 
                
                if len(covColumns[cov][key]) == 0:
                    print "nothing", key, cov
                    #assert cov <= 5
                    continue
                if len(covColumns[cov][key]) <= 1 and covColumns[cov][key][0][1] == '*':
                    F1st_2nd = (covColumns[cov][key][0][0], 0 ) 
                elif  covColumns[cov][key][0][1] == '*' or covColumns[cov][key][1][1] == '*':
                    F1st_2nd = ( covColumns[cov][key][0][0], covColumns[cov][key][1][0] ) 
                if F1st_2nd not in snp:
                    snp[ F1st_2nd ] = 0
                snp[ F1st_2nd ] = snp[ F1st_2nd ] + 1
            elif key not in multiplePos:
                print cov, key, covColumns[cov][key]
                if len(covColumns[cov][key]) <= 1 and covColumns[cov][key][0][1] == '*':
                    if len(covColumns[cov][key]) == 0:
                        print "nothing", key, cov
                        assert cov <= 5 
                        continue
                    F1st_2nd = (covColumns[cov][key][0][0], 0 )    
                elif len(covColumns[cov][key]) >=2 and (covColumns[cov][key][0][1] == '*' or covColumns[cov][key][1][1] == '*'):
                    if cov == 8 and covColumns[cov][key][1][0] >= 0.3:
                        print "non",key, covColumns[cov][key]
                    F1st_2nd = (covColumns[cov][key][0][0], covColumns[cov][key][1][0] )  
                else:
                    continue
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
    
