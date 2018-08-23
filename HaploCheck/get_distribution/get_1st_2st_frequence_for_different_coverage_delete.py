import sys
import read
import tools  

def read_columns(filename):
    f = open(filename, "r")
    covColumns ={}
    contentColumns = {}
    for line in f:
        words = line.strip().split(" ")
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
            #if words[0] != '*':    
            covColumns[coverage][pos].append((float(words[1])/coverage, words[0]))
    f.close()
    for cov in covColumns:
        for pos in covColumns[cov]:
            covColumns[cov][pos] = sorted(covColumns[cov][pos], reverse=True)
    return covColumns, contentColumns    
    

    

if __name__ == "__main__":
    # block1: 585989, 2702781,
    # block2: 2746291, 12954384,
    # block5: 29553836, 121757928
    #real_deletePosition, real_deleteContent = read.read_delete2(sys.argv[1], 585989, 2702781, -585989) # mutation_record
    #real_deletePosition, real_deleteContent = read.read_delete2(sys.argv[1], 2746291, 12954384, -2746291) # mutation_record


    if len(sys.argv) < 4:
        print ("python " + sys.argv[0] + " mutation_record *column which_block average_coverage")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3]+ " " + sys.argv[4])

    if sys.argv[3] == "1":
        start, end, base = 585989, 2702781, -585989
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553836
   
    ignoreLen = 10000
    real_deletePosition, real_deleteContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base+1) # mutation_record
    covColumns,contentColumns = read_columns(sys.argv[2]) # 
    '''  
    for key in real_deleteContent: 
        if key not in contentColumns: if key not in contentColumns, means coverage equal to 0
            print key
        assert key in contentColumns
    '''    
    
    #for cov in covColumns:
    

    averageCov = int(sys.argv[4])  
    for cov in range(2*averageCov - 1, 2*averageCov + 2):
        delete = {}
        nonDelete= {}
        for key in covColumns[cov]:
            case = covColumns[cov][key]
            if key in real_deletePosition:
                if real_deleteContent[key][1] != contentColumns[key-1] and real_deleteContent[key][1] != contentColumns[key-1].upper():
                    print "maybe something wrong"
                    print real_deleteContent[key], contentColumns[key] 
                assert real_deleteContent[key][1] == contentColumns[key-1] or real_deleteContent[key][1] == contentColumns[key-1].upper() 
                '''
                if cov == 8 and covColumns[cov][key][0] >= 0.9:
                    print "delete", key, covColumns[cov][key], real_deleteContent[key] 
                if len(case) == 0:
                    print "nothing", key, cov
                    assert cov <= 5
                    continue
                '''    
                assert len(case) >= 1
                if len(case) <= 1:
                    print case
                    assert case[0][1] == "*"
                    F1st_2nd = (case[0][0], 0 ) 
                else:    
                    print case
                    assert case[0][1] == "*" or case[1][1] == "*"
                    F1st_2nd = (case[0][0], case[1][0]) 
                if F1st_2nd not in delete:
                    delete[ F1st_2nd ] = 0
                delete[ F1st_2nd ] = delete[ F1st_2nd ] + 1
            else:
                if len(case) <= 1:
                    if len(case) == 0:
                        print "nothing", key, cov
                        assert cov <= 5 
                        continue
                    F1st_2nd = (case[0][0], 0 ) 
                else:   
                    '''
                    if cov == 8 and covColumns[cov][key][1] >= 0.3:
                        print "non",key, covColumns[cov][key]
                    '''    
                    F1st_2nd = (case[0][0], case[1][0] )  
                if F1st_2nd not in nonDelete:
                    nonDelete[ F1st_2nd ] = 0
                nonDelete[ F1st_2nd ] = nonDelete[ F1st_2nd ] + 1

        #print cov, len(delete), len(nonDelete)
        if len(covColumns[cov]) > 100000 and cov>5:
            foutDelete = open(str(cov)+"delete_F1st_2nd_frequence.txt","w")
            for (v1,v2) in delete:
                foutDelete.write('%.2f %.2f %d\n' % (v1, v2, delete[(v1, v2)]))
            foutDelete.close()
            
            foutNon = open(str(cov)+"nondelete_F1st_2nd_frequence.txt","w")
            for (v1,v2) in nonDelete:
                foutNon.write('%.2f %.2f %d\n' % (v1, v2, nonDelete[(v1, v2)]))
            foutNon.close()    
    
