import sys
import read
import tools  
import contig

def read_delete(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    starNumber = {}
    for line in f:  
        words = line.strip().split()
        if len(words) >= 5: #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            snpContent[pos] = (words[1], words[3])
            if words[1] == "*":
                starNumber[pos] = int(words[2])    
            if words[3] == "*":
                starNumber[pos] = int(words[4])
    return snpPosition, snpContent, starNumber
   


def write_delete_rate_distribution(FP, startNumber, cov, filename, averageCov):
    
    FPRateCount = {}
    fout = open(filename, "w")
    for c in FP:
        if c not in cov:
            print c, "not in cov"
            continue
        if c not in starNumber:
            print c, "not in starNumber"
            continue
        if cov[c] != averageCov:
            continue
        rate = round(starNumber[c]/float(cov[c]),2)
        if rate not in FPRateCount:
            FPRateCount[rate] = 0 
        FPRateCount[rate] += 1
    RC = sorted(FPRateCount.items())
    for (r,c) in RC:
        fout.write("%.2f %d\n" %(r, c))
    fout.close()    



if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print ("python " + sys.argv[0] + " delete_record *delete* which_block contig *_cov average_coverage")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3]+ " " + sys.argv[4] + " " + sys.argv[5] + " " +sys.argv[6])

    
    if sys.argv[3] == "10":
        start, end, base = 0, 10000000, 0 # first 10m
    if sys.argv[3] == "0":
        start, end, base = 0, sys.maxint, 0 # real and pre accord to each other now
    if sys.argv[3] == "1":
        start, end, base = 585989, 2702781, -585988
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
   
    ignoreLen = 10000
    real_deletePosition, real_deleteContent = read.read_snp2(sys.argv[1], start+ignoreLen, end-ignoreLen, base) # mutation_record
    deletePos, deleteContent, starNumber = read_delete(sys.argv[2]) # 

    contigs = contig.read_Contig(sys.argv[4])
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()
    FP, TN, TP = tools.print_TP(real_deletePosition, deletePos, contig)

    cov = read.read_cov(sys.argv[5])
    averageCov = int(sys.argv[6])
    print "write TP distribution"
    write_delete_rate_distribution(TP, starNumber, cov, "TP_delete_rate_distribution_" +str(averageCov), averageCov)

    print "write FP distribution"
    write_delete_rate_distribution(FP, starNumber, cov, "FP_delete_rate_distribution_"+str(averageCov), averageCov)

    '''
    averageCov = int(sys.argv[6])  
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
   ''' 
