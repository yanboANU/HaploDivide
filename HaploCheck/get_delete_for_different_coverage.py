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
                covColumns[coverage][pos] = 0
            if words[0] == '*':    
                covColumns[coverage][pos] +=  float(words[1])/coverage
    f.close()

    # two delete sigal add together
    # yu think push gap more stable
    '''
    for cov in covColumns:
        for pos in covColumns[cov]:
            if (pos+1 in contentColumns) and (contentColumns[pos] == contentColumns[pos+1]) and (pos+1 in covColumns[cov]):
                covColumns[cov][pos] += covColumns[cov][pos+1] 
                covColumns[cov][pos+1] = 0
    '''

    #print sorted(contentColumns.items())
    #print sorted(covColumns.items())
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

    if len(sys.argv) < 6:
        print ("python " + sys.argv[0] + " mutation_record *column chr1_multuiple_pos which_block averageCov")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3] + " "+ sys.argv[4]+ " "+ sys.argv[5])

    if sys.argv[4] == "1":
        start, end, base = 585989, 2702781, -585989
    if sys.argv[4] == "5":
        start, end, base = 29553836, 121757928 , -29553836
    
    # we ignore first and last 10k
    ignoreLen = 10000
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start + ignoreLen, end-ignoreLen, base+1) # mutation_record, record position of delete is previous position, so base need + 1
    covColumns,contentColumns = read_columns(sys.argv[2]) #

    '''
    for key in real_snpPosition:
        if key not in contentColumns: # if key not in contentColumns, means coverage equal to 0
            print key
        assert key in contentColumns
    '''
    multiplePos = read_multiple_pos(sys.argv[3], start + ignoreLen, end-ignoreLen, base)

    averageCov = int(sys.argv[5])  

    for cov in range(2*averageCov - 1, 2*averageCov + 2):
        snp = {}
        nonSnp= {}
        snpPoss = {}
        nonSnpPoss = {}
        for key in covColumns[cov]:
            if key in real_snpPosition:
                assert key-1 not in multiplePos
                if real_snpContent[key][1] != contentColumns[key-1] and real_snpContent[key][1] != contentColumns[key-1].upper():
                    print "maybe something wrong"
                    print real_snpContent[key], contentColumns[key-1], contentColumns[key] 
                assert real_snpContent[key][1] == contentColumns[key-1] or real_snpContent[key][1] == contentColumns[key-1].upper() 
                rate = covColumns[cov][key] 
                if cov == 2*averageCov and rate <= 0.1:
                    print "snp",key, rate, real_snpContent[key]     
                if rate not in snp:
                    snp[ rate ] = 0
                    snpPoss[ rate ] = []
                snp[ rate ] += 1
                snpPoss[ rate ].append(key)
            elif key-1 not in multiplePos:
                rate = covColumns[cov][key] 
                if rate >= 0.3 and cov == 2*averageCov:
                    print "non", key, rate
                if rate not in nonSnp:
                    nonSnp[ rate ] = 0
                    nonSnpPoss[ rate ] = []
                nonSnp[ rate ] += 1
                nonSnpPoss[ rate ].append(key)

        #print cov, len(snp), len(nonSnp)
        if len(covColumns[cov]) > 100000 and cov>5:
            foutSnp = open(str(cov)+"SNP_rate_frequence.txt","w")
            for v in snp:
                foutSnp.write('%.2f %d\n' % (v, snp[v]))
                #print (v, snpPoss[v])
                for pos in sorted(snpPoss[v]):
                    foutSnp.write("%s " % pos)
                foutSnp.write("\n")    
            foutSnp.close()
            
            foutNon = open(str(cov)+"nonSNP_rate_frequence.txt","w")
            for v in nonSnp:
                foutNon.write('%.2f %d\n' % (v, nonSnp[v]))
                #print (v, nonSnpPoss[v])
                for pos in sorted(nonSnpPoss[v]):
                    foutNon.write("%s " % pos)
                foutNon.write("\n")    
            foutNon.close()    
    
