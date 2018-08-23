from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import read
import tools  
import time

def read_Ins(filename):
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
                covColumns[coverage][pos] = [0,0,0,0]
            four = ['A', 'T', 'C', 'G']
            for ele in four:
                if words[0].count(ele) > 0:
                    covColumns[coverage][pos][four.index(ele)] = covColumns[coverage][pos][four.index(ele)] + int(words[1])
    f.close()
    for cov in covColumns:
        for pos in covColumns[cov]:
            #print cov, pos, covColumns[cov][pos]
            covColumns[cov][pos] = max(covColumns[cov][pos])/float(cov)
    return covColumns, contentColumns    
    
def get_cov_for_snp(filename, real_snpPosition):
    f = open(filename, "r")
    covColumns ={}
    for line in f:
        words = line.strip().split(" ")
        if len(words) <= 1:
            continue
        if words[0].startswith("ref"):
            pos = int(words[2])
            coverage = int(words[5])
            if pos in real_snpPosition:
                if coverage not in covColumns:
                    covColumns[coverage] = []
                covColumns[coverage].append(pos)    
    f.close()
    return covColumns


# based on reference, we can get multiple range/pos
# record index in reference AAAA/TTTT(len >= 3)  GCGCGCGC(len >= 8)

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
    
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 2746291, 12954384, -2746291) # mutation_record
    
    if len(sys.argv) < 5:
        print ("python " + sys.argv[0] + " mutation_record *Ins chr1_multiple_pos which_block *columns average_coverage")
        sys.exit() 
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4] + " " +sys.argv[5] + " " + sys.argv[6])
    
    if sys.argv[4] == "1":
        start, end, base = 585989, 2702781, -585989
    if sys.argv[4] == "2":
        start, end, base = 2746291, 12954384, -2746291
    if sys.argv[4] == "5":
        start, end, base = 29553836, 121757928 , -29553836

    time1 = time.clock() 

    multiplePos = read_multiple_pos(sys.argv[3], start, end, base) 

    time2 = time.clock()
    print ( "read multiple pos running %s Seconds" % (time2 - time1) )
    
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start, end, base) # mutation_record
    
    time3 = time.clock()
    print ( "read snp position %s Seconds" % (time3 - time2) )
    
    covColumns, contentColumns = read_Ins(sys.argv[2]) #

    time4 = time.clock()
    print ( "read insertion %s Seconds" % (time4 - time3) )
    # get coverage for all real snp
    # *column
    snpCov = get_cov_for_snp(sys.argv[5], real_snpPosition)
    print snpCov
 
    time5 = time.clock()
    print ( "read columns %s Seconds" % (time5 - time4) )
    averageCov = int(sys.argv[6])
    '''
    for key in real_snpContent:
        if key not in contentColumns:
            print "noReport", key, key-base, "have insert, but report nothing"
        #assert key in contentColumns
    ''' 
    for cov in range(2*averageCov - 1, 2*averageCov + 2):
        snp = {}
        nonSnp= {}
        snp[0] = 0
        snpPoss = {} #
        nonSnpPoss = {} #
        snpPoss[0] = []
        #print cov
        if cov not in snpCov:
            continue
        for key in snpCov[cov]:
            if key in covColumns[cov]:
                if key in multiplePos:
                    print key, key+start
                assert key not in multiplePos
                if real_snpContent[key][0] != contentColumns[key] and real_snpContent[key][0] != contentColumns[key].upper():
                    print "maybe something wrong"
                    print real_snpContent[key][0], contentColumns[key] 
                assert real_snpContent[key][0] == contentColumns[key] or real_snpContent[key][0] == contentColumns[key].upper() 
                
                rate = covColumns[cov][key]
                if cov == 2*averageCov and rate <= 0.1:
                    print "snp", key, rate, real_snpContent[key]
                if rate not in snp:
                    snp[ rate ] = 0
                    snpPoss[ rate ] =[]
                snp[ rate ] += 1
                snpPoss[rate].append(key)
            else:
                snp[0] += 1
                snpPoss[0].append(key)
                print "snp", key, 0, real_snpContent[key]
        for key in covColumns[cov]:
            if key not in snpCov[cov] and key not in multiplePos:
                #if key happen in same mutiple position
                #then don't consider
                #in the range or near range
                #we can find those postion based on ref
                rate = covColumns[cov][key]
                if cov == 2*averageCov and rate >= 0.5:
                    print "non", key, rate
                if rate not in nonSnp:
                    nonSnp[ rate ] = 0
                    nonSnpPoss[rate] = []
                nonSnp[ rate ] += + 1
                nonSnpPoss[rate].append(key)

        #print cov, len(snp), len(nonSnp)
        if len(covColumns[cov]) > 10000:
            foutSnp = open(str(cov)+"Ins_F1st_2nd_frequence.txt","w")
            for v1 in snp:
                foutSnp.write('%.2f %d\n' % (v1, snp[v1]))
                for pos in sorted(snpPoss[v1]):
                    foutSnp.write("%s " % pos)
                foutSnp.write("\n")    
            foutSnp.close()
            
            foutNon = open(str(cov)+"nonIns_F1st_2nd_frequence.txt","w")
            for v1 in nonSnp:
                foutNon.write('%.2f %d\n' % (v1, nonSnp[v1]))
                for pos in sorted(nonSnpPoss[v1]):
                    foutNon.write("%s " % pos)
                foutNon.write("\n")    
            foutNon.close()    
