from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import read
import tools  

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
            #A T C G
            four = ['A', 'T', 'C', 'G']
            for ele in four:
                if words[0].count(ele) > 0:
                    covColumns[coverage][pos][four.index(ele)] = covColumns[coverage][pos][four.index(ele)] + int(words[1])


            #print pos, columns[pos]
    f.close()
    for cov in covColumns:
        for pos in covColumns[cov]:
            #print cov, pos, covColumns[cov][pos]
            covColumns[cov][pos] = max(covColumns[cov][pos])/float(cov)
    return covColumns, contentColumns    
    

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
    
    #block1: 585989, 2702781, -585989
    #block2: 2746291, 12954384, -2746291
    #block5: 29553836, 121757928 , -29553836

    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 585989, 2702781, -585989) # mutation_record
    #real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], 2746291, 12954384, -2746291) # mutation_record
    
    if len(sys.argv) < 5:
        print ("python " + sys.argv[0] + " mutation_record *Ins chr1_multiple_pos which_block")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4])
    
    if sys.argv[4] == "1":
        start, end, base = 585989, 2702781, -585989
    if sys.argv[4] == "5":
        start, end, base = 29553836, 121757928 , -29553836



    multiplePos = read_multiple_pos(sys.argv[3], start, end, base) 
    
    real_snpPosition, real_snpContent = read.read_snp2(sys.argv[1], start, end, base) # mutation_record
    covColumns, contentColumns = read_Ins(sys.argv[2]) #
   
    #multiplePos = read_multiple_range(sys.argv[3], 585989, 2702781, -585989) 
   
    SNPCount = {}

    #get_multiple_pos(multipleRange)


    for key in real_snpContent:
        if key not in contentColumns:
            print "noReport", key, key-base, "have insert, but report nothing"
        #assert key in contentColumns


    for cov in covColumns:
        snp = {}
        nonSnp= {}
        for key in covColumns[cov]:
            if key in real_snpPosition:
                if key in multiplePos:
                    print key, key+start
                assert key not in multiplePos
                if real_snpContent[key][0] != contentColumns[key] and real_snpContent[key][0] != contentColumns[key].upper():
                    print "maybe something wrong"
                    print real_snpContent[key][0], contentColumns[key] 
                assert real_snpContent[key][0] == contentColumns[key] or real_snpContent[key][0] == contentColumns[key].upper() 
                
                if cov == 10:
                    print "snp", key, covColumns[cov][key], real_snpContent[key]
                 
                F1st_2nd = covColumns[cov][key]
                if F1st_2nd not in snp:
                    snp[ F1st_2nd ] = 0
                snp[ F1st_2nd ] = snp[ F1st_2nd ] + 1
            elif key not in multiplePos:
                #if key happen in same mutiple position
                #then don't consider
                #in the range or near range
                #we can find those postion based on ref
                F1st_2nd = covColumns[cov][key]
                if cov==10 and F1st_2nd >= 0.5:
                    print "non", key, F1st_2nd
                if F1st_2nd not in nonSnp:
                    nonSnp[ F1st_2nd ] = 0
                nonSnp[ F1st_2nd ] = nonSnp[ F1st_2nd ] + 1

        #print cov, len(snp), len(nonSnp)
        if len(covColumns[cov]) > 10000:
            foutSnp = open(str(cov)+"Ins_F1st_2nd_frequence.txt","w")
            for v1 in snp:
                foutSnp.write('%.2f %d\n' % (v1, snp[v1]))
            foutSnp.close()
            
            foutNon = open(str(cov)+"nonIns_F1st_2nd_frequence.txt","w")
            for v1 in nonSnp:
                foutNon.write('%.2f %d\n' % (v1, nonSnp[v1]))
            foutNon.close()    
