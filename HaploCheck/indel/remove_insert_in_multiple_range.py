import sys
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord


def read_multiple_pos(filename, start, end, base):
    multiplePos = [] 
    f = open(filename, "r")
    
    for line in f:
        words = line.strip().split(" ")
        val = int(words[0])
        if val <start:
            continue
        if val > end:
            break
        multiplePos.append(val+base)
    print "finish reading mutiple pos", len(multiplePos)
    return multiplePos    


def read_snp(filename, multiplePos):

    f = open(filename, "r")
    fout = open("filter_in_multiple_range_insert", "w")
    count1 = 0
    count2 = 0
    i = 0
    MPosLen = len(multiplePos)
    for line in f:  
        words = line.strip().split()
        if len(words) >= 5 or len(words)==3: # or len(words)==3: #for file contig_snp_mutation
            if words[0].startswith("in"):
                continue
            pos = int(words[0])
            count1 += 1
            '''
            if pos not in multiplePos:
                fout.write(line)
                count2 += 1
            '''
            #print i
            if pos == multiplePos[i]:
                i = i+1
                continue
            if i > MPosLen:
                fout.write(line)
                count2 += 1
                continue
            if i<MPosLen and pos > multiplePos[i]:
                while i<MPosLen and pos > multiplePos[i]:
                    i = i+1
                if i<MPosLen and pos == multiplePos[i]:
                    continue
                elif i<MPosLen and pos < multiplePos[i]:
                    fout.write(line)
                    count2 += 1
                    continue
            if i<MPosLen and pos < multiplePos[i]:
                fout.write(line)
                count2 += 1
    f.close()
    fout.close()
    print ("before filter delete in multiple range, insert number:", count1) 
    print ("after filter delete in multiple range, insert number:", count2)


if __name__ == "__main__":


    #multiple pos
    #*snp_mutation
    #*which block
    if sys.argv[3] == "0":
        #start, end, base = 0, sys.maxint, -1
        start, end, base = 0, sys.maxint, 0
    
    if sys.argv[3] == "10":
        start, end, base = 0, 10000000, 0
    if sys.argv[3] == "1":
        start, end, base = 585989, 2702781, -585988
    if sys.argv[3] == "2":
        start, end, base = 2746291, 12954384, -2746290 
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
    
    ignoreLen = 10000
    
    multiplePos = read_multiple_pos(sys.argv[1], start + ignoreLen, end-ignoreLen, base) # start with 1
    # sys.argv[4] is i or d
    #read_snp(sys.argv[2], multiplePos, sys.argv[4])
    
    read_snp(sys.argv[2], multiplePos) # start with 0
    

