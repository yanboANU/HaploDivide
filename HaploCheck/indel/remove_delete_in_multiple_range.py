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
    fout = open("filter_in_multiple_range_delete", "w")
    count1 = 0
    count2 = 0
    i = 0
    MPosLen = len(multiplePos)
    for line in f:  
        words = line.strip().split()
        if len(words) >= 5 or len(words)==2 or len(words)==3: #for file contig_snp_mutation
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
    print ("before filter delete in multiple range, delete number:", count1) 
    print ("after filter delete in multiple range, delete number:", count2)


if __name__ == "__main__":
    #multiple pos
    #*snp_mutation
    #*which block
    if len(sys.argv) < 3:
        print ("python " + sys.argv[0] + " mutiple_Pos pre_snp which_block")
        sys.exit()

    print ("python " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2] + " " + sys.argv[3])
    
    if sys.argv[3] == "0":
        #start, end, base = 0, sys.maxint, -1  #-2+1=-1,
        start, end, base = 0, sys.maxint, 0  #real and pre and multiple position is accord

    if sys.argv[3] == "10": # first 10M

        start, end, base = 0, 10000000, 0  #real and pre and multiple position is accord
        #start, end, base = 0, sys.maxint, -1  #-2+1=-1,
    if sys.argv[3] == "1":
        #start, end, base = 585989, 2702781, -585989 # 585989 reference start position 1
        start, end, base = 585989, 2702781, -585988 # 585989 reference start position 1
    if sys.argv[3] == "2":
        start, end, base = 2746291, 12954384, -2746290 
    if sys.argv[3] == "5":
        start, end, base = 29553836, 121757928 , -29553835
    
    ignoreLen = 10000

    #                        block 1 
    # reference 1 2 3 ....   585989 ...      2702781
    # read_snp               0          2702781 - 585989
    # some delete record exactly delete position
    # some delete record delete position -1, so need base+1

    # in real data

    # reference 1 2 3 ....   i ...      2702781
    # read_snp  0 1 2 ....  i-1
    # some delete record exactly delete position, base (snp138)
    # some delete record delete position -1, so need base+1(snp150)

    
    multiplePos = read_multiple_pos(sys.argv[1], start + ignoreLen, end-ignoreLen, base) 
    #multi_poe start with 1 , delete record = real-1( the result of base + 1 )
    #print multiplePos[0:10]
    #sys.exit()
    
    read_snp(sys.argv[2], multiplePos) #start with 0 => 1

    

