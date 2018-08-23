#!/usr/bin/python

import os
import sys
import string
import random
from itertools import ifilter,imap

def read_para(f):

    mutationRate = 0.01
    SNPRate = 0.5
    longestIndel = 10
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        words = line.split(':')
        if words[0] == "mutationRate":
            mutationRate = float(words[1])   
        if words[0] == "SNPRate":
            SNPRate = float(words[1])
        if words[0] == "longestIndel":
            longestIndel = int(words[1])
    f.close() 
    print mutationRate, SNPRate, longestIndel
    return mutationRate, SNPRate, longestIndel


def SNP(Nucleotide):
    alphabet = "ATCG"
    pos = random.randint(1,3) # random produce [1,3], every Nucleotide have equal chance to become other three Nucleotide
    #print "before SNP:", Nucleotide
    #print pos, (pos + alphabet.index(Nucleotide)) % 4
    Nucleotide = alphabet[(pos + alphabet.index(Nucleotide)) % 4]
    #print "after SNP:", Nucleotide
    return Nucleotide

def random_string(l):
    ans = ""
    alphabet = "ATCG"
    for i in range(l):
        pos = random.randint(0,3) # equal probability of creat A,T,C,G
        ans = ans + alphabet[pos] 
    return ans

def Indel_point_length(subS):
    #ll = [100, 100, 100, 200, 200, 200, 400, 800, 1600, 3200]
    ll = [25600]
    InLabel = random.randint(0,1)
    length = len(ll)
    if InLabel == 1:
       InLen = random.randint(0,length-1)
       return random_string(ll[InLen]) + subS
    else:
       DelLen = random.randint(0,length-1)
       #print "delete ", DelLen 
       return subS[ ll[DelLen]:]   


#Indel length random(1, longestIndel)
def Indel(subS):
    length = len(subS) 
    InLabel = random.randint(0,1)
    ###############################
    #equal chance insert or delete 
    #equal chance of different length
    ###############################
    if InLabel == 1:
       InLen = random.randint(1,length)
       #print "insert ", InLen
       return random_string(InLen) + subS
    else:
       DelLen = random.randint(1,length)
       #print "delete ", DelLen 
       return subS[DelLen:]   

def mutate(genome, mutationRate, SNPRate, longestIndel):

    ##############
    #SNP
    ##############
    SNPNumber = int(len(genome)*mutationRate*SNPRate)
    print "SNP number: ", SNPNumber
    
    mutations = {}
    for i in range(SNPNumber):       #[0,SNPNumber) 
        pos = random.randint(0, len(genome)-1)  #index [0, len(genome)-1] 
        while pos in mutations:
            pos = random.randint(0, len(genome)-1)  #index [0, len(genome)-1] 
        snp_c = SNP(genome[pos]) 
        mutations[pos] = (genome[pos], snp_c)
    #print (mutations)
    sorted_m = (sorted(mutations.items()))
    
    start = 0
    chromatid = ""
    for (key, v) in sorted_m:
        chromatid = chromatid + genome[start:key]
        chromatid = chromatid + v[1]
        start = key+1

    #print start
    #print sorted_m[-1]
    chromatid = chromatid + genome[start:len(genome)] 
    
    print "genome length:", len(genome)
    print "mutation ref length:", len(chromatid)

    assert len(genome) == len(chromatid)
    fout1 = open("mutation_record", "w")
    for c in sorted_m:
        fout1.write("%s %s %s\n" % (c[0], c[1][0], c[1][1]) )
    return chromatid 

if __name__ == "__main__":
    
    if(len(sys.argv) < 4):
        print "Usage: python give_mutation.py genome_file mutation_genome_file parameter"
        sys.exit()
#two kind mutation: SNP and small indel
#parameter include: (1),mutation rate(1%-5%) (2),SNP rate (3), longest indel

    print "Usage: python give_mutation.py", sys.argv[1], sys.argv[2], sys.argv[3]
    fpara = open(sys.argv[3], 'r')
    inputFileName = sys.argv[1]
    outputFileName = sys.argv[2]
    fIn = open(inputFileName, 'r')
    fOut = open(outputFileName, 'w')

    genome = ""  
    mutationRate, SNPRate, longestIndel = read_para(fpara)
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip,fIn)):        
        if(line.startswith('>')):
            if len(genome) != 0:
		print "genome name:", name
		print "genome length:", len(genome)
                # origal genome      
	        fOut.write(name +"\n")
	        fOut.write(genome + "\n")
                chromatid = mutate(genome, mutationRate, SNPRate, longestIndel)
                print "mutate genome length:", len(chromatid)
                fOut.write(name + "_dual\n")
	        fOut.write(chromatid + "\n")
            name = line
            genome = ""
        else:
            genome = genome + line

    print "genome name:", name
    print "genome length:", len(genome)

    fOut.write(name +"\n")
    fOut.write(genome + "\n")
    chromatid = mutate(genome, mutationRate, SNPRate, longestIndel)    
    
    print "mutate genome length:", len(chromatid)
    fOut.write(name + "_dual\n")
    fOut.write(chromatid + "\n")    

    fIn.close()
    fOut.close()
 
        

 
 




 
