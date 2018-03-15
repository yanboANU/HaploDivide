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
    #SNPNumber = 2
    
    chromatid = genome
    #print "fist genome: ", genome
    mutations = {}
    for i in range(SNPNumber):       #[0,SNPNumber) 
        #print "mutation number: ", i
        if i != 0:
            genome = chromatid
        pos = random.randint(0, len(genome)-1)  #index [0, len(genome)-1] 
        while pos in mutations:
            pos = random.randint(0, len(genome)-1)  #index [0, len(genome)-1] 
        snp_c = SNP(genome[pos]) 
        chromatid = genome[:pos] + snp_c +  genome[pos+1:]
        mutations[pos] = (genome[pos], snp_c)
        #print "mutation position: ", pos
        #print genome
        #print chromatid
    #print "final genome: ", chromatid
    genome = chromatid
    ##############
    #Indel
    ##############
   
    IndelNumber = int(len(genome)*mutationRate*(1-SNPRate))
    #IndelNumber = 2

    indels = {}
    for i in range(IndelNumber):
        if i != 0:
            genome = chromatid
        pos = random.randint(0, len(genome)-1)
        while pos in indels:
            pos = random.randint(0, len(genome)-1)
        indel_c = Indel(genome[pos:pos+longestIndel]) 
        chromatid = genome[:pos] + indel_c + genome[pos+longestIndel:]
        indels[pos] = (genome[pos:pos+longestIndel],indel_c)  
       
        #print "mutation position: ", pos
        #print genome
        #print chromatid

    print (mutations)
    sorted_m = (sorted(mutations.items()))
    sorted_i = (sorted(indels.items()))
 
    fout1 = open("mutation_record", "w")
    fout2 = open("indel_record", "w")
    for c in sorted_m:
        fout1.write("%s %s %s\n" % (c[0], c[1][0], c[1][1]) )
       
    for c in sorted_i:
        fout2.write("%s %s %s\n" % (c[0], c[1][0], c[1][1]) ) 
    return chromatid 

if __name__ == "__main__":
    
    if(len(sys.argv) < 4):
        print "Usage: python give_mutation.py genome_file mutation_genome_file parameter"
        sys.exit()
#two kind mutation: SNP and small indel
#parameter include: (1),mutation rate(1%-5%) (2),SNP rate (3), longest indel

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
                chromatid = mutate(genome, mutationRate, SNPRate, longestIndel)
                nameA = name.split("_")    
	        fOut.write(str(nameA[0]) + "_dual_" + str(nameA[1])  + "|chromatid\n")
	        fOut.write(chromatid + "\n")
            name = line
            genome = ""
        else:
            genome = genome + line

    print "genome name:", name
    print "genome length:", len(genome)
    chromatid = mutate(genome, mutationRate, SNPRate, longestIndel)    
    
    nameA = name.split("_")    
    fOut.write(str(nameA[0]) + "_dual_" + str(nameA[1])  + "|chromatid\n")
    fOut.write(chromatid + "\n")    

    fIn.close()
    fOut.close()
 
        

 
 




 
