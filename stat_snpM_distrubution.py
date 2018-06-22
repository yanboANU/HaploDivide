#!/usr/bin/env python3

import sys
import os 
import string


def read_file_one_by_one(filename):
    f = open(filename, "r")
    one_dis_snp = []
    for line in f:
        words = line.split()
        if len(words)>0 and words[0].isdigit():
            for i in range(len(words)-1):    
                one_dis_snp.append(int(words[i+1]) -int(words[i]))  

    return one_dis_snp

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        print "Usage: python snp_mutation_list"
        sys.exit()
        
    f = open(sys.argv[1], "r")
    dis_snp = []
    for line in f:
        one_dis_snp = read_file_one_by_one(line.strip())
        dis_snp.extend(one_dis_snp)
    
    dis_snp.sort()
    for cc in dis_snp:
        print (cc)
    #print (dis_snp)
    #fout = open()
        

    
