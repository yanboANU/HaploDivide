#!/usr/bin/python


import os
import sys
import random

snpRate = 0.001 #or snp 1
IndelRate = 0.0001
chrLen = 243000000 #

def read_positions(filename):
    f = open(filename, "r")
    pos = {}
    for line in f:
        words = line.split()
        pos[ int(words[0]) ] = (words[1], words[2])
    f.close()
    print len(pos)
    return pos;


def generate_random(pos, cnt):
    posLen = len(pos)
    if posLen < cnt:
        print "cannot generate too many snps"
        sys.exit()
    posChoose = set()
    while cnt>0:
        p = random.randint(0,posLen-1)
        while pos[p] in posChoose:
            p = random.randint(0,posLen-1)
        posChoose.add(pos[p])
        cnt = cnt-1
    return sorted(posChoose)   
 

if __name__ == "__main__":

    #
    if(len(sys.argv) < 2):
        print "Usage: python generate_snp_pos chrId_snps snp_position_filename "
        sys.exit()

    pos = read_positions(sys.argv[1])
    #cnt = snpRate*chrLen
    cnt = IndelRate*chrLen

    print "generate snps number: ", cnt
    posChoose = generate_random(pos.keys(), cnt)


    fout = open(sys.argv[2], "w")
    for k in posChoose:
        #print k, pos[k][0], pos[k][1]
        #print k, pos[k]
        fout.write("%s %s %s\n" % (k, pos[k][0], pos[k][1]))
    fout.close() 

