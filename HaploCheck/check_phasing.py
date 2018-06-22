#!/usr/bin/env python3

import sys


def read_result2(filename):
    f = open(filename, "r")
    phase0 = []
    phase1 = []
    count = 1  
    for line in f:
        if line.startswith('S'):
            if count%2 == 1:
                phase0 = line.strip().split(',')
            else:
                phase1 = line.strip().split(',')
            count += 1
        if len(phase0) >0 and len(phase1) >0:
            print ("phase0 size: ", len(phase0))
            print ("phase1 size: ", len(phase1))
            print ("check intersection: ", set(phase0).intersection(set(phase1)))
            small =[]
            large = []  
            for i in phase0:
                if i.startswith('S2'):
                    large.append(i)
                else:
                    small.append(i)
            print ("in phase0, small size: %s large size: %s" % (len(small), len(large)) )
            print ("small", small)
            print ("large", large)   

            small =[]
            large = []  
            for i in phase1:
                if i.startswith('S2'):
                    large.append(i)
                else:
                    small.append(i)
            print ("in phase1, small size: %s large size: %s" % (len(small), len(large)) )
 
            print ("small", small)
            print ("large", large)  
            print ("\n \n") 
            phase0 = []
            phase1 = []
    f.close()
    #print ('all right') 
 


def read_result(filename, divide):
    f = open(filename, "r")
    #fout = open(fileout, "w")
    phase0 = []
    phase1 = []
    mid = int(divide)  
    count = 1  
    for line in f:
        #if not line.startswith("read"):
        #    fout.write(line)
        #else:
        if line.startswith('S'):
            if count%2 == 1:
                phase0 = line.strip().split(',')
            else:
                phase1 = line.strip().split(',')
            count += 1
        if len(phase0) >0 and len(phase1) >0:
            print ("phase0 size: ", len(phase0))
            print ("phase1 size: ", len(phase1))
            print ("check intersection: ", set(phase0).intersection(set(phase1)))
            small =[]
            large = []  
            for i in phase0:
                c = i.split('_')[1]
                if int(c) > mid:
                    large.append(i)
                else:
                    small.append(i)
            print ("in phase0, small size: %s large size: %s" % (len(small), len(large)) )
            print ("small", small)
            print ("large", large)   

            small =[]
            large = []  
            for i in phase1:
                c = i.split('_')[1]
                if int(c) > mid:
                    large.append(i)
                else:
                    small.append(i)
            print ("in phase1, small size: %s large size: %s" % (len(small), len(large)) )
 
            print ("small", small)
            print ("large", large)  
            print ("\n \n") 
            phase0 = []
            phase1 = []
    f.close()
    #print ('all right') 
 
    #fout.close()

if __name__ == "__main__":

    #input: **phasing_result readID
    #check reads be classfied correctly or not    

    # small Id and larger Id
    #read_result(sys.argv[1], sys.argv[2])

    # Id 1 and Id 2

    read_result2(sys.argv[1])



