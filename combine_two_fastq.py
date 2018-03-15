#!/usr/bin/env python

import os
import sys
import string
from itertools import ifilter,imap

if __name__ == "__main__":
    
     
    if(len(sys.argv) <= 3):
        print "Usage: python file1 file2 combinefile"
        sys.exit()
    f1 = open(sys.argv[1], "r") 
    f2 = open(sys.argv[2], "r") 
    f3 = open(sys.argv[3], "w")
    

    count=0
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f1)):
        if (line.startswith('@')):
            count = count+1
        f3.write(line + "\n")
   
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f2)):
        if (line.startswith('@')): 
            count = count+1
            f3.write(line.split('_')[0]+"_"+str(count)+"\n")
        elif (line.startswith('+S')):
            f3.write(line.split('_')[0]+"_"+str(count)+"\n")
        else:
            f3.write(line + "\n")
    f1.close()
    f2.close()
    f3.close() 




