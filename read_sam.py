#!/usr/bin/env python3

import sys

def filter_SAM(filename, fileout):
    f = open(filename, "r")
    fout = open(fileout, "w")
    for line in f:
        '''
        if not line.startswith("read"):
            fout.write(line)
        else:
            words = line.split()
            if words[1] == "0" or words[1] == "16":
                fout.write(line)
        ''' 
        if line.startswith("read") or line.startswith('S'):
            words = line.split()
            if words[1] == "0" or words[1] == "16":
                fout.write(line)
        elif line.startswith("@"):    
            fout.write(line)
        else:
            continue
    f.close()
    fout.close()

if __name__ == "__main__":

    filter_SAM(sys.argv[1], sys.argv[2])
