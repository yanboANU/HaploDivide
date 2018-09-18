#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import time, threading
import contig
import tools
import column
import phasing


if __name__ == '__main__':

    # bam file, contig file   
    if len(sys.argv) < 3:
        print ("python3 " + sys.argv[0] + ".bam contig")

    print ("python3 " + sys.argv[0] + sys.argv[1] + sys.argv[2])

    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

    contigs = contig.read_Contig(sys.argv[2])
    #output Columns, Ins, Cov 
    for (contigName, contig) in contigs.items():
        columns = column.init_Columns(bamfile, contig, True, 0, contig._len)
        for (refPos, c) in columns.items():
            c._sort_map_content()
            if len(c._map_content) == 1:
                if c._map_content[0][0] == '*':
                    delete_rate = 1
            elif c._map_content[0][0] == '*':
                delete_rate = float(len(c._map_content[0][1]))/c._coverage
            elif c._map_content[1][0] == '*':
                delete_rate = float(len(c._map_content[1][1]))/c._coverage
            else:
                delete_rate = 0

            print (delete_rate)    




