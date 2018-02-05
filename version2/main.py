#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import time
import contig
import tools
import column


if __name__ == '__main__':
   
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

    contigs = contig.read_Contig(sys.argv[2])
    time1 = time.clock()
    columns = column.init_Columns(bamfile, contigs[0]._seq)
    time2 = time.clock()
    print ( "init columns running %s Seconds" % (time2 - time1) )

    mutation = []
    insert = []
    delete = []
    stable = []  
    for (refPos, c) in columns.items():
        c._set_Lable()
        #sys.exit()
        if c._is_insert == 1:
            insert.append(refPos)
        
        if c._is_mutation == 1:
            mutation.append(refPos)
        
        if c._is_delete == 1:
            delete.append(refPos)
        
        if c._is_stable == 1:
            stable.append(refPos)

