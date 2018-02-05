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

    for c in columns:
        c._is_MDS()
        #sys.exit()

