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



def loop(bamfile, contigs, a1, a2):

    for (contigName, contig) in contigs.items():
        print ("dealing ", contigName)
        time1 = time.clock()
        columns = column.init_Columns(bamfile, contig, a1, a2)
        time2 = time.clock()
        print ( "init columns running %s Seconds" % (time2 - time1) )

        p = phasing.Phasing(columns, contig)
        print ("step1")
        p._pre_Process()

        print ( "pre process" )
        p._phasing()

        print ( "finish phasing" )
        p._write_result() 

if __name__ == '__main__':

    # bam file, contig file   
    if len(sys.argv) < 3:
        print ("python3 " + sys.argv[0] + " .bam contig/scaffold/ref a1 a2")
        sys.exit()
   
    print ("python3 " + sys.argv[0] + " " + sys.argv[1]+ sys.argv[2] + " a1 a2")
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

    contigs = contig.read_Contig(sys.argv[2])
    print ( "finish reading contigs" )

    a1 = float(sys.argv[3])

    a2 = float(sys.argv[4]) 
    #loop(bamfile, contigs)
            
    t = threading.Thread(target=loop, name='LoopTHread', args=(bamfile,contigs,a1,a2,))
    t.start()
    t.join()
    

