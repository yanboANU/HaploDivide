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



def loop(bamfile, contigs):

    for (contigName, contig) in contigs.items():
        print ("dealing ", contigName)
        time1 = time.clock()
        columns = column.init_Columns(bamfile, contig)
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
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

    contigs = contig.read_Contig(sys.argv[2])
    print ( "finish reading contigs" )

     
    #loop(bamfile, contigs)
        
    t = threading.Thread(target=loop, name='LoopTHread', args=(bamfile,contigs,))
    t.start()
    t.join()
     

