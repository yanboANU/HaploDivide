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


DealLength = 5000000
Overlap = 150000

Divide = False

def loop(bamfile, contigs, a1, a2, obLen, averageCov):

    for (contigName, contig) in contigs.items():
        print ("dealing ", contigName)
        if (contig._len < DealLength) or (Divide == False):
            deal_one_part(bamfile, contig, 0, contig._len, a1, a2, obLen, averageCov)
        else:
            t2 = threading.Thread(target=loop2, name='Loop2THread', args=(bamfile,contig,a1,a2,obLen, averageCov))
            t2.start()
            t2.join()

def deal_one_part(bamfile, contig, start, end, a1, a2, obLen, averageCov):
    time1 = time.clock()
    

    columns = column.init_Columns(bamfile, contig, False, start, end, a1, a2)
    time2 = time.clock()
    print ( "init columns one part running %s Seconds" % (time2 - time1) )
     
    p = phasing.Phasing(columns, contig, obLen)
    p._pre_Process(start, end, contig._len, averageCov)
    
    time3 = time.clock()
    print ( "phasing pre_process one part running %s Seconds" % (time3 - time2) )
    print ( "pre process finish" )
    ''' 
    p._phasing()
    time4 = time.clock()
    print ( "phasing one part running %s Seconds" % (time4 - time3) )

    print ( "finish phasing" )
    p._write_result2(start, end) 

    time5 = time.clock()
    print ( "writing one part running %s Seconds" % (time5 - time4) )
    '''
def loop2(bamfile, contig, a1, a2, obLen, averageCov):
    start = 0
    
    while start < contig._len:
        end = min(start+DealLength, contig._len)
        deal_one_part(bamfile, contig, start, end, a1, a2, obLen, averageCov)
        start = start + DealLength - Overlap

    

if __name__ == '__main__':

    # bam file, contig file   
    if len(sys.argv) < 3:
        print ("python3 " + sys.argv[0] + " .bam contig/scaffold/ref a1(delete_rate) a2(insert_rate) obLen averageCov")
        sys.exit()
   
    print ("python3 " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4]+" "+ sys.argv[5] + " " + sys.argv[6])
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

    contigs = contig.read_Contig(sys.argv[2])
    print ( "finish reading contigs" )

    a1 = float(sys.argv[3])

    a2 = float(sys.argv[4]) 
    #loop(bamfile, contigs)

    obLen = int(sys.argv[5]) #defalt = 3
    averageCov = int(sys.argv[6])
    #mutiplePos = read.read_mutiplePos(sys.argv[7]) 
            
    t = threading.Thread(target=loop, name='LoopTHread', args=(bamfile,contigs,a1,a2,obLen, averageCov))
    t.start()
    t.join()
    

