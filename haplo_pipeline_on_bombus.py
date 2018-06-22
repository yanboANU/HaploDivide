#!/bin/python

import os, sys

#command1 = "flye --genome-size 1.5m --threads 8 --out-dir flye --pacbio-raw"
currentPath = "/home/yulin/liyanbo/Project/A_ten/flye/phasing/one_contig " 
command2 = "bwa index  " + sys.argv[1] # scaffolds.fasta"
command3 = "bwa mem -t 30 " + sys.argv[1] + " /home/yulin/bio/pacbio_aten/newId/Long5k/Long5k_reads.fasta >  coral.mem.sam"
command4 = "python /home/yulin/liyanbo/script/HaploDivide/read_sam.py coral.mem.sam 1_filter.sam"
command5 = "samtools view -@ 30 -S -b 1_filter.sam >1_filter.bam"

command6 = "samtools sort -@ 30 1_filter.bam -o 1_filter.sorted.bam"

command7 = "samtools index -@ 30 1_filter.sorted.bam"
command8 = "python3 /home/yulin/liyanbo/script/HaploDivide/filter_bam.py ../1_filter.sorted.bam"
command9 = "samtools index -@ 30 2_filter.bam"
command10 = "python3 /home/yulin/liyanbo/script/HaploDivide/version2/main.py 2_filter.bam ../" + sys.argv[1] +  " >log"


for k in range(1,2):
    '''
    mNum = "mutation"+str(k)
    readspath = "longest_yeast_" + mNum
    os.chdir(mNum)
    os.chdir(readspath)
    os.system("pwd")    
    os.system(command1 + " read_m" + str(k) + "_ref.fastq")
    os.chdir("flye")
    ''' 
    os.system("pwd")    
    os.system(command2)
     
    print (command3)
    os.system(command3)

    os.system(command4) 
    os.system(command5)
    os.system(command6) 
    os.system(command7)


    os.system("mkdir after_filter")
    os.chdir("after_filter")

    os.system("pwd")    
    os.system(command8)


    os.system(command9)

    os.system(command10)
     
