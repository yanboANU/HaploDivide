#!/bin/python

import os, sys

command1 = "flye --genome-size 1.5m --threads 8 --out-dir flye --pacbio-raw"
currentPath = "/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast/" 
command2 = "/home/admin-u6260133/Tools/bwa-0.7.17/bwa index scaffolds.fasta"
command3 = "/home/admin-u6260133/Tools/bwa-0.7.17/bwa mem -t 8 ../scaffolds.fasta "
command4 = "python /home/admin-u6260133/script/HaploDivide/read_sam.py yeast.mem.sam 1_filter.sam"
command5 = "/home/admin-u6260133/Tools/samtools/samtools/samtools view -@ 8 -S -b 1_filter.sam >1_filter.bam"

command6 = "/home/admin-u6260133/Tools/samtools/samtools/samtools sort -@ 8 1_filter.bam -o 1_filter.sorted.bam"

command7 = "/home/admin-u6260133/Tools/samtools/samtools/samtools index -@ 8 1_filter.sorted.bam"
command8 = "python3 /home/admin-u6260133/script/HaploDivide/filter_bam.py ../1_filter.sorted.bam"
command9 = "/home/admin-u6260133/Tools/samtools/samtools/samtools index -@ 8 2_filter.bam"
command10 = "python3 /home/admin-u6260133/script/HaploDivide/version2/main.py 2_filter.bam ../../scaffolds.fasta >log"

for k in range(1,2):
    mNum = "mutation"+str(k)

    readspath = "longest_yeast_" + mNum
    os.chdir(mNum)
    os.chdir(readspath)
    os.system("pwd")    

    #os.system(command1 + " read_m" + str(k) + "_ref.fastq")
    os.chdir("flye")

    os.system("pwd")    
    os.system(command2)
    os.system("mkdir haplo")
    os.chdir("haplo")
     
    os.system("pwd")    
    print (command3 + "../../read_m" + str(k) + "_ref.fastq > yeast.mem.sam")
    os.system(command3 + "../../read_m" + str(k) + "_ref.fastq > yeast.mem.sam")

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
    os.chdir("/media/admin-u6260133/Data1/Project/HaploDivide/longest_yeast")
     
