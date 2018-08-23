#!/bin/python

import os, sys

#reference = "/home/yanbo/bio/Project/Guofei/reference/chr1/chr1.fa"
reference = "/home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_1_10000000.fa"
reads ="sd_1_2.fq" 
command3 = "/home/yulin/liyanbo/Tools/bwa-0.7.17/bwa mem -t 30 " + reference + " " + reads + " >  chr1.mem.sam"
command4 = "samtools view -@ 30 -S -b chr1.mem.sam >chr1.bam"

command5 = "samtools sort -@ 30 chr1.bam -o chr1.sorted.bam"
command6 = "samtools index -@ 30 chr1.sorted.bam"

os.system(command3)
os.system(command4) 

os.system(command5)
os.system(command6)

