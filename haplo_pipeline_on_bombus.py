#!/bin/python

import os, sys

# command2 = bwa index
reference = "/home/yanbo/bio/Project/Poster/reference/Chr1/N_removed_chr1.fasta"
reads ="/home/yanbo/bio/Project/Poster/reads/Chr1/m0.004/readLength2k/15x/sd_1_2.fq" 
command3 = "/home/yulin/liyanbo/Tools/bwa-0.7.17/bwa mem -t 30 " + reference + " " + reads + " >  Chr1.mem.sam"
command4 = "/home/yulin/liyanbo/Tools/samtools/samtools-1.6/samtools view -@ 30 -S -b Chr1.mem.sam >Chr1.bam"

command5 = "/home/yulin/liyanbo/Tools/samtools/samtools-1.6/samtools sort -@ 30 Chr1.bam -o Chr1.sorted.bam"

command6 = "/home/yulin/liyanbo/Tools/samtools/samtools-1.6/samtools index -@ 30 Chr1.sorted.bam"
command7 = "python3 /home/yulin/liyanbo/script/HaploDivide/filter_bam.py ../Chr1.sorted.bam"
command8 = "samtools index -@ 30 2_filter.bam"
command9 = "python3 /home/yulin/liyanbo/script/HaploDivide/version2/main.py 2_filter.bam " + reference +  " 0 0 0 >log"

mutation = "/home/yanbo/bio/Project/Poster/reference/Chr1/mutation0.004/mutation_record"

command10 = "python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/calc_TP_switch_for_ref.py " + mutation + " ../*_snp_mutation* ../*_phasing_result* >TP_switch.result"

command11 = "python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/check_phasing.py ../*_phasing_result*  > reads_check.result"

command12 = "python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/get_1st_2st_frequence.py " + mutation + " ../*columns >get_12_frequency.log" 

print (command3)
os.system(command3)
os.system(command4) 
os.system("rm Chr1.mem.sam")
os.system(command5)
os.system("rm Chr1.bam")
os.system(command6)
os.system("mkdir after_filter")
os.chdir("after_filter")
os.system("pwd")    

os.system(command7)
os.system("rm ../Chr1.sorted.bam*")
os.system(command8)
os.system(command9)

os.system("mkdir check")
os.chdir("check")
os.system(command10)
os.system(command11)
os.system(command12)
     
     
