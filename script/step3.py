#!/bin/bash
import os
import sys



#rate = [(0.29, 0.39), (0.30, 0.40), (0.31, 0.41)]
rate = [(0.37, 0.47), (0.4, 0.5)]

for (i,j) in rate:
    filename = str(i) + "_" + str(j)
    os.system("mkdir " + filename) 
    os.chdir(filename)
    command1 = "nohup python3 /home/yulin/liyanbo/script/HaploDivide/version2/main.py ../chr1.sorted.bam /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_1_1*.fa " + str(i) + " " + str(j) + " 0 20 > main.log &"
    #command1 = "nohup sh ../../check.sh 20 >check.log &"
    os.system(command1)
    os.chdir("../")

'''
print "\n"
for (i,j) in rate:
    filename = str(i) + "_" + str(j)
    os.chdir(filename)
    os.system("tail -1 snp_result")
    os.chdir("../")

print "\n"
for (i,j) in rate:
    filename = str(i) + "_" + str(j)
    os.chdir(filename)
    os.system("tail -1 delete_filter.log")
    os.chdir("../")

print "\n"
for (i,j) in rate:
    filename = str(i) + "_" + str(j)
    os.chdir(filename)
    os.system("tail -1 delete_nofilter.log")
    os.chdir("../")

print "\n"
for (i,j) in rate:
    filename = str(i) + "_" + str(j)
    os.chdir(filename)
    os.system("tail -1 insert_filter.log")
    os.chdir("../")


print "\n"
for (i,j) in rate:
    filename = str(i) + "_" + str(j)
    os.chdir(filename)
    os.system("tail -1 insert_nofilter.log")
    os.chdir("../")
'''
