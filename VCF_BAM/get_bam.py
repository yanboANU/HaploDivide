#!/bin/python

'''
for i in range(501,624):
    command1 = "samtools view -@ 30 -S -b c1_" + str(i) + ".mem.sam > alignment_bam/c1_" + str(i) + ".bam"
    print command1
    command2 = "rm c1_" + str(i) + ".mem.sam"
    print command2
    command3 = "samtools sort -@ 30 alignment_bam/c1_" + str(i) +".bam >alignment_bam/c1_"+str(i)+"_sorted.bam"
    print command3
    command4 = "rm alignment_bam/c1_" + str(i) +".bam"
    print command4
'''

'''
command3 = "samtools merge 601_600.sorted.bam "
command4 = "rm "
for i in range(50,60):
    outbam = "c1_" + str(i*10+1) + "_" +str((i+1)*10) + ".bam "
    command1 = "samtools merge " + outbam
    command2 = "rm "
    for j in range(i*10+1,(i+1)*10+1):
        command1 = command1 + "c1_"+ str(j) + "_sorted.bam "
        command2 = command2 + "c1_"+ str(j) + "_sorted.bam "
    print command1
    print command2
    command3 = command3 + outbam + " "
    command4 = command4 + outbam + " "

print command3
print command4
'''

command1 = "samtools merge 601_623.sorted.bam "
command2 = "rm "
for i in range(601, 624):
    command1 = command1 + "c1_"+ str(i) + "_sorted.bam "
    command2 = command2 + "c1_"+ str(i) + "_sorted.bam "
print command1
print command2
