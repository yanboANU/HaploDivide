import os
import sys

l = ["10x","15x","20x","25x","30x"]

#l = ["15x","20x","25x","30x"]

#command = "/home/yulin/liyanbo/Tools/PBSIM-PacBio-Simulator/src/pbsim --depth "
#para = " --length-mean 20000 --length-min 3000 --length-max 50000 --accuracy-min 0.85 --accuracy-max 0.95 --accuracy-mean 0.90 --model_qc /home/yulin/liyanbo/Tools/PBSIM-PacBio-Simulator/data/model_qc_clr ~/bio/Project/Poster/reference/LongestYest/mutation0.001/ref_and_m0.001.fasta >pbsim.log "

nohup /home/yulin/liyanbo/Tools/PBSIM-PacBio-Simulator/src/pbsim --depth 15 --length-mean 2000 --accuracy-min 0.85 --accuracy-max 0.95 --accuracy-mean 0.90 --model_qc /home/yulin/liyanbo/Tools/PBSIM-PacBio-Simulator/data/model_qc_clr ~/bio/Project/Poster/reference/Chr1/mutation0.004/ref_and_m0.004.fasta >pbsim.log &


for f in l:
    os.system("mkdir "+ f)
    os.chdir(f)
    os.system(command + f[:-1] + para)
    os.system("cat sd_0001.fastq >sd_1_2.fq")
    os.system("cat sd_0002.fastq >>sd_1_2.fq")
    os.system("rm sd_000*")
    os.chdir("../.")

