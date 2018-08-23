import os
import sys


#coverage = [5, 10, 15, 20]
coverage = [25, 30]

command = "/home/yulin/liyanbo/Tools/PBSIM-PacBio-Simulator/src/pbsim --depth "
para = " --length-mean 10000 --accuracy-min 0.85 --accuracy-max 0.95 --accuracy-mean 0.90 --model_qc /home/yulin/liyanbo/Tools/PBSIM-PacBio-Simulator/data/model_qc_clr /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/first.fasta "
for i in coverage:
    f = str(i) + "x"
    os.system("mkdir " + f)
    os.chdir(f)
    os.system(command + str(i) + para)
    os.system("cat sd_0001.fastq >sd_1_2.fq")
    os.system("cat sd_0002.fastq >>sd_1_2.fq")
    os.system("rm sd_000*")
    os.system("nohup python ../haplo_pipeline_on_bombus.py > haplo.log &")
    os.chdir("../.")


'''
#step2
for i in coverage:
    f = str(i) + "x"
    os.chdir(f)
    os.system("nohup python ../haplo_pipeline_on_bombus.py > haplo.log &")
    os.chdir("../.")
'''


'''
#step4
neighbor = [0, 1, 2]
for i in coverage:
    f1 = str(i) + "x"
    os.chdir(f1)
    for nei in neighbor:
        f2 = "neighbor" + str(nei)
        print (f1, f2)
        os.chdir(f2)
        os.system("nohup sh ../../check.sh "+ str(i*2)  +" > check.log &")
        #os.chdir("../../.")
        
        #print "snp"
        #os.system("tail -1 snp_result")
        
        #print "delete"
        #os.system("tail -1 delete_nofilter.log")
        #os.system("tail -1 delete_filter.log")
        
        #print "insert"
        #os.system("tail -1 insert_nofilter.log")

        #os.system("head -7 insert_filter.log | tail -1")
        
        os.chdir("../.")
     
    os.chdir("../.")

'''
