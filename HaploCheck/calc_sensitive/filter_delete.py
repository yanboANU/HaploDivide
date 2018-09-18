import os
import sys



def read_delete(filename):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    threshold = 2
    fout = open("delete" + str(threshold), "w")
    #requre snp store in order
    for line in f:  
        words = line.strip().split()
        if len(words) == 5: #for file mutation record
            h1Cov = int(words[2])
            h2Cov = int(words[4])
            if h1Cov + h2Cov <= 35 or h2Cov+h2Cov >=60:
                continue
            if max(h1Cov, h2Cov) < threshold*min(h1Cov, h2Cov):
                fout.write(line)
    f.close()
    fout.close()


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print ("python " + sys.argv[0] + " delete_snp")
        sys.exit()

    print ("python " + sys.argv[0] + sys.argv[1])
    read_delete(sys.argv[1])




