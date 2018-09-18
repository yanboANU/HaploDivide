import sys


def sep_record(filename):
    f = open(filename, "r")
    fout = filename.split('_')
    fHete = open(fout[0] + '_' + fout[1] + "_hete", "w")
    fHomo = open(fout[0] + '_' + fout[1] + "_homo", "w")
    heteLen = 0
    homoLen = 0
    for line in f:
        words = line.split()
        homo = words[3] 
        deleteLen = max( len(words[1].split(',')[0]) -1, len(words[2].split(',')[0]) -1 )
        if homo == '1/1' or homo == '1|1' or homo == '1|2' or homo == '1/2':
            fHomo.write("%s" % line)
            homoLen += deleteLen
        elif homo == '1/0' or homo == '0/1' or homo == '1|0' or homo == '0|1':
            fHete.write("%s" % line)
            heteLen += deleteLen
        else:
            print line
    fHomo.write("total delete len: %d" % homoLen)
    fHete.write("totel delete len: %d" % heteLen)
    f.close()
    fHomo.close()
    fHete.close()




if __name__ == "__main__":
    if(len(sys.argv) <= 1):
        print "Usage: python *vcf"
        sys.exit()

    sep_record(sys.argv[1])    
    
           
