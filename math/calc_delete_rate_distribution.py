import os
import sys
import tools

'''
def read_input():

    for line in f:
        words = line.split("\t")
        localCov = int(words[0])
        slope = float(words[1])
'''


#1 113456269 0.718690399145
#2 29655782 0.187854985806
#3 9828722 0.0622601835891
#4 3221710 0.0204079692223
#5 1109865 0.0070304561121
#6 318519 0.0020176632747
#7 121852 0.000771873280239
#8 43930 0.000278275228974
#9 24984 0.000158261514243
#10 16346 0.000103543976617


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        print "Usage: python " + sys.argv[0] + " (local_coverage slope)/txt "
    
    cov = int(sys.argv[1])
    deleteRate = float(sys.argv[2])
    for i in range(0,cov):
        #xyz pattern
        Rate = tools.get_combination_number(cov,i)*(deleteRate**i)*((1-deleteRate)**(cov-i))*0.72
        #print ("delete rate: %.2f  probability: %.2f" % (float(i)/cov, Rate))
        #yy pattern
    #print "yy"
    #for i in range(0,cov):
        Rate += tools.get_combination_number(2*cov,i)*(deleteRate**i)*((1-deleteRate)**(2*cov-i))*0.19
        #print ("delete rate: %.2f  probability: %.2f" % (float(i)/cov, Rate))
        #yyy pattern
    #print "yyy"
    #for i in range(0,cov):
        Rate += tools.get_combination_number(3*cov,i)*(deleteRate**i)*((1-deleteRate)**(3*cov-i))*0.06
        Rate += tools.get_combination_number(4*cov,i)*(deleteRate**i)*((1-deleteRate)**(4*cov-i))*0.02
        Rate += tools.get_combination_number(5*cov,i)*(deleteRate**i)*((1-deleteRate)**(5*cov-i))*0.007
        Rate += tools.get_combination_number(6*cov,i)*(deleteRate**i)*((1-deleteRate)**(6*cov-i))*0.002
        Rate += tools.get_combination_number(7*cov,i)*(deleteRate**i)*((1-deleteRate)**(7*cov-i))*0.0008
        Rate += tools.get_combination_number(8*cov,i)*(deleteRate**i)*((1-deleteRate)**(8*cov-i))*0.0003
        Rate += tools.get_combination_number(9*cov,i)*(deleteRate**i)*((1-deleteRate)**(9*cov-i))*0.0002
        Rate += tools.get_combination_number(10*cov,i)*(deleteRate**i)*((1-deleteRate)**(10*cov-i))*0.0001
        print ("delete rate: %.2f  probability: %.2f" % (float(i)/cov, Rate))
