import os
import sys
import math


def get_factorial(n):
    if n <= 1:
        return 1
    else:
        return n * get_factorial(n-1)


def get_combination_number(n,m):
    first = get_factorial(n)
    second = get_factorial(m)
    third = get_factorial(n-m)
    return first/(second*third)

def get_FN_rate_no_error(filename):
    f = open(filename, "r")
    for line in f:
        words = line.split("\t")
        localCov = int(words[0])
        slope = float(words[1])
        FNRate = 0
        for i in range(0,localCov/2):
            if i*slope >= localCov-i:
                break
            #print i,FNRate
            FNRate += get_combination_number(localCov,i) / float(2**(localCov-1))
        print ("%.2f" % (FNRate))


def get_FN_rate(filename, errorRate):
    f = open(filename, "r")
    for line in f:
        words = line.split("\t")
        localCov = int(words[0])
        errorNumber = int(math.ceil(localCov * errorRate))  
        slope = float(words[1])
        FNRate = 0
        #for j in range(0, errorNumber):
        '''
        for j in range(0, localCov):
            sum2 = (localCov-j)
            #print j
            for i in range(0,sum2/2+1):

                #print i, sum2-i, i*slope
                if i*slope >= sum2 - i:
                    break
                #FNRate += get_combination_number(localCov,i) / float(2**(sum2-1)) * (errorRate)**(j) *(1-errorRate)**(sum2)
                FNRate += get_combination_number(localCov,j) * get_combination_number(sum2,i) / float(2**(sum2-1)) * ((errorRate)**(j)) * ((1-errorRate)**(sum2))

                #print i,FNRate
        print ("%.2f" % (FNRate))
            #print i,FNRate
        '''

                
        TPRate = 0  
        for j in range(0, localCov):
            sum2 = (localCov-j)
            #print j
            for i in range(0,sum2):
                if i>sum2-i:
                    i = sum2-i
                #print i, sum2-i, i*slope
                if i*slope >= sum2 - i:
                    TPRate += get_combination_number(localCov,j) * get_combination_number(sum2,i) / float(2**(sum2)) * ((errorRate)**(j)) * ((1-errorRate)**(sum2))
                    #print i, TPRate
            #print FNRate
        print ("%.2f" % (TPRate))
         

if __name__ == "__main__":

    #
    if(len(sys.argv) < 2):
        print "Usage: python " + sys.argv[0] + " local_coverage slope "
    filename = sys.argv[1]
    #get_FN_rate_no_error(filename)
    errorRate = 0.1
    print "with error"
    get_FN_rate(filename, errorRate)

