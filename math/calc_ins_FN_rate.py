import os
import sys

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

if __name__ == "__main__":

    #
    if(len(sys.argv) < 2):
        print "Usage: python " + sys.argv[0] + " (local_coverage slope)/txt "
    
    f = open(sys.argv[1], "r")
    for line in f:
        words = line.split("\t")
        localCov = int(words[0])
        slope = float(words[1])
        FNRate = 0
        for i in range(0,localCov/2+1):
            if i >= localCov*slope:
                break
            #print i,FNRate
            FNRate += get_combination_number(localCov,i) / float(2**(localCov))
        print ("%.2f" % (FNRate*100))    
