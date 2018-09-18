import sys
import read
import contig
import tools
 
def consistency_check(snpPosition, snpContent, seq, lens): # for delete

    pattern = [] 
    ''' # for delete
    for key in snpContent:
        pattern.append(seq[key-lens:key+lens+1]) #delete element in the middle
        if snpContent[key][1] != seq[key-1] and snpContent[key][1] != seq[key-1].upper():
            print key
            print real_snpContent[key]
            print seq[key-4:key+5]
        assert snpContent[key][1] == seq[key-1] or snpContent[key][1] == seq[key-1].upper()  
    '''   

    for key in snpContent:

        pattern.append(seq[key-lens : key-1] +  snpContent[key][1] + seq[key: key+lens]) #insert element in the middle
        if snpContent[key][0] != seq[key-1]:
            print key
            print real_snpContent[key]
            print seq[key-lens : key+lens+1]
        assert snpContent[key][0] == seq[key-1] 

    print "NIST real snp Content accord with this reference"
    return pattern 
 
def get_ref_identical_bases_pattern(ref, l):
    
    refLen = len(ref)
    pattern= []
    i=0
    while i+l<=refLen:
        c = ref[i:i+l]
        i+=1
        if c.count('N') == 0:
            pattern.append(c)
    return pattern        


def get_some_position_pattern(pp, ref):
    pattern = []
    for c in pp:
        s =  ref[c-9:c+10]
        if s.count('N') == 0:
            pattern.append(s)
    return pattern        


def get_ref_pattern(ref):
    ref = ref.upper()
    refLen = len(ref)
    lenCnt = {} #key: length value: count
    tempLen = 1
    i=0
    while i<refLen and ref[i] == "N":
        i += 1
    pre = ref[i]
    i += 1
    # must end with N, otherwise count a error
    while i<refLen:
        if ref[i] == pre:
            tempLen += 1
            i += 1
        else:
            if tempLen not in lenCnt:
                lenCnt[tempLen] = 1
            else:
                lenCnt[tempLen] += 1
            tempLen = 1
            while i<refLen and ref[i] == "N":
                i += 1
            if i<refLen:
                pre = ref[i]
                i += 1
     
    length = 0
    totalCnt = 0
    for key in lenCnt:
        length += key*lenCnt[key]
        totalCnt += lenCnt[key] 
    for key in lenCnt:    
        print key, lenCnt[key], float(lenCnt[key])/totalCnt
    print length
    assert length + ref.count("N") == refLen 

def same_aA(a,b):
    return a==b or a==b.upper() or a.upper() == b


def stat_pattern(pattern):
    
    xdy = []
    xddy = []
    xdddy = []
    xddddy = []
    xdddddy = []
    patternLen = float(len(pattern))
    mid = len(pattern[0])/2
    for c in pattern: # middle is delete position
        assert not same_aA(c[mid-1], c[mid])        
        if tools.is_same( c[mid : mid+5] ): 
            xdddddy.append(c)
        elif tools.is_same( c[mid : mid+4] ):
            xddddy.append(c)
        elif tools.is_same( c[mid : mid+3] ):
            xdddy.append(c)
        elif tools.is_same( c[mid : mid+2] ):
            xddy.append(c)
        else: 
            xdy.append(c)
    print "len xdy", len(xdy), len(xdy)/patternLen 
    print "len xddy", len(xddy), len(xddy)/patternLen 
    print "len xdddy", len(xdddy), len(xdddy)/patternLen 
    print "len xddddy", len(xddddy), len(xddddy)/patternLen 
    print "len xdddddy", len(xdddddy), len(xdddddy)/patternLen 
    knowPattern = ( len(xdy) + len(xddy) + len(xdddy) + len(xddddy) + len(xdddddy) )
    print "known pattern number:", knowPattern
    print "unknown pattern number:", patternLen - knowPattern


def stat_pattern_identical_bases(pattern):
    zero =[]
    one2 = []
    two2 = []
    one3 = []
    one4 = []
    one5 = []
    patternLen = len(pattern)
    for c in pattern:
        if tools.is_same(c):
            one5.append(c)
        elif tools.is_same(c[0:-1]) or tools.is_same(c[1:]):
            one4.append(c)
        elif tools.is_same(c[0:-2]) or tools.is_same(c[1:-1]) or tools.is_same(c[2:]):
            one3.append(c)
        elif (tools.is_same(c[0:-3]) and tools.is_same(c[-2:])) or (tools.is_same(c[1:-2]) and tools.is_same(c[-2:])) or (tools.is_same(c[0:-3]) and tools.is_same(c[-3:-1])): 
            two2.append(c)
        elif tools.is_same(c[0:-3]) or tools.is_same(c[1:-2]) or tools.is_same(c[2:-1]) or tools.is_same(c[3:]):
            one2.append(c)
        else:
            zero.append(c)
    assert len(zero) + len(one2) + len(one3) + len(two2) + len(one4) + len(one5) == patternLen      
    print "no consective identical bases", len(zero), len(zero)/float(patternLen)  
    print "one two consective identical bases", len(one2), len(one2)/float(patternLen) 
    print "two two consective identical bases", len(two2), len(two2)/float(patternLen)   
    print "one three consective identical bases", len(one3), len(one3)/float(patternLen) 
    print "one four consective identical bases", len(one4), len(one4)/float(patternLen)  
    print "one five consective identical bases", len(one5), len(one5)/float(patternLen) 
    



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print ("python " + sys.argv[0] + " real_snp contig")
        sys.exit()
   
    print ("python " + sys.argv[0] + " " + sys.argv[1]+" "+ sys.argv[2])

    real_snpPosition, real_snpContent = read.read_snp(sys.argv[1])
    contigs = contig.read_Contig(sys.argv[2]) # contig seq start with 0
    assert len(contigs) == 1
    contigName, contig = contigs.popitem()

    #get_ref_pattern(contig._seq)
    neighborLen = 10
    #pattern = get_ref_identical_bases_pattern(contig._seq.upper(), 2*neighborLen+1)
    pattern = consistency_check(real_snpPosition, real_snpContent, contig._seq.upper(), neighborLen)
    print "total pattern number:", len(pattern) 

    stat_pattern(pattern)

    #stat_pattern_identical_bases(pattern)
