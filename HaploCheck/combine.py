import os
import sys
import read
import tools


def find_common_start_end_index(s1, s2):

    #print "s1: ", s1
    #print "s2: ", s2
    a,len,s1_end,s2_end = tools.find_lcsubstr(s1, s2)
    
    if a.count(',') < 3:
        return -1,-1,-1,-1,-1,-1

    first = a.split(',')[1]
    last = a.split(',')[-2]
    #print "same content: ", a
    #print "same first position ", first
    #print "same last position ", last
    words = s1.split(',')
    index1_first = words.index(first)
    index1_last = words.index(last)

    words = s2.split(',')
    index2_first = words.index(first)
    index2_last = words.index(last)

    return index1_first, index1_last, index2_first, index2_last, s1_end, s2_end



def read_phasing_result(filename):
    
    lineNumber = 0
    f = open(filename, "r") 
    pos_binary =[] 
    for line in f:
        lineNumber += 1
        if lineNumber % 6 == 0: 
            #words = line.strip().split(',')
            pos = line.strip()
            print "pos line", pos
        elif lineNumber % 6 == 1 and lineNumber >1:
            binary = line.strip()
            print "binary line", binary
            pos_binary.append((pos, binary))

    return pos_binary

if __name__ == "__main__":  

    #filename = sys.argv[1] 
    os.system("ls *_phasing_result* > combine_file")
    files = read.read_file_list("combine_file")

    haplos = {}
    for f in files:
        pos_binary = read_phasing_result(f)
        for (pos, binary) in pos_binary:
            #print pos
            haplos[int(pos.split(',')[0])] = (pos, binary)
    sortedHaplos = haplos.items()
    sortedHaplos.sort()

    '''
    print "neibor"
    for s in sortedHaplos:
        print s[1][0]
    sys.exit()
    '''
    newPos = ""
    newBinary = ""
    print "before combine, segment number:", len(sortedHaplos)
    start_haplo = sortedHaplos[0][1] 
    count = 0
    for i in range(1, len(sortedHaplos)):
        curr_haplo = sortedHaplos[i][1]
        index1_first, index1_last, index2_first, index2_last, pos1_end, pos2_end = find_common_start_end_index(start_haplo[0], curr_haplo[0])
        if index1_first == -1:
            print "segment ", count 
            print start_haplo[0]
            print start_haplo[1]
            count += 1
            start_haplo = curr_haplo
            continue

        partBinary1 = start_haplo[1][index1_first:index1_last]
        partBinary2 = curr_haplo[1][index2_first:index2_last]
        print partBinary1
        print partBinary2
        len1 = len(partBinary1)
        #sys.exit() 
        if (partBinary1 == partBinary2):
            #if len(newPos) == 0: 
            print "same link"
            newPos = start_haplo[0][:pos1_end] + curr_haplo[0][pos2_end:]
            newBinary = start_haplo[1][:index1_last] + curr_haplo[1][index2_last:]
        elif tools.is_Bool_Reverse(partBinary1, partBinary2):
            print "reverse link"
            newPos = start_haplo[0][:pos1_end] + curr_haplo[0][pos2_end:]
            newBinary = start_haplo[1][:index1_last] + tools.bool_Reverse(curr_haplo[1][index2_last:])
        else:
            print "segment ", count 
            print start_haplo[0]
            print start_haplo[1]
            count += 1
            start_haplo = curr_haplo
            continue
        '''
        print start_haplo[0]
        print curr_haplo[0]
        print newPos
        print newBinary
        '''
        start_haplo = (newPos, newBinary)
         
    print "segment ", count 
    print start_haplo[0]
    print start_haplo[1]
