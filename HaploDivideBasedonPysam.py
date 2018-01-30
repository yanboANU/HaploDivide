#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import contig

def write_Insert(fileName, insert):
    fout = open(fileName,"w") 

    for (k,v) in sorted(insert.items()): 
       fout.write(("reference position: %d\n") % (k))
       for (l,r) in v:
           fout.write("%s %s\n" % (l,r) )
       fout.write("\n")
    fout.close()


def write_MutationOrDelete(fileName, m):
    fout = open(fileName,"w") 

    for (k,v) in sorted(m.items()): 
       fout.write(("reference position: %d\n") % (k))
       #print (v)
       #for (l,r) in v:
           #fout.write("%s %s\n" % (l,r) )
       write_Map(fout,v[0])
       write_Map(fout,v[1]) 
       fout.write("\n")
    fout.close()



def write_Map(fout, nuc):
    for (key,v) in sorted(nuc.items()):
        fout.write("%s %s\n"%(key, v)) 
    fout.write("\n")  

#'''
#def showObM(ObM, mPos):
#    for (a, b) in ObM:
#        i=a
#        seq1 = ""
#        seq2 = ""  
#        while i<=b:
#            (nuc, supportReadName) = mPos[i]
#'''
               
def show_ObSNP(ObSNP, mPos):
    readsFlag = {} 
    for i in range(len(ObSNP)):
        p = ObSNP[i] 
        (nuc, supportReadName) = mPos[p]
        supports = [] #  
        sortedNuc = sorted_Map_Value(nuc)
        a = sortedNuc[0][0]
        b = sortedNuc[1][0]
        for (k, v) in supportReadName.items():
            #print (k,v)
            if k != a and k != b:
                supports.extend(v)
        #print (supports)
        #print ("bbb")
        #sys.exit("stop")
        for readId in supportReadName[a]:
            if readId not in readsFlag:
                readsFlag[readId] = [3]*len(ObSNP)
            readsFlag[readId][i] = 0
 
        for readId in supportReadName[b]:
            if readId not in readsFlag:
                readsFlag[readId] = [3]*len(ObSNP)
            readsFlag[readId][i] = 1

        for readId in supports:
            if readId not in readsFlag:
                readsFlag[readId] = [3]*len(ObSNP)
            readsFlag[readId][i] = 2
    sortFlag = sorted_Map_Value(readsFlag, False)
    for ele in sortFlag:
        print (ele)
 
    return readsFlag       

def phasing_Reads(readsFlag, ObSNP):
    phase1 = set()
    phase2 = set()
    length = len(ObSNP)

#    for i in range(int(length/3)-1):
#        phases = {}
#        for (r, k) in readsFlag.items():
#            f = ''.join( str(j) for j in k[3*i :3*(i+1)] )
#            if f not in phases:
#                phases[f] = 0
#            phases[f] += 1
#        #print (phases)   
#        print (sorted_Map_Value(phases))
#        #break
    
    for i in range(length-1):

        readsSupport = {}
        phases = {}
        for (r, k) in readsFlag.items():
            f = ''.join( str(j) for j in k[ i :i+3] )
            if f not in phases:
                phases[f] = 0
                readsSupport[f] = []
            phases[f] += 1
            readsSupport[f].append(r)
        #print (phases)   
        
        print (ObSNP[i:i+3])
        for ele in (sorted_Map_Value(phases)):
            if ele[0][0] != '3' and ele[0][1] != '3' and ele[0][2] != '3':
                print (ele, end="")
                print (readsSupport[ele[0]])
        print ("\n")    
        #break          


def get_unStableRange(insert, mutationOrDeletePos, mergeLen, ll, rr):
   
    diffPos = set()
    for key in insert:
        if key >ll and key <rr:
            diffPos.add(key)
    for key in mutationOrDeletePos:
        if key >ll and key <rr: 
            diffPos.add(key)
    diffPos = sorted(list(diffPos))   
    checkRange = []
    i = 0
    while i < len(diffPos)-1: 
        s = diffPos[i]
        assert ( (i +1) < len(diffPos) )
        while (diffPos[i+1] <= (diffPos[i] + mergeLen)):
            i = i+1
            if i >= len(diffPos)-1:
                break
            #if i > len(diffPos):
        i = min(i, len(diffPos)-1)  
        e = diffPos[i]
        checkRange.append((s,e))
        i = i+1
    print ("position of dismatch:", diffPos)
    print ("dismatch range:", checkRange)
    print ("number of dismatch range:", len(checkRange))
    return checkRange  

# two position stable and between them don't have insert ==> merge to ==> stable range
def pos_2_Range(stablePos, insert):
    stableRange = []
    #for i in range(len(stablePos)-1): 
    i = 0
    while i < len(stablePos)-1:
        s = stablePos[i]
        e = stablePos[i]
        while (stablePos[i]+1 == stablePos[i+1]) and (stablePos[i] not in insert):
            i = i+1
            e = stablePos[i]
            if (i >= len(stablePos)-1):
                break 
        i = i+1
        stableRange.append((s,e)) 
 #   ''' 
 #   for n in stablePos:
 #       if not stableRange or n > stableRange[-1][-1] + 1:
 #           stableRange.append([])
 #       stableRange[-1][1:] = n, #      
 #   return [','.join(map(str,r)) for r in stableRange]
 #   '''
    return stableRange



#the number of reads  who support insert must larger than cov*0.8
# and consider some long insert will cause different position insert, and insert number is samll
# so insert length is considered
def filter_Insert(insert, pileupcolumns):
    
    for p in pileupcolumns:
        if p.pos in insert:
            insertLen = 0
            for (name, seq) in insert[p.pos]:
                insertLen += len(seq)
            if len(insert[p.pos]) < p.n * 0.2 and insertLen < p.n*0.5:
            #if len(insert[p.pos]) < p.n * 0.3:
                insert.pop(p.pos)
            '''
            else:
                if p.pos not in diff:
                    diff[p.pos] = []
                diff[p.pos].append('I')
                print (p.pos, "insert")
            ''' 
    return insert

def get_StableRange(samfile, contigs, insert, ll, rr):
    mutationOrDeletePos = {}
    stablePos = []  # ref same with all reads position
    pileupcolumns = []
    fout = open("statColumns_"+str(ll)+"_"+str(rr),"w") 
    for pileupcolumn in samfile.pileup(): 
        pileupcolumns.append(pileupcolumn) 
        pp = pileupcolumn.pos
        cov = pileupcolumn.n  
        #print (pp)
        #print (len(contigs[0].Seq))
        if pp > len(contigs[0].Seq):
            sys.exit()
        nucleotide = contigs[0].Seq[pp]
         
        if (pileupcolumn.pos < ll):# and pileupcolumn.pos <2050):
            continue
        if (pileupcolumn.pos > rr):
            break
        if not is_Mutation(pileupcolumn, nucleotide, mutationOrDeletePos, fout):
            stablePos.append(pp)
 
    fout.close()

    write_MutationOrDelete("statMutation_"+str(ll)+"_"+str(rr),mutationOrDeletePos) 


    insert = filter_Insert(insert, pileupcolumns)
    write_Insert("statInsert_"+str(ll)+"_"+str(rr), insert) 
    
    stableRange = pos_2_Range(stablePos, insert)

    mutationPos = sorted(mutationOrDeletePos)
    mutationRange = pos_2_Range(mutationPos ,insert)


    print ("stableRange:", stableRange) 
    print ("len(stableRange):", len(stableRange))
    fout2 = open("stableRange_"+str(ll)+"_"+str(rr), "w")
    for (l,r) in stableRange:
        fout2.write("%s %s\n" % (l,r))
    fout2.close() 

    print ("MutationRange:", mutationRange) 
    print ("len(MutationRange):", len(mutationRange))

    #ObM = get_Ob_MutationRange(mutationRange, stableRange, insert)
    ObSNP = get_Ob_SNP(mutationPos, stableRange, insert)
    ################
    #showObM(ObM,pileupcolumns)  
    #showObM(ObM, mutationOrDeletePos)  
    readsFlag = show_ObSNP(ObSNP, mutationOrDeletePos) 
    phasing_Reads(readsFlag, ObSNP)
    return stableRange, mutationOrDeletePos

#stableRange is sorted
def is_SubRange(i, j, stableRange):
    for (a,b) in stableRange:
        if i >= a and j <= b:
            return True
        if a >= j:
            return False
    return False

def get_Ob_SNP(mutationPos, stableRange, insert):
    ObSNP = [] 
    for a in mutationPos:
        if (a-1 not in insert) and (a+1 not in insert):
            if is_SubRange(a-5, a-1,stableRange) and is_SubRange(a+1,a+5, stableRange):
                ObSNP.append(a)

    print ("ObSNP:", ObSNP) 
    print ("len(ObSNP):", len(ObSNP))
    fout = open("SNPdiff","w")
    for i in range(len(ObSNP)-1):
        fout.write("%s\n" % (ObSNP[i+1]-ObSNP[i]))

    fout.close()
    return ObSNP


def get_Ob_Mutation(mutationRange, stableRange, insert):
    # [a,b] is sub-range in stableRange
    ObM = [] 
    for (a, b) in mutationRange:
        if (a-1 not in insert) and (b+1 not in insert):
            if is_SubRange(a-5, a-1,stableRange) and is_SubRange(b+1,b+5, stableRange):
                ObM.append((a,b))

    print ("ObM:", ObM) 
    print ("len(ObM):", len(ObM))

    return ObM    



def get_Insert(samfile):
   
    Insert = {} # id: reference pos key: insert content between pos, pos+1
     
    for read in samfile.fetch():
        #print (read.reference_end)
        #print (read.query_alignment_end)
        alignedPairs = read.get_aligned_pairs()
        readName = read.query_name  
        # include indel
        #print (alignedPairs)
        i = 0
        while i < len(alignedPairs):
            (readPos, referPos) = alignedPairs[i]
            if isinstance(readPos, int) and not type(referPos) == int:
                (readPos_pre, referPos_pre) = alignedPairs[i-1]
                
                if not type(referPos_pre) == int:
                    i = i+1
                    continue  
                if referPos_pre >= read.reference_end-1:
                    break   
                if referPos_pre not in Insert:
                    Insert[referPos_pre] = [] 
                insertSeq = read.query_sequence[readPos] 
                while i+1 < len(alignedPairs): 
                    i = i+1
                    (readPos, referPos) = alignedPairs[i]
                    if isinstance(readPos, int) and not type(referPos) == int:
                        insertSeq = insertSeq + read.query_sequence[readPos] 
                    else:
                        break
                Insert[referPos_pre].append((readName, insertSeq)) 
            i = i+1
        #break  
    #print (Insert)
    #print (sorted(Insert.items())) # diff in python2 and python3
    #print (Insert)
    return Insert   
        #for key, value in sorted
        #break
''' 
class Haplo(object):
    def __init__(supporReadsId, seq):


class Diplo(object):
    def __init__(h1, h2, originalReferenceStart, originalReferenceEnd):
'''
def get_Cov(pos, pileupcolumns):
    recordS = pileupcolumns[0][0].pos
    p = pileupcolumns[pos-recordS]
    #print (p[0].pos, pos)
    assert p[0].pos == pos 
    return p[0].n



def is_Mutation(oneColumn, nucleotide, mutationOrDeletePos, fout):

    pp = oneColumn.pos
    cov = oneColumn.n
    nuc = {}
    supportReadName = {}   
    for pread in oneColumn.pileups:
        if pread.is_del:
            if '*' not in nuc:
                nuc['*'] = 1
                supportReadName['*'] = [] 
            else:
                nuc['*'] += 1  
            supportReadName['*'].append(pread.alignment.query_name) 
        elif pread.is_refskip:
            print ('\tbase in read %s = .(is_refskip)' % (pread.alignment.query_name)) 
            if '.' not in nuc:
                nuc['.'] = 1
                supportReadName['.'] = [] 
            else:
                nuc['.'] += 1 
            supportReadName['.'].append(pread.alignment.query_name) 
        else:
            cc = pread.alignment.query_sequence[pread.query_position]
            if cc.upper() not in nuc:
                nuc[cc.upper()] = 1
                supportReadName[cc.upper()] = [] 
            else:
                nuc[cc.upper()] += 1  
            supportReadName[cc.upper()].append(pread.alignment.query_name)  
    fout.write(("reference position: %d coverage: %s\n") % (pp,cov))
    write_Map(fout, nuc)
    sortedNuc = sorted_Map_Value(nuc)
    #res.append( (sortedNuc[0][0], supportReadName[sortedNuc[0][0]]) )

    if (sortedNuc[0][0] == nucleotide and sortedNuc[0][1] >= cov*0.7):
        return False  
        
    if (sortedNuc[0][1] < cov*0.7 and sortedNuc[1][1] >= cov*0.3): 
        mutationOrDeletePos[pp] = (nuc, supportReadName)
        #print ("SNP mutation")
        #res.append((sortedNuc[1][0], supportReadName[sortedNuc[1][0]]))
        return True
    #unfinish, one case: reference is wrong, most >= 0.7 but not nucleotide  
    return False


'''
def check_Mutation(pos, pileupcolumns):
    print ("mutation check")
    recordS = pileupcolumns[0][0].pos
    p = pileupcolumns[pos-recordS]
    #print (p[0].pos, pos)
    assert p[0].pos == pos 
    nuc = {}
    supportReadName = {}
     
    for pread in p[1]:
        if pread.is_del:
            #print ('\tbase in read %s = *(is_del)' % (pread.alignment.query_name))
            if '*' not in nuc:
                nuc['*'] = 1
                supportReadName['*'] = [] 
            else:
                nuc['*'] += 1  

            supportReadName['*'].append(pread.alignment.query_name) 
        elif pread.is_refskip:
            print ('\tbase in read %s = .(is_refskip)' % (pread.alignment.query_name)) 
            if '.' not in nuc:
                nuc['.'] = 1
                supportReadName['.'] = [] 
            else:
                nuc['.'] += 1 
            supportReadName['.'].append(pread.alignment.query_name) 
        else:
            cc = pread.alignment.query_sequence[pread.query_position]

            #print ('\tbase in read %s = %s' % (pread.alignment.query_name, cc))
            if cc.upper() not in nuc:
                nuc[cc.upper()] = 1
                supportReadName[cc.upper()] = [] 
            else:
                nuc[cc.upper()] += 1  
            supportReadName[cc.upper()].append(pread.alignment.query_name) 
 
    #print (nuc)
    sortedNuc = sorted_Map_Value(nuc)
    #print (sortedNuc) 
    #print (supportReadName)
    res = []
    
    res.append( (sortedNuc[0][0], supportReadName[sortedNuc[0][0]]) )
    #if (sortedNuc[1][1] >= max(p[0].n*0.2, sortedNuc[0][1]*0.5)): # >0.25*cov and > 0.5*max_character

    if (sortedNuc[1][1] >= p[0].n*0.2): # >0.2*cov
        print ("SNP mutation")
        res.append((sortedNuc[1][0], supportReadName[sortedNuc[1][0]]))
    else:
        print ("stable")
    return res
'''
def sorted_Map_Value(m, R=True):
    sortedM = []
    for k, v in [(k, m[k]) for k in sorted(m, key=m.get, reverse=R)]:
        sortedM.append((k,v))  
    return sortedM  


def check_Insert(insertList):
# rule2: mostly insert number > 0.5*len(insertList) 
    ins = {}
    supportReadName = {}
    for (readId, content) in insertList:
        if content not in ins:
            ins[content] = 1
            supportReadName[content] = [] 
        else:
            ins[content] += 1
        supportReadName[content].append(readId)

    sortedIns = sorted_Map_Value(ins)
    #print (sortedIns)
    #print (supportReadName)
    res = []
    '''
    if (sortedIns[0][1] >= len(insertList)*0.5):
        res.append((sortedIns[0][0], supportReadName[sortedIns[0][0]]))
    '''

    if (len(set(supportReadName[sortedIns[0][0]])) >= max(len(insertList)*0.5,3)):  # same read support mutiple time
        res.append((sortedIns[0][0], supportReadName[sortedIns[0][0]]))

    return res
   
'''
def check_Insert(insertList, cov):
    # maybe insert content not same
    # count mostly insert number > 0.2*cov 
    ins = {}
    supportReadName = {}
    for (readId, content) in insertList:
        if content not in ins:
            ins[content] = 1
            supportReadName[content] = [] 
        else:
            ins[content] += 1
        supportReadName.append(readId)

    sortedIns = sorted_Map_Value(ins)
    #print (sortedIns)
    #print (supportReadName)
    res = []
    
    if (sortedIns[0][1] >= cov*0.2):
        print ("SNP mutation")
             
   
    return res
'''   

def update_SupportRead(supportSeq1, supportSeq2, curAB, seq1, seq2):

    if len(supportSeq1) == 0 or (len(supportSeq1.intersection(curAB[0][1])) > len(supportSeq1.intersection(curAB[1][1]))
			      and len(supportSeq2.intersection(curAB[0][1])) < len(supportSeq2.intersection(curAB[1][1]))):
        seq1 = seq1 + curAB[0][0]    
        seq2 = seq2 + curAB[1][0] 
        supportSeq1 = set(curAB[0][1])
        supportSeq2 = set(curAB[1][1]) 
    elif (len(supportSeq1.intersection(curAB[0][1])) < len(supportSeq1.intersection(curAB[1][1])) 
       and len(supportSeq2.intersection(curAB[0][1])) > len(supportSeq2.intersection(curAB[1][1]))):
        seq1 = seq1 + curAB[1][0]
        seq2 = seq2 + curAB[0][0]    
        supportSeq1 = set(curAB[1][1])
        supportSeq2 = set(curAB[0][1]) 
    else:
        print ("canot decide")
        print ("last position support:")
        print (supportSeq1)
        print (supportSeq2)
        print ("current position support:")
        print (curAB[0][1])
        print (curAB[1][1])

        print (supportSeq1.intersection(curAB[0][1]) , supportSeq1.intersection(curAB[1][1]))
        print (supportSeq2.intersection(curAB[0][1]) , supportSeq2.intersection(curAB[1][1]))
        sys.exit()   

    return supportSeq1, supportSeq2, seq1, seq2

     
def haplo_Divide(stableRange, insert, contigs, pileupcolumns, notInsert):
    print ("DC")
    #DC
    n = len(stableRange)
    
    #if n > 10:
    #    d1  = haplo_Divide(stableRange[:n/2], insert, contigs, pileupcolumns)    
    #    d2  = haplo_Divide(stableRange[-n/2:], insert, contigs, pileupcolumns)    
    #    d =  merge(d1, d2)
    #    return d
    
    originalRS = stableRange[0][0]
    originalRE = stableRange[-1][1]
    seq1 = ""
    seq2 = ""
    supportSeq1 = set() ### or using a map store every have two choices position
    supportSeq2 = set() ###
    #for (l, r) in stableRange:
    i = 0
    while i+1 < len(stableRange):
        (fs, fe) = stableRange[i]    # (fs, fe) (ss, se)
        (ss, se) = stableRange[i+1]  # first start, first end, second start, second end  
        seq1 = seq1 + contigs[0].Seq[fs:fe+1] 
        seq2 = seq2 + contigs[0].Seq[fs:fe+1]
        print ("(%s, %s) (%s, %s)" % (fs, fe, ss, se))
        if fe in insert:
            #insertStr = (insert[fe])
            #print (insert[fe]) 
            #check_Insert(insert[fe], get_Cov(fe, pileupcolumns))
            '''
            ii = check_Insert(insert[fe])
            if len(ii) > 0:

                print ("Insert happen pos: %s insert %s" % (fe, ii[0][0]))
                print ("support insert reads name: " , ii[0][1])
                print ("support not insert reads name: ", notInsert[fe] )
                ii.append(("",notInsert[fe]))
                supportSeq1, supportSeq2, seq1, seq2 = update_SupportRead(supportSeq1, supportSeq2, ii, seq1, seq2)
            '''
            print (fe, "insert")
        for j in range(fe+1, ss):
            if j in insert:
                ''' 
                #print (insert[j])
                ii = check_Insert(insert[fe])
                if len(ii) > 0:
                    print ("Insert happen pos: %s insert %s" % (j, ii[0][0]))
                    print ("support insert reads name: " , ii[0][1])
                    print ("support not insert reads name: ", notInsert[j])
                    ii.append(("",notInsert[j]))
                    supportSeq1, supportSeq2, seq1, seq2 = update_SupportRead(supportSeq1, supportSeq2, ii, seq1, seq2)
                '''  
                print (j, "insert")
            rr = check_Mutation(j, pileupcolumns) 
            if len(rr) == 1:
                seq1 = seq1 + rr[0][0]
                seq2 = seq2 + rr[0][0] 
            elif len(rr) == 2:
                '''
                seq1 = seq1 + rr[0][0]
                seq2 = seq2 + rr[1][0]  
                '''
                print ("Mutation happen pos: %s Mutation (%s, %s)" % (j, rr[0][0], rr[1][0]))
                print ("support one mutation reads name: " , rr[0][1])
                print ("support other mutation reads name: ", rr[1][1])

                supportSeq1, supportSeq2, seq1, seq2 = update_SupportRead(supportSeq1, supportSeq2, rr, seq1, seq2)
        print ("seq1: ", seq1)
        print ("seq2: ", seq2)        
        i = i+1  
        #break 

def get_Not_Insert_Read_Support(insert, pileupcolumns):
    
    notInsert = {}
    for (p, reads) in pileupcolumns:
        if p.pos in insert:
            notInsert[p.pos] = []
            rr = set()
            for (readName, content) in insert[p.pos]:
                rr.add(readName)
            for read in reads:
                if read.alignment.query_name not in rr:
                    notInsert[p.pos].append(read.alignment.query_name)
    return notInsert
 
if __name__ == "__main__":

    # sorted bam
    samfile = pysam.AlignmentFile(sys.argv[1], "rb")
   
    # contig
    contigs = contig.read_Contig(sys.argv[2])
    insert = get_Insert(samfile)
    #insert = {} 
    stableRange, mutationOrDeletePos = get_StableRange(samfile,contigs,insert,20000, 80000)
    '''
    print ("check: ",pileupcolumns[0][0].pos) 
    print ("check: ",pileupcolumns[0][0].n)
    print ("check: ",pileupcolumns[0][0].reference_name)
    print ("check: ",pileupcolumns[0][0].reference_id)
    '''
    #notInsert = get_Not_Insert_Read_Support(insert, pileupcolumns)
    #print ("check: ",pileupcolumns[0].pileups) pileups is iterator  
    #print (sorted(insert.items())) # during get_StableRange, insert change and these change keep 
    #haplo_Divide(stableRange, insert, contigs, pileupcolumns, notInsert)

