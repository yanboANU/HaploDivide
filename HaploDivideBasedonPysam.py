#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import contig


# two position stable and between them don't have insert ==> merge to ==> stable range
def pos_2_Range(stablePosition, insert):
    stableRange = []
    #for i in range(len(stablePosition)-1): 
    i = 0
    while i < len(stablePosition)-1:
        s = stablePosition[i]
        e = stablePosition[i]
        while (stablePosition[i]+1 == stablePosition[i+1]) and (stablePosition[i] not in insert):
            i = i+1
            e = stablePosition[i]
            if (i >= len(stablePosition)-1):
                break 
        i = i+1
        stableRange.append((s,e)) 
    ''' version 2
    for n in stablePosition:
        if not stableRange or n > stableRange[-1][-1] + 1:
            stableRange.append([])
        stableRange[-1][1:] = n, #      
    return [','.join(map(str,r)) for r in stableRange]
    '''
    return stableRange



#the number of reads  who support insert must larger than cov*0.8
# and consider some long insert will cause different position insert, and insert number is samll
# so insert length is considered
def filter_Insert(insert, pileupcolumns, diff):
    
    for (p, reads) in pileupcolumns:
        if p.pos in insert:
            insertLen = 0
            for (name, seq) in insert[p.pos]:
                insertLen += len(seq)
            if len(insert[p.pos]) < p.n * 0.2 and insertLen < p.n*0.2:
            #if len(insert[p.pos]) < p.n * 0.3:
                insert.pop(p.pos)
            else:
                if p.pos not in diff:
                    diff[p.pos] = []
                diff[p.pos].append('I')
                print (p.pos, "insert")
    
    return insert

def get_StableRange(samfile, contigs, insert, ll, rr):
    pileupcolumns = []
    stablePosition = []  # ref same with all reads position
    diff = {}
    for pileupcolumn in samfile.pileup(): 
        pp = pileupcolumn.pos
        # ref pos
        cov = pileupcolumn.n  
        # n is coverage
        nucleotide = contigs[0].Seq[pileupcolumn.pos]
        #print ("\n coverage at base %s = %s" % (pp , cov))  
        #print ("nucleotide at base %s = %s" % (pp, nucleotide)) # ref pos, n is coverage
         
        if (pileupcolumn.pos < ll):# and pileupcolumn.pos <2050):
            continue
        if (pileupcolumn.pos > rr):
            break
        
        count = 0
        reads = [] 
        for pileupread in pileupcolumn.pileups: # pileups is iterator, cann't go back
            reads.append(pileupread)
            if pileupread.is_del or pileupread.is_refskip:
                continue
 
            if pileupread.alignment.query_sequence[pileupread.query_position] == nucleotide:
                count += 1    
            # query position is None if is_del or is_refskip is set.
        pileupcolumns.append((pileupcolumn,reads))  # pileupcolumn's pileups during append lost
        #if count >= cov*0.8:  # 80% same with reference

        if count >= cov*0.7 or len(is_Mutation(pileupcolumn))==1:  # 80% same with reference
            stablePosition.append(pp)
        ''' 
        else:
            if pp not in diff:
                diff[pp] = []
            diff[pp].append('M') 
            print (pp, "mutation")
        '''
    '''
    insert = filter_Insert(insert, pileupcolumns, diff)   # have more than 20% insert
    print ("diff:")
    for (key, value) in sorted(diff.items()):
        print (key, value) 

    print ("filter insert" , sorted(insert.items())) # diff in python2 and python3, sorted by key
    print ("stablePosition:", stablePosition)
    '''
    stableRange = pos_2_Range(stablePosition, insert)

    print ("stableRange:", stableRange)
    
    print ("len(stableRange):", len(stableRange))
    '''
    print ("check: ",pileupcolumns[0][0].pos) 
    print ("check: ",pileupcolumns[0][0].n)
    print ("check: ",pileupcolumns[0][0].reference_name)
    print ("check: ",pileupcolumns[0][0].reference_id)
    print ("check: ",type(pileupcolumns[0][1])) 
    '''
    return stableRange, pileupcolumns, diff

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
    print (sorted(Insert.items())) # diff in python2 and python3
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
    print (p[0].pos, pos)
    assert p[0].pos == pos 
    return p[0].n



def is_Mutation(oneColumn, nucleotide):
    '''
    count = 0
    for pileupread in oneColumn.pileups: # pileups is iterator, cann't go back
        if pileupread.is_del or pileupread.is_refskip:
            continue

        if pileupread.alignment.query_sequence[pileupread.query_position] == nucleotide:
      	    count += 1    
    # query position is None if is_del or is_refskip is set.
    #  pileupcolumns.append((pileupcolumn,reads))  # pileupcolumn's pileups during append lost
    #if count >= cov*0.8:  # 80% same with reference

    if count >= cov*0.7:  # 70% same with reference
        return stablePosition.append(pp)
        return False
    '''
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
 
    sortedNuc = sorted_Map_Value(nuc)
    
    res.append( (sortedNuc[0][0], supportReadName[sortedNuc[0][0]]) )

    if (sortedNuc[0][1] < p.n*0.7 and sortedNuc[1][1] >= p[0].n*0.2): 
        print ("SNP mutation")
        res.append((sortedNuc[1][0], supportReadName[sortedNuc[1][0]]))
        return False
    
        print ("stable")
    return True



def check_Mutation(pos, pileupcolumns):
    print ("mutation check")
    recordS = pileupcolumns[0][0].pos
    p = pileupcolumns[pos-recordS]
    print (p[0].pos, pos)
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

def sorted_Map_Value(m):
    sortedM = []
    for k, v in [(k, m[k]) for k in sorted(m, key=m.get, reverse=True)]:
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
    filter_Insert(insert)
    stableRange, pileupcolumns, diff = get_StableRange(samfile,contigs,insert)
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

