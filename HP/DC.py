#########################################################################
# File Name: DC.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Wed 24 Oct 2018 15:38:51 AEDT
#########################################################################
#!/bin/bash

import sys
import tools
import copy
import numpy as np
import operator
    

from libprism.local.prepare import verify_complexity, blocks_from_clouds, contigs_from_clouds
import cloud
import haplotype


def print_aligned_region(start, end, clouds_at_index):

    print("aligned detail region: %d %d " % (start, end))
    cs = cloud.get_clouds(start, end, clouds_at_index)

    aligned_detail = []
    for c in cs:
        aligned_detail.append(c[start:end+1]) 

    for a in sorted(aligned_detail):
        pattern = -1
        rep = ['-' if x== pattern else x for x in a]
        #fout.write("%s %s\n" % (c.start, c.end))
        print("%s" % ''.join(map(str, rep)))     
        #fout.write("%s\n" % ''.join( map(str,rep) ))
        

class DC(object):


    def find_position_no_reads_span(self, start, end, clouds_at_index):
        assert end > start
        pos = []
        for i in range(start, end-1):
            j=i+1
            ij = clouds_at_index[i].intersection(clouds_at_index[j])
            #print i,j, "were spaned by ", len(ij), "reads"
            if len(ij) == 0:
                pos.append(i)
        return pos        

    def generate_blocks(self, clouds_at_index):

        blocks = list() 
        indexLen = len(clouds_at_index) - 1
        pos = self.find_position_no_reads_span(0, indexLen, clouds_at_index)
        #sys.exit()
        l = len(pos)
        for i in range(l-1):
            j=i+1
            if pos[j] -pos[i] > 1:
                blocks.append( (pos[i]+1, pos[j] ) )
        if indexLen > pos[-1]:
            blocks.append( (pos[-1]+1, indexLen) )
        return blocks        

    #def __init__(self, k, data_start, data_end, all_clouds, clouds_at_index, hmm): # all_clouds, a cloud list
    def __init__(self, k, data_start, data_end, all_clouds, clouds_at_index, contigs_at_index): # all_clouds, a cloud list
        self.length = data_end - data_start + 1
        self.data_start = data_start # snp start no
        self.data_end = data_end     # snp end no     
        self.matrix = list()
        self.k = int(k)
        self.y_states = ((0,1), (1,0))
        self.k_mers = {} # key is snp no, value is k_mer # k_mer is a map, key is tuple (0,-1,0) value is a list [readID1, readID2]
        self.k_mers_01 = {} # key is snp no, value tuple list
        # cloud is objective, libprism/local/cloud.py, cloud.name, cloud.seq(100-110), cloud[i], cloud[i:i+3]
        # cloud.start, cloud.end
        # clouds_at_index, dict, key is snp No., value is cloud who cover this SNP by 0/1 
        # contigs_at_index, dict, key is snp np, value is cloud ID who span this SNP 0/1/-, ???
        #contigs_at_index = contigs_from_clouds(all_clouds, self.data_start, self.data_end)  # bug, because a read multiple flag 0/16
        #self.hmm = hmm
        self.enumerate_length_threshold = 5 # trade-off between accuracy and time 

        ''' 
        print contigs_at_index[3220], len(clouds_at_index[3220])

        for c in clouds_at_index[3220]:
            print c.name
            
        print contigs_at_index[3221], len(clouds_at_index[3221])
        print contigs_at_index[3222], len(clouds_at_index[3222])
        print contigs_at_index[3223], len(clouds_at_index[3223])
        sys.exit()
        ''' 


        # optional   
        '''
        good_snps = self.pick_snp(clouds_at_index)
        self.set_k_mers_only_for_good_snps(clouds_at_index, good_snps)
        '''
        
        self.set_k_mers(clouds_at_index)
        #self.blocks = blocks_from_clouds(all_clouds, clouds_at_index) # probhap generate block
        #self.blocks = self.generate_blocks(clouds_at_index)
        self.blocks = self.generate_blocks(contigs_at_index)

        fout = open("_block2", "w") 
        fout.write( "block number: %d\n" % len(self.blocks) )
        for block in self.blocks:
            fout.write("%s %s\n" % (block[0], block[-1]) )


        '''
        uncovered_positions = set()
        # detect positions that have no coverage:
        for j in range(self.data_start, self.data_end+1):
            if len(clouds_at_index[j]) == 0:
                uncovered_positions.add(j)

        self.uncovered_positions = uncovered_positions
        print "uncovered number", len(self.uncovered_positions)
        print sorted(self.uncovered_positions)
        sys.exit()
        '''

    def enumerate_seq(self, start, end, clouds_at_index):

        print "enumerate seq, length:", end - start + 1 
        minMEC = 10000
        bestHaps = []
        dis = end - start + 1

        if dis >= self.enumerate_length_threshold: # need improve
            return haps
        allPosSeq = tools.enumerate_01_list(dis)    
        l = len(allPosSeq)/2
        haps = list()
         
        for seq in allPosSeq[0:l]:
            h = haplotype.Haplotype(start)        
            h.set_seq(seq)
            MEC = self.updata_clouds_support_and_calculate_MEC(h, clouds_at_index)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = []
                bestHaps.append(h)
            elif MEC == minMEC:
                bestHaps.append(h)
        l = len(bestHaps)
        if l == 1:
            print "enumerate block susscess"
            return bestHaps  # clouds keep 
        else:
            print "enumerate seq more than 1"
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq, 
                seqs.append(h.seq)
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) 
               h.assign_clouds(clouds_at_index) # updata clouds
               if len(a) - a.count(-1) >= 2:  
                   haps.append(h)
            return haps   

    def enumerate_block_prefix(self, start, firstH, clouds_at_index):
        
        dis = firstH.start - start
        haps = list()
        if dis > self.enumerate_length_threshold or dis <= 0:
            haps.append(firstH)
            return haps
        allPosSeq = tools.enumerate_01_list(dis)    
        minMEC = 10000
        for seq in allPosSeq:
            tempSeq = seq + firstH.seq
            newH = haplotype.Haplotype(start)        
            newH.set_seq(tempSeq)
            MEC = self.updata_clouds_support_and_calculate_MEC(newH, clouds_at_index)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = []
                bestHaps.append(newH)
            elif MEC == minMEC:
                bestHaps.append(newH)
        l = len(bestHaps)
        if l == 1:
            print "enumerate prefix susscess"
            return bestHaps # clouds updata in updata_clouds_support_and_calculate_MEC
        else:
            print "enumerate prefix seq more than 1"
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq
                seqs.append(h.seq)
            #sys.exit()
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            #print segments
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) - a.count(-1) >= 2:  
                   haps.append(h)
            return haps   


    def enumerate_block_suffix(self, lastH, end, clouds_at_index):
        dis = end - lastH.end
        haps = list()
        if dis > self.enumerate_length_threshold or dis <= 0:
            haps.append(lastH)
            return haps
        start = lastH.start
        allPosSeq = tools.enumerate_01_list(dis)   

        minMEC = 10000
        for seq in allPosSeq:
            tempSeq = lastH.seq + seq
            newH = haplotype.Haplotype(start)        
            newH.set_seq(tempSeq)
            MEC = self.updata_clouds_support_and_calculate_MEC(newH, clouds_at_index)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = []
                bestHaps.append(newH)
            elif MEC == minMEC:
                bestHaps.append(newH)
        l = len(bestHaps)
        if l == 1:
            print "enumerate suffix susscess"
            return bestHaps
        else:
            print "enumerate suffix seq more than 1"
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq
                seqs.append(h.seq)
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            #print segments
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) - a.count(-1) >= 2:  
                   haps.append(h)
            return haps   



    '''
    def dp(self, haps):
        addHaps = []
        l = len(haps)
        for i in range(l-1):

            h1 = haps[i] 
            h2 = haps[i+1]
            start = h1.end + 1 
            end = h2.start - 1
            if end <= start + 10:
                continue
            b = range(start, end+1) 
            h = haplotype.Haplotye(start)
            s = run_viterbi(b)
    '''        

    def generate_distribution(self, haps): # distance distribution
        
        l = len(haps)
        dis = {}
        unphased = 0
        for i in range(l-1):
            h1 = haps[i] 
            h2 = haps[i+1]
            d = h2.start - 1 - h1.end
            if d >= 10:
                unphased += d
            if d not in dis:
                dis[d] = 1
            else:
                dis[d] += 1
        fout = open("distance_distribution", "w")
        for (d, count) in sorted(dis.items()):
            fout.write("%s %s\n" % (d, count))

        fout.write("%s %.2f\n" % (unphased, float(unphased)/self.length*100))
        fout.close()


 
    def run_dc_in_block(self, start, end, clouds_at_index): 
        
        label = 2
        haps = list()
        while (len(haps) == 0 and label >= 1):
            print "check bit by label", label
            haps = self.generate_haps_from_k_mers(start, end, label, clouds_at_index)
            label -= 1 # loose bit-complement constraint
 
        print "after generate haps from kmers, haps number ", len(haps)
        hapsNum = len(haps)
        if hapsNum == 0:
            haps = self.enumerate_seq(start, end, clouds_at_index)
        elif hapsNum >= 1:
            haps = self.generate_newhaps_loose_bit_constraint_locally(haps, start, end, label, clouds_at_index)

        #enumerate prefix and subfix
        if len(haps) >= 1:
            #for h in haps:
                #print h.start, h.end, h.seq
            prefixHaps = self.enumerate_block_prefix(start, haps[0], clouds_at_index) 
            haps.remove(haps[0])
            haps = prefixHaps + haps
            suffixHaps = self.enumerate_block_suffix(haps[-1], end, clouds_at_index)
            haps.remove(haps[-1])
            haps = haps + suffixHaps
            #print "after enumerate prefix and subfix"
            #for h in haps:
                #print h.start, h.end, h.seq 

        #enumerate gaps        
        if len(haps) > 1: 
            a = 1000000
            count = 1
            while len(haps) < a:
                a = len(haps)
                print "merge count", count
                haps = self.enumerate_gap_among_haps2(haps, clouds_at_index)
                #haps = self.enumerate_gap_between_haps(haps, clouds_at_index)
                count += 1
                 
            #haps = self.enumerate_gap_between_haps(haps, clouds_at_index)
            #haps = self.enumerate_gap_among_haps2(haps, clouds_at_index)  
            #haps = self.enumerate_gap_between_haps(haps, clouds_at_index)
            #haps = self.enumerate_gap_among_haps2(haps, clouds_at_index)
            #haps = self.enumerate_gap_between_haps(haps, clouds_at_index)
        print "block region", start, end, "has", len(haps), "haplotype segements"   
        return haps

    def generate_newhaps_loose_bit_constraint_locally(self, haps, start, end, label, clouds_at_index):
        
        print "old haps in subfunction, number", len(haps)
        subHaps = {}
        hapsNum = len(haps)
        for i in range(hapsNum+1): 
            s = start if i==0 else haps[i-1].end + 1   
            e = end if i==hapsNum else haps[i].start - 1
            subLabel = label
            hs = list()
            if (e - s + 1 >= self.enumerate_length_threshold) and subLabel >= 1: # e>s, error rate high, 
                while (len(hs) == 0 and subLabel >= 1):
                    hs = self.generate_haps_from_k_mers(s, e, subLabel, clouds_at_index)
                    subLabel -= 1
                if len(hs) > 0:
                    subHaps[i] = hs
        
        if len(subHaps) > 0:    
            newHaps = list()        
            for i in range(hapsNum+1):
                if i in subHaps:
                    newHaps.extend(subHaps[i])
                if i < hapsNum:
                    newHaps.append(haps[i]) 
            print "new haps in subfunction, number", len(newHaps)
            haps = newHaps        
        return haps


        

    def run_dc(self, clouds_at_index): # 12, Dec

        #print self.find_position_no_reads_span(0, len(clouds_at_index), clouds_at_index)
        #for block in self.blocks:
            #print block
        #sys.exit() 
        #print_aligned_region(8218, 8230, clouds_at_index)
        #print_aligned_region(17885, 17895, clouds_at_index)
        #sys.exit()
        
        allHaps = []
        count = 1
        print "There are", len(self.blocks), "blocks"
        for block in self.blocks:  # blo/cks break at no reads cover region
            print "phasing block: ", count, ", region ", block[0] , block[-1]
            allHaps.extend(self.run_dc_in_block(block[0], block[-1], clouds_at_index))
            count += 1

        #self.generate_distribution(allHaps)
        #sys.exit()
        #print "before finall enmuerate ", len(allHaps)    
        #allHaps = self.enumerate_gap_between_haps(allHaps, clouds_at_index) # only need work for between block, need to updata
                                                                  # for two consective hap, local region , two MEC
                                                                  # longer , maybe have MEC
                                                                  
        '''
        fout1 = open("hap.out" , "w")
        for i, block in enumerate(allHaps):
            fout1.write( "%d %d\n" % (block.start, block.end) )
            for k,y in enumerate(block.seq):
                j = block.start + k
                if y != -1:
                    fout1.write("%d\t%d\t%d\t\n" % (j, y, tools.int_reverse(y)))
            fout1.write("********\n")
        fout1.close()
        ''' 
        return allHaps  





    def enumerate_gap_between_haps(self , haps, clouds_at_index):

        MEC = 0
        l = len(haps)
        fout = open("enumerate_fail" ,"w")
        s = 0
        while True:
            if len(haps) == 1 or s+1 >= len(haps):
                break
            h1 = haps[s]
            h2 = haps[s+1]
            haps.remove(h1)
            haps.remove(h2)

            if h1.end > h2.start -1:  # have overlap

                print "check intersection"
                newHaps = self.from_two_haps(h1, h2, clouds_at_index) # satisfy reads intersection and overlap 
            elif h1.end == h2.start -1: # no gap, no overlap    
                newHaps = self.from_two_haps(h1, h2, clouds_at_index)
                if len(newHaps) == 2:
                    newHaps = self.enumerate_two_haps(h1, h2, clouds_at_index, fout)
                #if len(newHaps1) < len(newHaps2):
                    #newHaps = newHaps1
                #else:
                    #newHaps = newHaps2
            else: # have gap 
                print "enumrate gap"
                newHaps = self.enumerate_two_haps(h1, h2, clouds_at_index, fout)
      
            haps = haps[0:s] + newHaps + haps[s:]
            s += len(newHaps) - 1

        print "after enumerate gap bwtween haps, number of haplotype", len(haps)
        for h in haps:
            print h.start, h.end, h.seq
        return haps


    # enumerate gap between two haps, not optimal solution
    # dp is optimal solution
    def enumerate_gap_among_haps(self, haps, clouds_at_index): # every time include 3 haps  
       
        MEC = 0
        #l = len(haps)
        s = 0
        for h in haps:
            h.get_clouds(clouds_at_index)
        while True:    
            if len(haps) < 3 or s+2 >= len(haps):
                break
            
            newHaps = list()
            h1 = haps[s]
            h2 = haps[s+1]
            h3 = haps[s+2]

            intersectionNum = len( h1.clouds.intersection(h3.clouds) )
            
            print h1.start, h1.end, h3.start, h3.end, "intersectionNum:", intersectionNum
            if intersectionNum > 2: # only one reads maybe error
                if h2.start > h1.end and h3.start > h2.end:
                    dis = h2.start - h1.end + h3.start - h2.end - 2 + 2
                    if dis < self.enumerate_length_threshold:
                        print "dis:", dis
                        haps.remove(h1)
                        haps.remove(h2)
                        haps.remove(h3)
                        newHaps = self.enumerate_three_haps(h1, h2, h3, clouds_at_index)
                        haps = haps[0:s] + newHaps + haps[s:]
                        if len(newHaps) == 1:
                            s += len(newHaps) - 1          
                        elif len(newHaps) > 1:
                            s += len(newHaps) - 2          
            if len(newHaps) == 0:                    
                s += 1        

        return haps      


               
    def enumerate_gap_among_haps2(self, haps, clouds_at_index): # every time may include more than 2 haps  
       
        MEC = 0
        s = 0
        for h in haps:
            h.get_clouds(clouds_at_index)
        print "debug, orignal haps number", len(haps)    
        while True:    
            if len(haps) < 2 or s+1 >= len(haps):
                break 
            newHaps = list()
            waitMerge = list()
            dis = haps[s+1].start - haps[s].end
            print haps[s].start, haps[s].end, haps[s+1].start, haps[s+1].end
            assert haps[s].end < haps[s+1].start
            intersectionNum = len( haps[s].clouds.intersection(haps[s+1].clouds) )
            #if intersectionNum > 0:
            waitMerge.append(haps[s])
            waitMerge.append(haps[s+1])
            i = 2
            while dis < self.enumerate_length_threshold:
                if s+i < len(haps):
                    intersectionNum = len( haps[s].clouds.intersection(haps[s+i].clouds) )
                else:
                    break
                if intersectionNum < 1:
                    break
                if i > 3:
                    break # at most enumerate three haps
                assert haps[s+i].start > haps[s+i-1].end
                dis = dis + haps[s+i].start - haps[s+i-1].end
                if dis > self.enumerate_length_threshold:
                    dis = dis - (haps[s+i].start - haps[s+i-1].end)
                    break
                waitMerge.append(haps[s+i])
                print haps[s+i].start, haps[s+i].end
                i += 1
            print "debug merge haps number ", len(waitMerge), "intersectionNum:", intersectionNum, "dis",  dis
            
            for h in waitMerge:
                haps.remove(h)
            newHaps = self.enumerate_haps(waitMerge, clouds_at_index)
            haps = haps[0:s] + newHaps + haps[s:]
            if len(newHaps) == 1:
                s += len(newHaps) - 1          
            elif len(newHaps) > 1:
                s += 1          
            #if len(newHaps) == 0:                    
                #s += 1        

        return haps      
        

    ''' 
    def enumerate_clouds(start, end, cloudsMid):
        leftClouds = []
        rightClouds = []
        cloudsNum = len(cloudsMid)
        for i in range(0, pow(2, cloudsNum-1)):
            # unfinish
            #
    '''       

    def enumerate_two_haps(self, h1, h2, clouds_at_index, fout):

        minMEC = 1000000
        hapsFromEnumerateTwoHaps = list()
        bestHaps = []
        dis = h2.start - h1.end - 1
        print h1.start, h1.end, h2.start, h2.end
        print "enumerate gap between two haps, dis:", dis 
        start = h1.start 
        #if h1.end >= h2.start: # only work for k>=3
       
        assert h1.end <= h2.start-1 # no gap to enumerate
        
        '''    
        if h1.end >= h2.start-1: # no gap to enumerate
            print "two segments have overlap"
            return self.from_two_haps(h1, h2, clouds_at_index) # satisfy reads intersection and overlap 
                                              # change return 
        '''    
        if dis >= self.enumerate_length_threshold: # need improve
            print "dis is", dis, "too long to enumerate"   # (1) enumerate cloud (2) run viterbi/dp
            cloudsMid = list(cloud.get_clouds(h1.end+1, h2.start-1, clouds_at_index))
            cloudsNum = len(cloudsMid)
            print cloudsNum 
            '''
            b = range(h1.end+1, h2.start-1)
            h = self.hmm.run_viterbi(b)
            '''
            hapsFromEnumerateTwoHaps.append(h1)
            hapsFromEnumerateTwoHaps.append(h2)
            return hapsFromEnumerateTwoHaps
        
        allPosSeq = tools.enumerate_01_list(dis)    
        for seq in allPosSeq:
            for i in range(0,2): # debug maybe h1.seq + seq + tools.reverse(h2.seq)
                if i == 0:
                    temp_seq = h1.seq + seq + h2.seq 
                elif i == 1:
                    temp_seq = h1.seq + seq + tools.list_reverse(h2.seq) # debug maybe h1.seq + seq + tools.reverse(h2.seq)

                h3 = haplotype.Haplotype(h1.start)
                h3.set_seq(temp_seq)
                MEC = self.updata_clouds_support_and_calculate_MEC(h3, clouds_at_index)
                if MEC < minMEC:
                    minMEC = MEC
                    bestHaps = []
                    bestHaps.append(h3)
                elif MEC == minMEC:
                    bestHaps.append(h3)
        #print "min MEC", minMEC
        #print "best seq", best_seq

        print "h1.MEC + h2.MEC", h1.MEC + h2.MEC
        l = len(bestHaps)
        if l == 0:
            hapsFromEnumerateTwoHaps.append(h1)
            hapsFromEnumerateTwoHaps.append(h2)
            return hapsFromEnumerateTwoHaps
        if l == 1:
            print "only one minimal MEC, fill and merge"
            #h1 = copy.deepcopy(bestHaps[0])  # not right way to use
            print bestHaps[0].MEC
            #bestHaps[0].assign_clouds(clouds_at_index) #no need
            hapsFromEnumerateTwoHaps.append( bestHaps[0] )
            return hapsFromEnumerateTwoHaps
        elif l >= 2:   # can check with reads intersection
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq, h.MEC
                seqs.append(h.seq)

            print "ll"
            for c in h1.left_clouds.intersection(h2.left_clouds):
                print c.name, c.seq 
            print "lr"    
            for c in h1.left_clouds.intersection(h2.right_clouds):
                print c.name, c.seq
            
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            print segments
            haps = list()
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) -a.count(-1) >= 2:
                   hapsFromEnumerateTwoHaps.append(h)
            return hapsFromEnumerateTwoHaps 
        else:
            print "some error happen", l, minMEC


    def generate_merge_seq(self, hs, label, originalLen, temp_seq, ans):
        if label == originalLen:
            ans.append(temp_seq)
            return
        
        dis = hs[label].start - hs[label-1].end - 1
        allPosSeq = tools.enumerate_01_list(dis)  
        for seq in allPosSeq:
            for j in range(0,2):
                if j == 0:
                    self.generate_merge_seq(hs, label+1, originalLen, temp_seq + seq + hs[label].seq, ans)
                if j == 1:
                    self.generate_merge_seq(hs, label+1, originalLen, temp_seq + seq + tools.list_reverse(hs[label].seq), ans)
        return            

    



    def enumerate_three_haps(self, h1, h2, h3, clouds_at_index):
        
        minMEC = 1000000
        hapsFrom3Haps = list()
        bestHaps = []
        dis12 = h2.start - h1.end - 1
        dis23 = h3.start - h2.end - 1
        
        start = h1.start 

        print "enumerate 3 haps", h1.start, h1.end, h2.start, h2.end, h3.start, h3.end 
        allPosSeq12 = tools.enumerate_01_list(dis12)    
        allPosSeq23 = tools.enumerate_01_list(dis23)    
        for seq12 in allPosSeq12:
            for seq23 in allPosSeq23:
                for i in range(0,4): # debug maybe h1.seq + seq + tools.reverse(h2.seq)
                    if i == 0:
                        temp_seq = h1.seq + seq12 + h2.seq + seq23 + h3.seq
                    elif i == 1:
                        temp_seq = h1.seq + seq12 + tools.list_reverse(h2.seq) + seq23 + h3.seq
                    elif i == 2:
                        temp_seq = h1.seq + seq12 + h2.seq + seq23 + tools.list_reverse(h3.seq)
                    elif i == 3:
                        temp_seq = h1.seq + seq12 + tools.list_reverse(h2.seq) + seq23 + tools.list_reverse(h3.seq)
                    #print "debug", temp_seq  
                    h4 = haplotype.Haplotype(h1.start)
                    h4.set_seq(temp_seq)
                    MEC = self.updata_clouds_support_and_calculate_MEC(h4, clouds_at_index)
                    if MEC < minMEC:
                        minMEC = MEC
                        bestHaps = []
                        bestHaps.append(h4)
                    elif MEC == minMEC:
                        bestHaps.append(h4)

        l = len(bestHaps)
        if l == 0:
            hapsFrom3Haps.append(h1)
            hapsFrom3Haps.append(h2)
            hapsFrom3Haps.append(h3)
            return hapsFrom3Haps
        if l == 1:
            print "only one minimal MEC, fill and merge", minMEC
            #bestHaps[0].assign_clouds(clouds_at_index)
            hapsFrom3Haps.append( bestHaps[0] )
            #sys.exit()
            return hapsFrom3Haps

        elif l >= 2:
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq
                seqs.append(h.seq)
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            print segments
            haps = list()
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) -a.count(-1) >= 2:
                   hapsFrom3Haps.append(h)
            return hapsFrom3Haps 
        else:
            print "some error happen", l, minMEC


    def enumerate_haps(self, hs, clouds_at_index):
       
        print "enumerate haps number",len(hs) 
        minMEC = 1000000
        hapsFromHs = list()
        bestHaps = []
        start = hs[0].start 
        originalLen = len(hs)    
        assert originalLen >= 2
        temp_seq = hs[0].seq  
        ans = list()
        for h in hs:
            print h.start, h.end, h.seq
        self.generate_merge_seq(hs, 1, originalLen, temp_seq, ans)
        for temp_seq in ans:
            #print "debug", temp_seq  
            h4 = haplotype.Haplotype(hs[0].start)
            h4.set_seq(temp_seq)
            MEC = self.updata_clouds_support_and_calculate_MEC(h4, clouds_at_index)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = []
                bestHaps.append(h4)
            elif MEC == minMEC:
                bestHaps.append(h4)

        l = len(bestHaps)
        if l == 0:
            return hs
        if l == 1:
            print "only one minimal MEC, fill and merge", minMEC
            hapsFromHs.append( bestHaps[0] )
            return hapsFromHs
        elif l >= 2:
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq, h.MEC
                seqs.append(h.seq)
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            print segments
            haps = list()
            l = len(segments)
            if l >= 2:
                print "enumerate haps more than 1"
                hapsFromHs = self.deal_with_same_MEC(hs[0], hs[-1], bestHaps, clouds_at_index)
                if len(hapsFromHs) == 1:
                    print "same MEC success"
                    return hapsFromHs
                else:
                    hapsFromHs = list()
                #return self.deal_with_same_MEC(hs[0], hs[-1], clouds_at_index)
                '''
                print "ll"
                hapsFromHs = self.from_two_haps(hs[0], hs[-1], clouds_at_index)
                if len(hapsFromHs) == 1:
                    return hapsFromHs
                else:
                    hapsFromHs = list()

                cloudsSupportTwoPos = clouds_at_index[hs[0].end].intersection(clouds_at_index[hs[-1].start])
                ss = ""
                for c in cloudsSupportTwoPos:
                    print "support two pos reads", c.name, c.seq
                    print hs[0].end, hs[-1].start
                    ss = c[hs[0].end : hs[-1].start+1]
                    print "ss", ss
                #sys.exit()
            
                # merge haplotype according to cloudsSupportTwoPos
                if len(ss) >0 and hs[0].end + 1 == hs[-1].start:
                    if (ss[0] == hs[0][-1]  and ss[-1]== hs[-1][0]) or (ss[0] != hs[0][-1]  and ss[-1] != hs[-1][0]):
                        temp_seq = hs[0].seq + hs[-1].seq
                    else:
                        temp_seq = hs[0].seq + tools.list_reverse(hs[-1].seq)
                        
                    h = haplotype.Haplotype(hs[0].start)        
                    h.set_seq(temp_seq) # cloud not updata
                    h.assign_clouds(clouds_at_index)
                    hapsFromHs.append(h)
                    print "merge haplotype according to cloudsSupportTwoPos"
                    return hapsFromHs
                '''    
            print "again into enumerate"
            #sys.exit()
            for i in range(l):
                a = [ -1 if e==-2 else e for e in segments[i][1] ]
                h = haplotype.Haplotype(segments[i][0])        
                h.set_seq(a) # cloud not updata
                h.assign_clouds(clouds_at_index)
                if len(a) -a.count(-1) >= 2:
                    hapsFromHs.append(h)
            return hapsFromHs 
        else:
            print "some error happen", l, minMEC

    def deal_with_same_MEC(self, h1, h2, bestHaps, clouds_at_index):
        
        print "deal with same MEC"
         
        print h1.end, h2.start
        cloudsSupportTwoPos = clouds_at_index[h1.end].intersection(clouds_at_index[h2.start])
        seq = []  
        minMEC = 1000000
        bestBestHaps = []
        print h1.end, h2.start
        if h1.end + 1 == h2.start:
            for h in bestHaps:
                print h.seq
                tempSeq = h[h.start : h1.end] + [-1,-1] + h[h2.start+1 : h.end+1]
                print tempSeq
                newH = haplotype.Haplotype(h.start)
                newH.set_seq(tempSeq)
                MEC = newH.assign_clouds(clouds_at_index)
                print "new MEC", MEC
                if MEC < minMEC:
                    minMEC = MEC
                    bestBestHaps = []
                    bestBestHaps.append(newH)
                elif MEC == minMEC:
                    bestBestHaps.append(newH)
            if len(bestBestHaps) >= 2:
                print "same MEC try one fail"
                bestBestHaps = list()
                for c in cloudsSupportTwoPos:
                    print "support two pos reads", c.name, c.seq
                    ss = c[h1.end : h2.start+1]
                    print "ss", ss
                    seq.append(ss)
                if len(seq) >= 1:
                    ss = tools.most_frequency_list(seq)
                    print "most ss", ss
                    if len(ss) > 0:
                        if (ss[0] == h1[-1]  and ss[-1]== h2[0]) or (ss[0] != h1[-1]  and ss[-1] != h2[0]):
                            temp_seq = h1.seq + h2.seq
                        else:
                            temp_seq = h1.seq + tools.list_reverse(h2.seq) 
                        h = haplotype.Haplotype(h1.start)        
                        h.set_seq(temp_seq) # cloud not updata
                        h.assign_clouds(clouds_at_index)
                        print "same MEC try two success"
                        bestBestHaps.append(h)
                        return bestBestHaps        
                print "same MEC try two to be continue"
        return bestBestHaps  
        #sys.exit()


    #def updata_clouds_support_and_calculate_MEC(self, h1, clouds):
    def updata_clouds_support_and_calculate_MEC(self, h1, clouds_at_index):
        MEC = 0
        h1.clouds = set()
        #h2 = haplotype.Haplotype(h1.start)
        #h2.set_seq(tools.list_reverse(h1.seq))
        #print h1.seq, h2.seq
        h1.assign_clouds(clouds_at_index)
        return h1.MEC 
        
        
    def fill_gap_phasing_unsure_clouds(self, haps): # useless

        MEC = 0
        count = 0
        for i, h in enumerate(haps):
            print "before fill gap"
            h.printH()
            #h.fill_gap_inside()
            h.fill_gap_inside_outside()
            print "after fill gap and before deal unsure_clouds"
            h.printH()
            if len(h.unsure_clouds) !=0 :
                h.deal_unsure_clouds()
            MEC += h.calc_MEC()
            count +=1

            print "after deal unsure_clouds"
            h.printH()


    def pick_snp(self, clouds_at_index): # useless
        # some snp all have 0, first round ignore them
        print "total snp", self.length
        good_snp = list()
        for j in range(self.data_start, self.data_end+1):
            #print "snp no", j
            zero_label = False
            one_label = False
            for c in clouds_at_index[j]:
                if c[j] == 0:
                    zero_label = True
                elif c[j] == 1:
                    one_label = True
                if zero_label and one_label:
                    good_snp.append(j)
                    break
        print "good snp", len(good_snp)
        return good_snp


    def set_k_mers_only_for_good_snps(self, clouds_at_index, good_snps): # useless and if use this, have bug about consistency between k_mer and its' index
        lenSNP =  len(good_snps)
        for j in range(lenSNP-self.k):
            k_mer = {}
            cover_k_mer_clouds = set()

            #get all clouds cover j or j+1 or ... j + self.k
            for jj in range(j,j+self.k):
                for c in clouds_at_index[ good_snps[jj] ]:
                    cover_k_mer_clouds.add(c)

            for c in cover_k_mer_clouds:
                temp = []
                for jj in range(j,j+self.k):
                    temp.append( c[ good_snps[jj] ] )
                temp = tuple(temp)  # kmer is tuple type , easy translate with list j 
                if temp not in k_mer:
                    k_mer[temp] = []
                k_mer[temp].append(c.name)   
            #print "position: ", j, k_mer
            #self.k_mers[j] = tools.sorted_map_value_len(k_mer)
            self.k_mers[j] = k_mer
            self.k_mers_01[j] = []

            sorted_k_mer = tools.sorted_map_value_len(k_mer)
            for (key, v) in sorted_k_mer:
                #if key.find("-1") == -1: # for string
                if key.count(-1) == 0: # for tuple
                    self.k_mers_01[j].append(key)

    def set_k_mers(self, clouds_at_index):   # idea: skip some bad snp
        #print "last kmer pos",self.data_end-self.k
        #sys.exit()
        for j in range(self.data_start, self.data_end+2-self.k): # all SNP position
            k_mer = {}
            cover_k_mer_clouds = set()

            #get all clouds cover j or j+1 or ... j + self.k
            for jj in range(j,j+self.k):
                for c in clouds_at_index[jj]:
                    cover_k_mer_clouds.add(c)

            for c in cover_k_mer_clouds:
                #print "name c[j] c[j+1] c[j+2]", c.name, c[j], c[j+1], c[j+2]
                temp = tuple(c[j:j+self.k])                          # kmer is tuple type , easy translate with list j 
                #temp = ''.join( str(e) for e in c[ j : j+self.k ] ) # kmer is string type
                if temp not in k_mer:
                    k_mer[temp] = []
                #k_mer[temp].append(c.name)  
                k_mer[temp].append(c)   
            #print "position: ", j, k_mer
            #self.k_mers[j] = tools.sorted_map_value_len(k_mer)
            self.k_mers[j] = k_mer
            self.k_mers_01[j] = []

            sorted_k_mer = tools.sorted_map_value_len(k_mer)
            for (key, v) in sorted_k_mer:
                #if key.find("-1") == -1: # for string
                if key.count(-1) == 0: # for tuple
                    self.k_mers_01[j].append(key)


    def generate_haps_from_k_mers(self, start, end, label, clouds_at_index):  # label control level of bit-complement
        
        haps = list()
        temp_seq = list()
        for j in range(start, end - self.k + 2): # ?? 
            boolKmers = self.k_mers_01[j]
            if self.check_bit_complement(boolKmers, self.k_mers[j], label):
                if len(temp_seq) == 0:
                    h = haplotype.Haplotype(j) 
                    temp_seq.extend(boolKmers[0])
                elif temp_seq[-self.k+1:] == list(boolKmers[0])[:-1]: # may have reads conflict 
                    temp_seq.append(boolKmers[0][-1])
                elif tools.is_bool_reverse( temp_seq[-self.k+1:], boolKmers[0][:-1] ):
                    temp_seq.append( tools.int_reverse( boolKmers[0][-1] ) )
                else:
                    #print j, "kmer break type 2, not satisfy consistency"
                    assert len(temp_seq) != 0
                    if len(temp_seq) != 0:
                        h.set_seq(temp_seq)
                        haps.append(h)
                        temp_seq = list()
            else:
                if len(temp_seq) != 0:
                    h.set_seq(temp_seq)
                    haps.append(h)
                    temp_seq = list()
                #if len(boolKmers) == 0:
                #print j, "kmer break type 3, no reads cover all position in k_mer"
                #else:
                #print j, "kmer break type 1, not satisfy bit-complement"
        if len(temp_seq) != 0:
            h.set_seq(temp_seq)
            haps.append(h)

        #self.unsure_clouds = set()
        for h in haps:
            h.assign_clouds(clouds_at_index)  
            #if h.check_clouds_intersection() == False:
                #h.unsure_clouds.update( h.remove_intersection() ) 
        return haps  
    
    ############################################### case study
    #some case left, first and third bit-complement
    #second = first
    ############################################
    # very low coverage like 3
    


    #07/Jan./2019
    def check_bit_complement2(self, boolKmers, kmers, label): # after talking code with yu, maybe more reasonable 
        lenKmer = len(boolKmers) # bool_k_mers a list of kmer which not include -1
        curCov = 0
        #for k in kmers:
        for k in boolKmers:
            curCov += len(kmers[k])
        if lenKmer >= 2:
            kmer0cov = len( kmers[boolKmers[0]] )
            kmer1cov = len( kmers[boolKmers[1]] )
            if lenKmer >= 3:
                kmer2cov = len( kmers[boolKmers[2]] )

        if label == 1:
            if lenKmer == 1:
                return True
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1]) and kmer0cov + kmer1cov > 0.5*curCov ):
                return True
            if lenKmer >= 3 and kmer1cov == kmer2cov :
                if ( tools.is_bool_reverse( boolKmers[0], boolKmers[2]) and kmer0cov + kmer2cov > 0.5*curCov ):
                    return True
                if kmer0cov == kmer1cov:
                    if ( tools.is_bool_reverse( boolKmers[1], boolKmers[2]) and kmer1cov + kmer2cov > 0.5*curCov ):
                        return True


        if label == 2:        
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1]) and kmer0cov + kmer1cov > 0.66*curCov ):
                return True
            if lenKmer >= 3 and kmer1cov == kmer2cov : # 3 2 2
                if ( tools.is_bool_reverse( boolKmers[0], boolKmers[2]) and kmer0cov + kmer2cov > 0.66*curCov ):
                    return True
                if kmer0cov == kmer1cov: # 3 3 3
                    if ( tools.is_bool_reverse( boolKmers[1], boolKmers[2]) and kmer1cov + kmer2cov > 0.66*curCov ):
                        return True
        return False    


    def check_bit_complement(self, boolKmers, kmers, label): 

        # most loose 
        lenKmer = len(boolKmers) # bool_k_mers a list of kmer which not include -1
        if label == 1:
            if ( ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1] ) ) or 
                    lenKmer == 1 or 
                    ( lenKmer==3 and  (tools.is_bool_reverse(boolKmers[0], boolKmers[1]) or 
                                         tools.is_bool_reverse(boolKmers[0], boolKmers[2])) ) ):
                return True
            return False

        if label == 2:
            #a bit higher coverage, need improve 
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1])  ):
                return True
            return False

        if label == 3:
            curCov = 0
            for k in boolKmers:
                #print "kmer detail", k, len(self.k_mers[j][k])
                curCov += len(kmers[k])
            if lenKmer >= 2:
                kmer1cov = len( kmers[boolKmers[0]] )
                kmer2cov = len( kmers[boolKmers[1]] )
                
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1]) and kmer1cov + kmer2cov > 0.5*curCov ):
                    #print kmer1cov
                    #print kmer2cov
                    #print "cur cov", cur_cov
                    # for pacbio and nanopore, should I use different para
                    return True
            return False    
            

    def fill_gap_between_haps(self, haps):  # useless
        
        l = len(haps)
        print "before check, number of haplotype", l
        i = 0
        while True:
            while ( l>=2 and self.fill_two_haps(haps[i], haps[i+1]) ):
                haps.remove(haps[i+1])
                l = len(haps) - i
            i += 1
            if i+1 >= len(haps):
                break
        print "after fill gap bwtween haps, number of haplotype", len(haps)
    
           
    def fill_two_haps(self, h1, h2): #useless
        assert h1.start <= h2.start and h1.end <= h2.end
        #assert h1.end < h2.start
        print "xx", h1.start, h1.end
        print "xx", h2.start, h2.end
        h1.fill_gap_next(h2.end)
        h2.fill_gap_pre(h1.start)
        print "xx fill two haps"
        if h1.end < h2.start:
            return False # no overlap 
        if h1.end > h2.start:
            print h1.start, h1.end, h2.start, h2.end
            assert h1.end <= h2.end
            s1 = h1[ h2.start : h1.end+1 ]
            s2 = h2[ h2.start : h1.end+1 ]
            print s1
            print s2
            #h1.printH()
            #h2.printH()
            lenS1 = tools.count(s1)
            lenS2 = tools.count(s2)
            lenS = min(lenS1, lenS2)
            print lenS1, lenS2
            similarRate = float( tools.similar_distance(s1, s2) ) / lenS
            diffRate = float( tools.hamming_distance(s1, s2) ) / lenS
            print "similarRate diffRate ", similarRate, diffRate 
            print "after check two"
            if (similarRate > 0.7 and lenS >= 5) or (similarRate > 0.9 and lenS >= 3 and diffRate < 0.01) : # need improve
                temp_seq = h1.seq + h2[ h1.end+1: h2.end + 1 ]
                h1.set_seq( temp_seq ) 
                h1.left_clouds.update(h2.left_clouds)
                h1.right_clouds.update(h2.right_clouds)
                print "fill and merge"
                h1.printH()
                return True
            if (diffRate > 0.7 and lenS >= 5) or (diffRate > 0.9 and lenS >= 3 and similarRate < 0.01) : # need improve
                temp_seq = h1.seq + tools.list_reverse( h2[ h1.end+1 : h2.end+1] )
                h1.set_seq( temp_seq )
                h1.left_clouds.update(h2.right_clouds)
                h1.right_clouds.update(h2.left_clouds)
                print "fill and merge"
                h1.printH()
                return True
            halfOverlap = len(s1)/2
            #print halfOverlap
            #print h1.seq
            temp_seq = h1.seq [ : - halfOverlap]
            while len(temp_seq) > 1 and temp_seq[-1] == -1:
                temp_seq = temp_seq[:-1]
            #print temp_seq
            h1.set_seq(temp_seq)
            
            temp_seq = h2.seq [ halfOverlap : ]
            h2.start = h2.start + halfOverlap
            while len(temp_seq) > 1 and temp_seq[0] == -1:
                temp_seq = temp_seq[1:]
                h2.start = h2.start + 1

            h2.set_seq(temp_seq)

            print "xx", h1.start, h1.end
            print "xx", h2.start, h2.end
            h1.printH()
            h2.printH()
            return False    


    def merge_haps(self, haps): # useless
        
        l = len(haps)
        print "before merge, number of haplotype", l
        i = 0
        while True:
            while ( l>=2 and self.from_two_haps(haps[i], haps[i+1]) ):
                haps.remove(haps[i+1])
                l = len(haps) - i
            i += 1
            if i+1 >= len(haps):
                break
        print "after merge, number of haplotype", len(haps)
        
    def from_two_haps(self, h1, h2, clouds_at_index):  # part2 change
        assert h1.start < h2.start and h1.end < h2.end
        
        if h1.check_clouds_intersection() == False:
            print "h1 imposible"
            for c in h1.left_clouds:
                print c, 
            print "\n"    
            for c in h1.right_clouds:
                print c,


        assert h1.check_clouds_intersection() # should not have
        assert h2.check_clouds_intersection() # intersection
       

        llR = h1.left_clouds.intersection(h2.left_clouds)
        lrR = h1.left_clouds.intersection(h2.right_clouds)
        rlR = h1.right_clouds.intersection(h2.left_clouds)
        rrR = h1.right_clouds.intersection(h2.right_clouds)

        print len(llR), len(lrR), len(rrR), len(rlR)
        llR2 = set()
        lrR2 = set()
        rrR2 = set()
        rlR2 = set()
        for c in llR:
            if len(c.seq) != 2:
                llR2.add(c)

        for c in lrR:
            if len(c.seq) != 2:
                lrR2.add(c)

        for c in rlR:
            if len(c.seq) != 2:
                rlR2.add(c)
            
        for c in rrR:
            if len(c.seq) != 2:
                rrR2.add(c)
        ll = len(llR2)
        lr = len(lrR2) 
        rl = len(rlR2) 
        rr = len(rrR2) 

        print ll, lr, rr, rl
        hapsFromConnectTwoHaps = list() 
        #check reads intersection and overlap
        #if (ll > lr and rr > rl) or ( (ll>0 or rr>0) and (lr==0 or rl==0) ):
        # need improve
        #if (ll > lr and rr > rl) or  (ll == lr and rr > rl) or (ll > lr and rr == rl) :
        if (ll+rr > lr+rl): # according to MEC
            #print "merge two ", h1.start, h1.end, h2.start, h2.end
            #print h1.seq,h2.seq
            #print ll, lr, rr, rl
            if h1[ h2.start : h1.end+1 ] !=  h2[ h2.start : h1.end+1 ] : # need update, no need exactly =
                print 'reads intersection conflict with overlap info, need imporve, overlap no need exactly correct '
                hapsFromConnectTwoHaps.append(h1)
                hapsFromConnectTwoHaps.append(h2)
                return hapsFromConnectTwoHaps
            temp_seq = h1.seq + h2[ h1.end+1: h2.end + 1 ]
            h1.set_seq( temp_seq )
            print "after merge two ", h1.start, h1.end
            print h1.seq
            h1.clouds =set()
            h1.assign_clouds(clouds_at_index)
            hapsFromConnectTwoHaps.append(h1)
            return hapsFromConnectTwoHaps

        #elif (lr > ll and rl > rr) or ( lr == ll and rl > rr) or ( lr > ll and rl == rr ):
        elif (lr+rl > ll+rr):
            print "merge two ", h1.start, h1.end, h2.start, h2.end
            print h1.seq, h2.seq
            print ll, lr, rr, rl
            ##################################### 
            #   72      73     74  75 76
            #   h1.s          h1.e 
            #           h2.s          h2.e
            ####################################
            print "check", h1[ h2.start : h1.end+1 ] , h2[ h2.start : h1.end+1 ]
            print tools.list_reverse( h2[ h2.start : h1.end+1 ]) 
            if h1[ h2.start : h1.end+1 ] != tools.list_reverse( h2[ h2.start : h1.end+1 ] ):
                print 'reads intersection conflict with overlap, need imporve, overlap no need exactly correct'
                hapsFromConnectTwoHaps.append(h1)
                hapsFromConnectTwoHaps.append(h2)
                return hapsFromConnectTwoHaps
            temp_seq = h1.seq + tools.list_reverse( h2[ h1.end+1 : h2.end+1] )
            h1.set_seq( temp_seq )
            print "after merge two ", h1.start, h1.end 
            print h1.seq
            h1.clouds =set()
            h1.assign_clouds(clouds_at_index)
            hapsFromConnectTwoHaps.append(h1)
            return hapsFromConnectTwoHaps
        else:
            print "have overlap, can not merge haps ll, lr, rr, rl", ll, lr, rr, rl

            
            print "ll" 
            for c in llR:
                print c.name, c.seq

            print "lr" 
            for c in lrR:
                print c.name, c.seq

            print "rr" 
            for c in rrR:
                print c.name, c.seq

            print "rl" 
            for c in rlR:
                print c.name, c.seq
                
            temp_seq = h1[h1.start: h2.start]
            h1.set_seq( temp_seq ) # have overlap, cann't connect, remove the tail of h1
            h1.clouds =set()
            h1.assign_clouds(clouds_at_index) 
            hapsFromConnectTwoHaps.append(h1)
            hapsFromConnectTwoHaps.append(h2)
            print "reads intersection provide no info to merge, but overlap may helpful, need imporve "
            return hapsFromConnectTwoHaps
                
