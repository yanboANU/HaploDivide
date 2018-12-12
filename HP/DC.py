#########################################################################
# File Name: DC.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Wed 24 Oct 2018 15:38:51 AEDT
#########################################################################
#!/bin/bash

import sys
import tools
from libprism.local.prepare import verify_complexity, blocks_from_clouds, contigs_from_clouds
import cloud
import copy
import numpy as np
import operator
class Haplotype: # struct
    def __init__(self, start):
        self.start = start
        self.seq = list()
        self.len = 0
        self.end = 0
        self.left_clouds = set() # element is clouds name => clouds
        self.right_clouds = set()
        self.unsure_clouds = set()


    def __len__(self):
        return len(self.seq)

    def __getitem__(self,index): # directly []
        start, end = self.start, self.end
        assert self.len == len(self.seq)
        if index < start or index > end:
            return -1
        else:
            return self.seq[index-start]
       
       
    def __setitem__(self,index, val):

        self.seq[index-self.start] = val

    def __getslice__(self,s,e): # directly [s,e), yanbo
        ans = []
        index = s
        while index < e:
            ans.append(self.__getitem__(index))
            index += 1
        return ans        

     
    def find_gaps(self):
        gaps = {}
        mid_gap_pos = []
        for i in range(self.start, self.end+1):
            if self[i] == -1:
                mid_gap_pos.append(i)
        #print "mid gap pos", mid_gap_pos

        for c in self.left_clouds:
            # before
            for i in range(c.start, self.start):
                if i not in gaps and c[i] != -1:
                    gaps[i] = []
                if c[i] != -1:    
                    gaps[i].append(c[i])        
            # after
            for i in range(self.end+1, c.end+1):
                if i not in gaps and c[i] != -1:
                    gaps[i] = []
                if c[i] != -1:    
                    gaps[i].append(c[i])        

            for i in mid_gap_pos:
                if c[i] != -1:
                    if i not in gaps:
                        gaps[i] = []
                    gaps[i].append(c[i])
                

        for c in self.right_clouds:
            for i in range(c.start, self.start):
                if i not in gaps and c[i] != -1:
                    gaps[i] = []
                if c[i] != -1:    
                    gaps[i].append(tools.int_reverse(c[i]))        
            for i in range(self.end+1, c.end+1):
                if i not in gaps and c[i] != -1:
                    gaps[i] = []
                if c[i] != -1:
                    gaps[i].append(tools.int_reverse(c[i]))      
            for i in mid_gap_pos:
                if c[i] != -1:
                    if i not in gaps:
                        gaps[i] = []
                    gaps[i].append(tools.int_reverse(c[i]))
        return gaps
    

    def fill_gap_next(self, stop):
        print "aa"
        gaps= self.find_gaps()
        s = self.start
        e = self.end
        sortedGapsKey = sorted (gaps.keys())
        l = len(gaps)
        if l == 0:
            return
        print sortedGapsKey 
        if sortedGapsKey[l-1] > e: # find new end
            i = l-1
            key = e
            while i >= 0 and sortedGapsKey[i] > stop:
                i -= 1
            while i >= 0 and sortedGapsKey[i] > e:
                if tools.most_frequency( gaps [ sortedGapsKey[i] ] ) != -1:
                    key = sortedGapsKey[i]
                    break
                i -= 1
            if key > e:     
                add = [ tools.most_frequency(gaps[key]) ]
                temp_seq = self.seq + (key-self.end-1) * [-1] + add
                self.set_seq(temp_seq)
        for key in gaps:
            if key <= self.end and key >= self.start:
                self[key] = tools.most_frequency(gaps[key])
        print "new start, new end", self.start, self.end    
        print self[self.start]
        print self[self.end]
        assert self[self.start] != -1 and self[self.end] != -1

    def fill_gap_pre(self, stop):
        gaps= self.find_gaps()
        s = self.start
        e = self.end
        sortedGapsKey = sorted (gaps.keys())
        l = len(gaps)
        if l == 0:
            return
        if sortedGapsKey[0] < s:  # find new start
            i = 0
            key = s

            while i < l and sortedGapsKey[i] < stop:
                i += 1     
            while i < l and sortedGapsKey[i] < s:
                if tools.most_frequency( gaps [ sortedGapsKey[i] ]  ) != -1:
                    key = sortedGapsKey[i]
                    break
                i += 1
            if key < s:     
                add = [ tools.most_frequency(gaps[key]) ]
                temp_seq = add + (self.start-key-1) * [-1] + self.seq
                self.start = key 
                self.set_seq(temp_seq)

          
        for key in gaps:
            if key <= self.end and key >= self.start:
                self[key] = tools.most_frequency(gaps[key])
        print "new start, new end", self.start, self.end    
        print self[self.start]
        print self[self.end]
        assert self[self.start] != -1 and self[self.end] != -1


    def fill_gap_inside_outside(self):
        gaps= self.find_gaps()
        print "gaps", gaps
        #h.fill_seq(gaps)
        s = self.start
        e = self.end
        sortedGapsKey = sorted (gaps.keys())
        l = len(gaps)
        if l == 0:
            return
        print self.start, self.end          
        print sortedGapsKey
        # find smallest and biggest postion of gaps
        # samllest and biggest should coles to orignal start and end
        if sortedGapsKey[0] < s:  # find new start
            i = 0
            key = s
            #while i < l and sortedGapsKey[i] < s-10:
                #i += 1     
            while i < l and sortedGapsKey[i] < s:
                if tools.most_frequency( gaps [ sortedGapsKey[i] ]  ) != -1:
                    key = sortedGapsKey[i]
                    break
                i += 1
            if key < s:     
                add = [ tools.most_frequency(gaps[key]) ]
                temp_seq = add + (self.start-key-1) * [-1] + self.seq
                self.start = key 
                self.set_seq(temp_seq)

          
        if sortedGapsKey[l-1] > e: # find new end
            i = l-1
            key = e
            #while i >= 0 and sortedGapsKey[i] > e+10:
                #i -= 1
            while i >= 0 and sortedGapsKey[i] > e:
                if tools.most_frequency( gaps [ sortedGapsKey[i] ] ) != -1:
                    key = sortedGapsKey[i]
                    break
                i -= 1
            if key > e:     
                add = [ tools.most_frequency(gaps[key]) ]
                temp_seq = self.seq + (key-self.end-1) * [-1] + add
                self.set_seq(temp_seq)

        #sys.exit()
        for key in gaps:
            if key <= self.end and key >= self.start:
                self[key] = tools.most_frequency(gaps[key])
        print "new start, new end", self.start, self.end    
        print self[self.start]
        print self[self.end]
        assert self[self.start] != -1 and self[self.end] != -1

        
        '''
        
        #print "cc", key, self.start, self.end
        # fill mid and also extend two sides 
        if key < s-10 or key > e+10: # should not fill gap too much
            continue
        if key < self.start:
            #print "dd", key, self.start, self.end
            add = [ tools.most_frequency(gaps[key]) ]
            if add[0] == -1: # add only has one element
                continue
            temp_seq = add + (self.start-key-1) * [-1] + self.seq
            self.start = key
            self.set_seq(temp_seq)
        elif key > self.end:
            add = [ tools.most_frequency(gaps[key]) ]
            if add[0] == -1:
                continue
            temp_seq = self.seq + (key-self.end-1) * [-1] + add
            #self.end = key # no need, every time set_seq, will calculate end every time
            self.set_seq(temp_seq)
        else:
            assert key <= self.end and key >= self.start
            #print "bb", key, self.start, self.end, gaps[key]
            self[key] = tools.most_frequency(gaps[key])
            #print self[key]
        
        # only fill mid gap
        if key <= self.end and key >= self.start:
            self[key] = tools.most_frequency(gaps[key])
        '''    
    def fill_gap_inside(self):
        gaps= self.find_gaps()
        #print "gaps", gaps
        #h.fill_seq(gaps)
        for key in gaps:
            # only fill mid gap
            if key <= self.end and key >= self.start:
                self[key] = tools.most_frequency(gaps[key])





    def calc_MEC(self):

        #print self.start
        #print self.seq
        ans = 0
        for c in self.left_clouds:
            s1 = c.seq
            s2 = self[c.start : c.end+1]
            #print s1, s2
            ans += tools.hamming_distance(s1,s2)  
        print "left MEC", ans    
        for c in self.right_clouds:
            s1 = tools.list_reverse(c.seq)
            s2 = self[c.start : c.end+1]
            #print s1, s2
            ans += tools.hamming_distance(s1,s2)  
        print "right MEC", ans 
        return ans
   
    def deal_unsure_clouds(self):
        print "deal unsure clouds" 

        sure = list() 
        for c in self.unsure_clouds:
            s0 = c.seq
            s1 = tools.list_reverse(s0)
            s2 = self[c.start : c.end+1] # haplotype seq 
            print s0, s1, s2
            if tools.hamming_distance(s0, s2) < tools.hamming_distance(s1, s2):
                self.left_clouds.add(c)    
                sure.append(c)
            elif tools.hamming_distance(s0, s2) > tools.hamming_distance(s1, s2): 

                self.right_clouds.add(c)    
                sure.append(c) 
        #print sure
        #sys.exit()
        for c in sure:
            self.unsure_clouds.remove(c)

    def printH(self):
        print "start", self.start, "end", self.end
        print "seq", self.seq
        print "left_clouds", len(self.left_clouds)
        print "right_clouds", len(self.right_clouds)
        print "unsure_clouds", len(self.unsure_clouds)
        print "\n"


    def update_clouds(self, left_kmer, kmers):
        right_kmer = tools.bool_reverse(left_kmer)
        #print left_kmer, right_kmer
        #print kmers
        for key in kmers:
            if tools.hamming_distance(key, left_kmer) < tools.hamming_distance(key, right_kmer):
                self.left_clouds.update(kmers[key])
            elif tools.hamming_distance(key, left_kmer) > tools.hamming_distance(key, right_kmer):
                self.right_clouds.update(kmers[key])
        #print self.left_clouds
        #print self.right_clouds

    def set_seq(self, seq):
        self.seq = seq
        self.len = len(seq)
        self.end = self.start + self.len - 1

    # left clouds should not have intersection with right
    def check_clouds_intersection(self):    
        if len( self.left_clouds.intersection(self.right_clouds) ) == 0:
            return True
        return False

    def remove_intersection(self):
        #print "cc"
        inter_clouds = self.left_clouds.intersection(self.right_clouds)
        #print "inter clouds", inter_clouds
        #print "len inter clouds", len(inter_clouds)
        if len(inter_clouds)>= 1:
            #print "dd inter clouds", inter_clouds
            self.left_clouds = self.left_clouds-inter_clouds
            self.right_clouds = self.right_clouds-inter_clouds
        return inter_clouds    
    


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
        #fout.write("%s\n" % ''.join( map(str, c[h1.start:h2.end+1]) ))
        #sys.exit()
        

class DC(object):

    def generate_blocks(self, clouds_at_index):

        blocks = list() 
        indexLen = len(clouds_at_index)
        pos = self.find_position_no_reads_span(0, indexLen, clouds_at_index)
        l = len(pos)
        for i in range(l-1):
            j=i+1
            if pos[j] -pos[i] > 1:
                blocks.append( (pos[i]+1, pos[j] ) )
        if indexLen > pos[-1]:
            blocks.append( (pos[-1], indexLen) )
        return blocks        

    def __init__(self, k, data_start, data_end, all_clouds, clouds_at_index): # all_clouds, a cloud list
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
        # contigs_at_index = contigs_from_clouds(all_clouds, self.data_start, self.data_end)

        # optional   
        '''
        good_snps = self.pick_snp(clouds_at_index)
        self.set_k_mers_only_for_good_snps(clouds_at_index, good_snps)
        '''
        
        self.set_k_mers(clouds_at_index)
        #self.blocks = blocks_from_clouds(all_clouds, clouds_at_index) # probhap generate block
        self.blocks = self.generate_blocks(clouds_at_index)

        ''' 
        uncovered_positions = set()
        # detect positions that have no coverage:
        for j in range(self.data_start, self.data_end+1):
            if len(clouds_at_index[j]) == 0:
                uncovered_positions.add(j)

        self.uncovered_positions = uncovered_positions
        '''

    def run_dc_v1(self, clouds_at_index): # first version

        haps = self.generate_haps_from_k_mers() # slide window, connect by overlap statify (1)bit-complement
        #pre_start = -1
        #pre_end = -1
        self.merge_haps(haps) # using reads support two haps connect haplotype and overlap
        self.fill_gap_phasing_unsure_clouds(haps) # clouds is reads
                                                  # two step of fill gap
                                                  # fill gap inside of hap: finish
                                                  # fill gap between hap: not yet
        ######################################################################                                          
        #self.connect_by_unsure_clouds(haps) # using first round unsure clouds
        return haps

    def enumerate_seq(self, start, end, clouds_at_index):

        print "enumerate block, length:", end - start + 1 
        minMEC = 10000
        best_hap = []
        dis = end - start + 1
        allPosSeq = tools.enumerate_01_list(dis)    
        l = len(allPosSeq)/2
        clouds = cloud.get_clouds(start, end, clouds_at_index)
        
        for seq in allPosSeq[0:l]:
            h = Haplotype(start)        
            h.set_seq(seq)
            MEC = self.calculate_MEC(h, clouds)
            if MEC < minMEC:
                minMEC = MEC
                best_hap = []
                best_hap.append(h)
            elif MEC == minMEC:
                best_hap.append(h)
        #print "min MEC", minMEC
        #print "best hap", best_seq
        l = len(best_hap)
        if l == 1:
            #print "enumerate block susscess"
            return best_hap
        else:
            #print "enumerate block fail"
            seqs = [] 
            for h in best_hap:
                #print h.seq
                seqs.append(h.seq)
            arr = np.sum(seqs, axis=0)    
            temp_seq = []
            for ele in arr:
                if ele == 0:
                    temp_seq.append(0)
                elif ele == l:
                    temp_seq.append(1)
                else:
                    temp_seq.append(-1)
            #print temp_seq
            h = Haplotype(start)        
            h.set_seq(temp_seq)
            return [h]
            

    def run_dc_in_block(self, start, end, clouds_at_index):

        haps = list() 
        haps = self.generate_haps_from_k_mers(start, end)

        #return haps  
         
        print "after generate haps from kmers ", len(haps)
        if len(haps) == 0:
            haps = self.enumerate_seq(start, end, clouds_at_index)
            if len(haps) > 0:    
                print "block no haps, so enumerate:", haps[0].start, haps[0].end, haps[0].seq
            return haps    
        elif len(haps) == 1:
            print haps[0].start, haps[0].end, haps[0].seq
            #haps[0].fill_gap_pre(start)
            #haps[0].fill_gap_next(end)
            return haps
        else:
            #print_aligned_region(start, end, clouds_at_index)
            self.enumerate_gap_between_haps(haps, clouds_at_index)
            #print "after enumerate gap ", len(haps)
            #self.fill_gap_between_haps(haps)
            #print "after fill gaps ", len(haps)
            #haps[0].fill_gap_pre(start)
            #haps[-1].fill_gap_next(end)
            return haps
        

    def run_dc(self, clouds_at_index): # 12, Dec

        #print self.find_position_no_reads_span(0, len(clouds_at_index), clouds_at_index)
        #for block in self.blocks:
            #print block
        #sys.exit() 

        allHaps = []
        count = 1
        print "There are", len(self.blocks), "blocks"
        for block in self.blocks:  # blo/cks break at no reads cover region
            print "phasing block: ", count, ", region ", block[0] , block[-1]
            allHaps.extend(self.run_dc_in_block(block[0], block[-1], clouds_at_index))
            count += 1
 
        self.enumerate_gap_between_haps(allHaps, clouds_at_index) # only need work for between block, need to updata
        allHaps[0].fill_gap_pre(-1)
        allHaps[-1].fill_gap_next(300000000)
        return allHaps  


    def run_dc_v3(self, clouds_at_index): # 6, Dec

        haps = self.generate_haps_from_k_mers() # slide window, connect by overlap statify (1)bit-complement

        # write
        
        fout1 = open("step1.out" , "w")
        for i, block in enumerate(haps):
            fout1.write("BLOCK: offset: %d len: %d positions: %d\n" % (block.start, block.len, -1))
            for k,y in enumerate(block.seq):
                j = block.start + k
                if y != -1:
                    fout1.write("%d\t%d\t%d\t\n" % (j, y, tools.int_reverse(y)))
            fout1.write("********\n")
        fout1.close()
        



        '''
        self.add_haps(haps)
        fout2 = open("add.out" , "w")
        for i, block in enumerate(haps):
            fout2.write("BLOCK: offset: %d len: %d positions: %d\n" % (block.start, block.len, -1))
            for k,y in enumerate(block.seq):
                j = block.start + k
                if y != -1:
                    fout2.write("%d\t%d\t%d\t\n" % (j, y, tools.int_reverse(y)))
            fout2.write("********\n")
        fout2.close()
        #sys.exit()
        '''  

        #enumerate
        l = len(haps)
        print "before enumerate, number of haplotype", l
        self.enumerate_gap_between_haps(haps, clouds_at_index)

        
        fout2 = open("enumerate.out" , "w")
        for i, block in enumerate(haps):
            fout2.write("BLOCK: offset: %d len: %d positions: %d\n" % (block.start, block.len, -1))
            for k,y in enumerate(block.seq):
                j = block.start + k
                if y != -1:
                    fout2.write("%d\t%d\t%d\t\n" % (j, y, tools.int_reverse(y)))
            fout2.write("********\n")
        fout2.close()
        
        
        self.fill_gap_between_haps(haps)
        haps[0].fill_gap_pre(-1)
        print "fill last haps "
        haps[-1].fill_gap_next(300000000)
        for h in haps:
            h.fill_gap_inside()

        l = len(haps)
        print "after fill gap, number of haplotype", l

        return haps

    def enumerate_gap_between_haps(self , haps, clouds_at_index):

        MEC = 0
        i = 0
        l = len(haps)
        fout = open("enumerate_fail" ,"w")
        while True:
            while ( l>=2 and self.enumerate_two_haps(haps[i], haps[i+1], clouds_at_index, fout) ): # gap >=0, then enumarate 
                                                                                             # gap <0. then have overlap, merge    
                haps.remove(haps[i+1])
                print "enumerate two haps and merge happen", haps[i].start, haps[i].end
                l = len(haps) - i
            i += 1
            if i+1 >= len(haps):
                break
        print "after enumerate gap bwtween haps, number of haplotype", len(haps)

    
    def add_haps(self , haps):

        print "add haps" 
        l = len(haps)
        add_haps = list()
        i = 0
        #while True:
        while (l>=2):
            dis =  haps[i+1].start - haps[i].end - 1 
            if dis >= 10 and dis < 20:
                h = Haplotype(haps[i].end + dis/2)
                h.set_seq([1])
                haps.insert(i+1, h)
                print "add", i+1, h.start
            elif dis >= 20:
                h = Haplotype(haps[i].end + 10)
                h.set_seq([1])
                haps.insert(i+1, h)
                print "add", i+1, h.start
            l = len(haps) - i
            i += 1
            if i+1 >= len(haps):
                break

    '''  
    def find_position_no_reads_span(self, h1, h2, clouds_at_index): # between two haplotype
        assert h2.start > h1.end
        pos = []

        for i in range(h1.end, h2.start):
            j=i+1

            #for c in clouds_at_index[i]:
                #print c[i]
            #sys.exit()    
            ij = clouds_at_index[i].intersection(clouds_at_index[j])
            print i,j, "were spaned by ", len(ij), "reads"
            if len(ij) == 0:
                pos.append(i)
        return pos        
    '''    
    
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

       
    def enumerate_two_haps(self, h1, h2, clouds_at_index, fout):

        minMEC = 10000
        best_hap = []
        print "enumerate gap between two haps, dis:", h2.start - h1.end -1   
        print h1.start, h1.end, h2.start, h2.end
        
        if h1.end >= h2.start:
            print "two segments have overlap"
            return self.from_two_haps(h1, h2) # satisfy reads intersection and overlap
            
        dis = h2.start - h1.end - 1
        clouds = cloud.get_clouds(h1.start, h2.end, clouds_at_index)
        #aligned_detail.append(c[h1.start : h2.end+1])
        
        if dis >= 21: # need improve
            print "dis is", dis, "too long to enumerate"
            return False
        
        print h1.start, h2.end, "this region supported reads number", len(clouds)    
        allPosSeq = tools.enumerate_01_list(dis)    

        for seq in allPosSeq:
            for i in range(0,2): # debug maybe h1.seq + seq + tools.reverse(h2.seq)
                if i == 0:
                    temp_seq = h1.seq + seq + h2.seq 
                elif i == 1:
                    temp_seq = h1.seq + seq + tools.list_reverse(h2.seq) # debug maybe h1.seq + seq + tools.reverse(h2.seq)

                h3 = Haplotype(h1.start)
                h3.set_seq(temp_seq)
                MEC = self.calculate_MEC(h3, clouds)
                if MEC < minMEC:
                    minMEC = MEC
                    best_hap = []
                    best_hap.append(h3)
                elif MEC == minMEC:
                    best_hap.append(h3)
        #print "min MEC", minMEC
        #print "best seq", best_seq

        l = len(best_hap)
        if l == 1:
            print "only one minimal MEC, fill and merge"
            #h1 = copy.deepcopy(best_hap[0])  # ??
            h1.set_seq(best_hap[0].seq)
            h1.left_clouds = best_hap[0].left_clouds
            h1.right_clouds = best_hap[0].right_clouds
            h1.unsure_clouds = best_hap[0].unsure_clouds

            return True
        elif l == 2:
            for h in best_hap:
                print h.seq


            idx = map(operator.eq, best_hap[0].seq, best_hap[1].seq )
            ###########
            #000100
            #000000
            #TTTFTT
            ###########
            unsure_pos = idx.index(False)
            idxF = idx.count(False)
            if idxF == 1 or idxF < len(idx)*0.1:
                temp_seq = best_hap[0].seq[ : unsure_pos]
                temp_seq = temp_seq + [-1]
                temp_seq = temp_seq + best_hap[0].seq[unsure_pos+1 : ]
                h1.set_seq(temp_seq)
                print "2 enumerate", h1.seq
                return True

            # two situations
            ###########
            #000100   000100-10   
            #000011   000011-11
            #TTTFFF   TTTFFF TF
            ###########
            left_idx = idx[:unsure_pos]
            right_idx = idx[unsure_pos:]
            left_temp = best_hap[0].seq[ : unsure_pos]
            right_temp =  best_hap[0].seq[unsure_pos : ]
            if left_idx.count(True) == len(left_idx) and (right_idx.count(False) == len(right_idx) or
                                                        right_idx.count(False) + right_temp.count(-1)  == len(right_idx)   ) :
                h1.set_seq(left_temp)
                h2.start = h1.start + unsure_pos
                h2.set_seq(right_temp)
                print "3 enumerate", h1.seq, h1.start
                print "3 enumerate", h2.seq, h2.start
                return False

           
            print "enumerate between gaps, unconsider case:"
            for h in best_hap:
                print h.seq
            print h1.start, h1.end, h2.start, h2.end 
            print h1.seq
            print h2.seq        
            print_aligned_region(h1.start, h2.end, clouds_at_index)
            #print temp_seq        
            
            return False
        else:
            print "enumerate between gaps, more than two MEC"
            return False
            

     
        return False



    def calculate_MEC(self, h1, clouds):
        MEC = 0
        h2 = Haplotype(h1.start)
        h2.set_seq(tools.list_reverse(h1.seq))
        #print h1.seq, h2.seq 
        for c in clouds:
            s1 = c.seq
            s2 = h1[c.start : c.end + 1]
            s3 = h2[c.start : c.end + 1]
            # this read has same distance to two haplotypes, should use this read to calculate MEC??
            d12 = tools.hamming_distance(s1, s2)
            d13 = tools.hamming_distance(s1, s3)
            if d12 < d13:
                h1.left_clouds.add(c)
                #MEC += d12
            elif d12 > d13:
                h1.right_clouds.add(c)
                #MEC += d13
            else:
                h1.unsure_clouds.add(c)
            MEC += min(d12, d13)  # compare with do not considering unsure reads
                                  # this one better                
        #print "MEC:", MEC    
        return MEC 
        

 
            
        
        


    def run_dc_v2(self, clouds_at_index): # second version: 5, Dec, new version pipeline of dc

        haps = self.generate_haps_from_k_mers() # slide window, connect by overlap statify (1)bit-complement
       
        fout1 = open("step1", "w")
        for h in haps:
             fout1.write("%s %s %s\n"   % (h.start, h.end, h.seq))
        
        #"before merge no inside gap"    
        self.merge_haps(haps) # using reads support two haps connect haplotype and overlap

        fout2 = open("step2", "w") 
        for h in haps:
            h.fill_gap_inside()
            fout2.write("%s %s %s\n"   % (h.start, h.end, h.seq))
            

        # check haps is sorted or not
        i = 0
        l = len(haps)
        while i+1 < l:
            h1 = haps[i]
            h2 = haps[i+1]
            assert h1.start < h2.start
            i+=1
 
        self.fill_gap_between_haps(haps)
        haps[0].fill_gap_pre(-1)
        print "fill last haps "
        #print len(haps[-1].left_clouds)
        #print len(haps[-1].right_clouds)
        #print len(haps[-1].unsure_clouds)
        haps[-1].fill_gap_next(300000000)
        for h in haps:
            h.fill_gap_inside()
            #print "check haps again is sorted or not", h.start, h.end    
        return haps

       

        
    def fill_gap_phasing_unsure_clouds(self, haps):

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

    



    def pick_snp(self, clouds_at_index): 
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



    def connect_by_unsure_clouds(self, haps):  # phasing unsure clouds and try connect haplotype
        print "connect by unsure clouds"
        sum = 0
        for h in haps:
            sum += len(h.unsure_clouds)
        print sum     

    def set_k_mers_only_for_good_snps(self, clouds_at_index, good_snps):
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
        for j in range(self.data_start, self.data_end+1-self.k): # all SNP position
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


     



    def generate_haps_from_k_mers(self, start, end):  # need to updata, different coverage, different fix_rule
        
        haps = list()
        temp_seq = list()
        #for j in range(self.data_start, self.data_end+1-self.k):
        for j in range(start, min(end + 2 - self.k, 22799)): 
            curr_k_mer = self.k_mers_01[j]
            len_k_mer = len(curr_k_mer)
            #print j, curr_k_mer
            cur_cov = 0
            for k in self.k_mers[j]:
                #print "kmer detail", k, len(self.k_mers[j][k])
                cur_cov += len(self.k_mers[j][k])
            ############################################### case study
            #some case left, first and third bit-complement
            #second = first
            ############################################
            # very low coverage like 3
            #if ( ( len_k_mer >= 2 and tools.is_bool_reverse(curr_k_mer[0], curr_k_mer[1]) ) or len_k_mer == 1 or 
                #( len_k_mer==3 and  (tools.is_bool_reverse(curr_k_mer[0], curr_k_mer[1]) or 
                #                     tools.is_bool_reverse(curr_k_mer[0], curr_k_mer[2])) ) ):
            
            
            # most strict
            '''
            if len_k_mer >= 2:
                kmer1cov = len(self.k_mers[j][curr_k_mer[0]])
                kmer2cov = len(self.k_mers[j][curr_k_mer[1]])
            
            if ( len_k_mer >= 2 and tools.is_bool_reverse(curr_k_mer[0], curr_k_mer[1]) and kmer1cov + kmer2cov > 0.5*cur_cov ):
                #print kmer1cov
                #print kmer2cov
                #print "cur cov", cur_cov
                # for pacbio and nanopore, should I use different para
            '''

            #a bit higher coverage, need improve 
            if ( len_k_mer >= 2 and tools.is_bool_reverse(curr_k_mer[0], curr_k_mer[1])  ):
                if len(temp_seq) == 0:
                    h = Haplotype(j) 
                    temp_seq.extend(curr_k_mer[0])
                    h.update_clouds(curr_k_mer[0], self.k_mers[j])
                elif temp_seq[-self.k+1:] == list(curr_k_mer[0])[:-1]: # may have reads conflict 
                    temp_seq.append(curr_k_mer[0][-1])
                    h.update_clouds(curr_k_mer[0], self.k_mers[j])
                elif tools.is_bool_reverse( temp_seq[-self.k+1:], curr_k_mer[0][:-1] ):
                    temp_seq.append( tools.int_reverse( curr_k_mer[0][-1] ) )
                    h.update_clouds( tools.bool_reverse( curr_k_mer[0] ), self.k_mers[j] )
                else:
                    print j, "break type 2, the end of haplotype is not consistent"
                    #print temp_seq
                    #print temp_seq[-self.k+1:], list(curr_k_mer[0])[:-1]
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
                if len_k_mer == 0:
                    print j, "break type 3, no reads cover k_mer"
                    #print "check no reads cover ", self.k_mers_01[j]
                    #print "check no reads cover ", self.k_mers[j].keys()
                else:
                    print j, "break type 1, first and second most frequency k_mer not bit-complement"
                    #print self.k_mers_01[j]
                    #print self.k_mers[j].keys()
        if len(temp_seq) != 0:
            h.set_seq(temp_seq)
            haps.append(h)

        #self.unsure_clouds = set()
        for h in haps:
            if h.check_clouds_intersection() == False:
                #print "have intersection between left and right clouds"
                #h.printH()
                h.unsure_clouds.update( h.remove_intersection() )

        return haps  

    
    def fill_gap_between_haps(self, haps):
        
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
    
           
    def fill_two_haps(self, h1, h2):
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


    def merge_haps(self, haps):
        
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
        '''
        for h in haps:
            h.printH()
        '''
        
    def from_two_haps(self, h1, h2):
        assert h1.start < h2.start and h1.end < h2.end
        #assert h1.end < h2.start
        assert h1.check_clouds_intersection() # should not have
        assert h2.check_clouds_intersection() # intersection
        '''
        if h1.check_clouds_intersection() == False:
            print "h1 imposible"
            h1.printH()
        '''
        ll = len(h1.left_clouds.intersection(h2.left_clouds))
        lr = len(h1.left_clouds.intersection(h2.right_clouds)) 
        rl = len(h1.right_clouds.intersection(h2.left_clouds)) 
        rr = len(h1.right_clouds.intersection(h2.right_clouds)) 

        # one link: using reads overlap
        #assert h1.start + h1.len <= h2.start
        #if (ll > lr and rr > rl) or ( (ll>0 or rr>0) and (lr==0 or rl==0) ):
        # need improve
        if (ll > lr and rr > rl) or  (ll == lr and rr > rl) or (ll > lr and rr == rl) :
            print "merge two ", h1.start, h1.end, h2.start, h2.end
            print h1.seq,h2.seq
            print ll, lr, rr, rl
            if h1.end >= h2.start:
                if h1[ h2.start : h1.end+1 ] !=  h2[ h2.start : h1.end+1 ] : # need update, no need exactly =
                    print 'reads intersection conflict with overlap info, need imporve, overlap no need exactly correct '
                    return False
                temp_seq = h1.seq + h2[ h1.end+1: h2.end + 1 ]
            elif h1.end < h2.start:    
                temp_seq = h1.seq + [-1]*(h2.start-h1.start-h1.len) + h2.seq
            h1.set_seq( temp_seq )
            print "after merge two ", h1.start, h1.end
            print h1.seq
            h1.left_clouds.update(h2.left_clouds)
            h1.right_clouds.update(h2.right_clouds)
            if h1.check_clouds_intersection() == False:
                h1.unsure_clouds.update( h1.remove_intersection() )
            return True

        elif (lr > ll and rl > rr) or ( lr == ll and rl > rr) or ( lr > ll and rl == rr ):
            print "merge two ", h1.start, h1.end, h2.start, h2.end
            print h1.seq, h2.seq
            print ll, lr, rr, rl
            if h1.end >= h2.start: # assert h1.end < h2.end
                #   72      73     74  75 76
                #   h1.s          h1.e 
                #           h2.s          h2.e
                print "check", h1[ h2.start : h1.end+1 ] , h2[ h2.start : h1.end+1 ]
                print tools.list_reverse( h2[ h2.start : h1.end+1 ]) 
                if h1[ h2.start : h1.end+1 ] != tools.list_reverse( h2[ h2.start : h1.end+1 ] ):
                    print 'reads intersection conflict with overlap, need imporve, overlap no need exactly correct'
                    return False
                temp_seq = h1.seq + tools.list_reverse( h2[ h1.end+1 : h2.end+1] )
            elif h1.end < h2.start:    
                temp_seq = h1.seq + [-1]*(h2.start-h1.start-h1.len) + tools.list_reverse(h2.seq)
            
            h1.set_seq( temp_seq )
            print "after merge two ", h1.start, h1.end
            h1.left_clouds.update(h2.right_clouds)
            h1.right_clouds.update(h2.left_clouds)
            if h1.check_clouds_intersection() == False:
                h1.unsure_clouds.update( h1.remove_intersection() )
            return True
        else:
            # have overlap can't merge, break short part
            print "cant merge haps ll, lr, rr, rl", ll, lr, rr, rl
            if h1.end >= h2.start:
                temp_seq = h1[h1.start: h2.start]
                h1.set_seq( temp_seq ) # have overlap, cann't connect, remove the tail of h1
                print "reads intersection provide no info to merge, but overlap may helpful, need imporve "

        return False
                
