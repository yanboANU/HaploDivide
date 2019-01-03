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
            blocks.append( (pos[-1], indexLen) )
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
        self.enumerate_length_threshold = 11 # trade-off between accuracy and time 

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


        #sys.exit()
        
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
        best_hap = []
        dis = end - start + 1
        allPosSeq = tools.enumerate_01_list(dis)    
        l = len(allPosSeq)/2
        clouds = cloud.get_clouds(start, end, clouds_at_index)
        
        for seq in allPosSeq[0:l]:
            h = haplotype.Haplotype(start)        
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
            print "enumerate block more than 1"
            seqs = [] 
            for h in best_hap:
                print h.seq
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
            h = haplotype.Haplotype(start)        
            h.set_seq(temp_seq)
            return [h]

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

    def generate_distribution(self, haps):
        
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
            haps = self.generate_haps_from_k_mers(start, end, label)
            label -= 1 # loose bit-complement constraint
 
        print "after generate haps from kmers, haps number ", len(haps)
        hapsNum = len(haps)
        if hapsNum == 0:
            haps = self.enumerate_seq(start, end, clouds_at_index)
        elif hapsNum >= 1:
            #self.dp(haps)
            print "old haps number", len(haps)
            haps = self.generate_newhaps_loose_bit_constraint_locally(haps, start, end, label)
            print "new haps number", len(haps)
        # optial 
        if haps[0].start - 1 > start:   # enumerate start
            haps = self.enumerate_seq(start, haps[0].start-1, clouds_at_index) + haps
        if haps[-1].end + 1 < end:      # enumerate end
            haps = haps + self.enumerate_seq(haps[-1].end+1, end, clouds_at_index)
        
        if len(haps) > 1:
            self.enumerate_gap_between_haps(haps, clouds_at_index)
        print "block region", start, end, "has", len(haps), "haplotype segements"   
        return haps

    def generate_newhaps_loose_bit_constraint_locally(self, haps, start, end, label):
        
        print "old haps in subfunction, number", len(haps)
        '''
        for h in haps:
            print h.start, h.end
        '''    
        subHaps = {}
        hapsNum = len(haps)
        for i in range(hapsNum+1): 
            s = start if i==0 else haps[i-1].end + 1   
            e = end if i==hapsNum else haps[i].start - 1
            subLabel = label
            hs = list()
            if (e - s + 1 >= self.enumerate_length_threshold) and subLabel >= 1: # e>s, error rate high, 
                while (len(hs) == 0 and subLabel >= 1):
                    hs = self.generate_haps_from_k_mers(s, e, subLabel)
                    subLabel -= 1
                if len(hs) > 0:
                    subHaps[i] = hs
        
        if len(subHaps) > 0:    
            '''
            print "subHaps"
            for k in subHaps:
                for h in subHaps[k]:
                    print h.start, h.end
            '''     
            newHaps = list()        
            for i in range(hapsNum+1):
                if i in subHaps:
                    newHaps.extend(subHaps[i])
                if i < hapsNum:
                    newHaps.append(haps[i])
            
            print "new haps in subfunction, number", len(newHaps)
            '''
            for h in newHaps:
                print h.start, h.end
            '''
            haps = newHaps        
        return haps


        

    def run_dc(self, clouds_at_index): # 12, Dec

        #print self.find_position_no_reads_span(0, len(clouds_at_index), clouds_at_index)
        #for block in self.blocks:
            #print block
        #sys.exit() 

        '''
        print_aligned_region(3695, 3705, clouds_at_index)
        print_aligned_region(8380, 8390, clouds_at_index)
        print_aligned_region(14955, 14965, clouds_at_index)
        #print_aligned_region(3670, 3680, clouds_at_index)
        sys.exit()
        '''
        allHaps = []
        count = 1
        print "There are", len(self.blocks), "blocks"
        for block in self.blocks:  # blo/cks break at no reads cover region
            print "phasing block: ", count, ", region ", block[0] , block[-1]
            allHaps.extend(self.run_dc_in_block(block[0], block[-1], clouds_at_index))
            count += 1

        self.generate_distribution(allHaps)
        #sys.exit()


        print "before finall enmuerate ", len(allHaps)    
        self.enumerate_gap_between_haps(allHaps, clouds_at_index) # only need work for between block, need to updata
                                                                  # for two consective hap, local region , two MEC
                                                                  # longer , maybe have MEC
                                                                  
        print "after finall enmuerate ", len(allHaps)    

        fout1 = open("hap.out" , "w")
        for i, block in enumerate(allHaps):
            fout1.write( "%d %d\n" % (block.start, block.end) )
            '''
            for k,y in enumerate(block.seq):
                j = block.start + k
                if y != -1:
                    fout1.write("%d\t%d\t%d\t\n" % (j, y, tools.int_reverse(y)))
            fout1.write("********\n")
            '''
        fout1.close()

        return allHaps  


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
        
        if dis >= self.enumerate_length_threshold: # need improve
            print "dis is", dis, "too long to enumerate"
            print len(clouds)

            cloudsMid = list(cloud.get_clouds(h1.end+1, h2.start-1, clouds_at_index))
            cloudsNum = len(cloudsMid)
            '''               
            if cloudsNum < 10:
                print "mid clouds"
                for c in cloudsMid:
                    print c.start, c.end, c.seq
                    
                h = enumrate_clouds(h1.end+1, h2.start-1, cloudsMid)
                
                print h1.end+1, h2.start-1
                #temp_seq = cloudsMid[0][h1.end+1 : h2.start-1]
                for i in range(1, cloudsNum):
                    #if is_bool_equal(cloudsMid[i][h1.end+1 : h2.start-1]) or is_bool_reverse2()   
                sys.exit() 
            '''    


            '''
            b = range(h1.end+1, h2.start-1)
            h = self.hmm.run_viterbi(b)
            print "after run viterbi"
            print h
            sys.exit()
            '''
            return False
        
        print h1.start, h2.end, "this region supported reads number", len(clouds)    
        allPosSeq = tools.enumerate_01_list(dis)    

        for seq in allPosSeq:
            for i in range(0,2): # debug maybe h1.seq + seq + tools.reverse(h2.seq)
                if i == 0:
                    temp_seq = h1.seq + seq + h2.seq 
                elif i == 1:
                    temp_seq = h1.seq + seq + tools.list_reverse(h2.seq) # debug maybe h1.seq + seq + tools.reverse(h2.seq)

                h3 = haplotype.Haplotype(h1.start)
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
            unsure_pos = idx.index(False)
            idxF = idx.count(False)


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
            ###########
            #000100
            #000000
            #TTTFTT
            ###########
            if idxF == 1 or idxF <= min(len(idx)*0.1, 2):
                temp_seq = best_hap[0].seq[ : unsure_pos]
                temp_seq = temp_seq + [-1]
                temp_seq = temp_seq + best_hap[0].seq[unsure_pos+1 : ]
                h1.set_seq(temp_seq)
                print "2 enumerate", h1.seq
                return True
           
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
            #fout.write()
            start = best_hap[0].start
            end = best_hap[0].end + 1
            print best_hap[0].seq, best_hap[0].start, best_hap[0].end

            ##########
            #case 1
            ########## 
            seqs = [] 
            for h in best_hap:
                print h.seq
                seqs.append(h.seq)
            arr = np.sum(seqs, axis=0)    
            temp_seq = []
            count = 0
            for ele in arr:
                if ele == 0:
                    temp_seq.append(0)
                elif ele == l:
                    temp_seq.append(1)
                else:
                    temp_seq.append(-1)
                    count += 1
            print temp_seq, float(count)/len(best_hap[0].seq)
            #if (count == 1 or (count ==2 and float(count)/len(best_hap[0].seq)<0.05) ):
            #if (count == 1 or  float(count)/len(best_hap[0].seq)<0.07) :
            if (temp_seq[0] != -1 and temp_seq[-1] != -1):
                h1.set_seq(temp_seq)
                return True

            #########
            #case 2
            #########
            for i in range(h2.start-1, h1.end-1, -1):
                temp_seq = best_hap[0][start:i]
                #print "temp_seq", temp_seq
                for j in range(1,l):
                    #print best_hap[j][start:i]
                    if temp_seq == best_hap[j][start:i] or tools.is_bool_reverse2(temp_seq, best_hap[j][start:i]):
                        continue
                    else:
                        break
                if j >= l-1 and (temp_seq == best_hap[-1][start:i] or tools.is_bool_reverse2(temp_seq, best_hap[-1][start:i]) ):
                    #print i
                    break
            #print "left seq"
            #print temp_seq # this right seq
            h1.set_seq(temp_seq)

            for i in range(h1.end+1, h2.start+1):
                temp_seq = best_hap[0][i:end]
                #print "temp_seq", temp_seq
                for j in range(1,l):
                    #print best_hap[j][i:end]
                    if temp_seq == best_hap[j][i:end] or tools.is_bool_reverse2(temp_seq, best_hap[j][i:end]):
                        continue
                    else:
                        break
                if j >= l-1 and (temp_seq == best_hap[-1][i:end] or tools.is_bool_reverse2(temp_seq, best_hap[-1][i:end]) ):
                    #print i
                    break
            #print "right seq"
            #print temp_seq # this right seq
            h2.start = i
            h2.set_seq(temp_seq)
            print h1.seq, h2.seq
            return False
            
        return False



    def calculate_MEC(self, h1, clouds):
        MEC = 0
        h2 = haplotype.Haplotype(h1.start)
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


    def generate_haps_from_k_mers(self, start, end, label):  # need to updata, different coverage, different fix_rule
        
        haps = list()
        temp_seq = list()
        #for j in range(self.data_start, self.data_end+1-self.k):
        for j in range(start, min(end + 2 - self.k, 22799)): 
            boolKmers = self.k_mers_01[j]
            #print j, boolKmers
    
            if self.check_bit_complement(boolKmers, self.k_mers[j], label):
                if len(temp_seq) == 0:
                    h = haplotype.Haplotype(j) 
                    temp_seq.extend(boolKmers[0])
                    h.update_clouds(boolKmers[0], self.k_mers[j])
                elif temp_seq[-self.k+1:] == list(boolKmers[0])[:-1]: # may have reads conflict 
                    temp_seq.append(boolKmers[0][-1])
                    h.update_clouds(boolKmers[0], self.k_mers[j])
                elif tools.is_bool_reverse( temp_seq[-self.k+1:], boolKmers[0][:-1] ):
                    temp_seq.append( tools.int_reverse( boolKmers[0][-1] ) )
                    h.update_clouds( tools.bool_reverse( boolKmers[0] ), self.k_mers[j] )
                else:
                    print j, "kmer break type 2, not satisfy consistency"
                    #print temp_seq
                    #print temp_seq[-self.k+1:], list(boolKmers[0])[:-1]
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
                if len(boolKmers) == 0:
                    print j, "kmer break type 3, no reads cover all position in k_mer"
                    #print "check no reads cover ", self.k_mers_01[j]
                    #print "check no reads cover ", self.k_mers[j].keys()
                else:
                    print j, "kmer break type 1, not satisfy bit-complement"
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
    
    ############################################### case study
    #some case left, first and third bit-complement
    #second = first
    ############################################
    # very low coverage like 3
    
    def check_bit_complement(self, boolKmers, kmers, label): # most loose

        lenKmer = len(boolKmers) # bool_k_mers a list of kmer which not include -1
        if label == 1:
            if ( ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1] ) ) or 
                    lenKmer == 1 or 
                    ( lenKmer==3 and  (tools.is_bool_reverse(boolKmers[0], boolKmers[1]) or 
                                         tools.is_bool_reverse(boolKmers[0], boolKmers[2])) ) ):
                return True
            return False

        #def check_bit_complement2(self, boolKmers, kmers):
        if label == 2:
            #a bit higher coverage, need improve 
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1])  ):
                return True
            return False

        #def check_bit_complement3(self, boolKmers, kmers): # most strict 
        if label == 3:
            curCov = 0
            for k in kmers:
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
                
