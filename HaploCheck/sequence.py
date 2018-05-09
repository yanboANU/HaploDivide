#!/usr/bin/env python

import os
import sys
import string
from itertools import ifilter,imap
import collections


# for alignment
class Sequence(object):

    def __init__(self, name, length, s, e, seq, direction):
        self._name = name
        self._length = int(length)
        self._s = int(s)
        self._e = int(e)
        self._seq = seq   # ATCG- 
        self._seq_pos = [] 
        self._dir = direction # + or - 
        # _seq_pos[i] = j , the position in self sequence(have no -) j is reposible to i in alignSeq(self.seq) 
    
    def _print_seq(self):
        print (self._seq)
    
    def _generate_seq_pos(self):
        #align seqence to position sequence in original seq 
        if self._dir == '+':
            ss = self._s
            for ii in self._seq: 
                if ii !=  '-':
                    self._seq_pos.append(ss)
                    ss += 1 
                else: 
                    self._seq_pos.append(-1)

        if self._dir == '-':
            ss = self._s
            for ii in self._seq:
                if ii !=  '-':
                    self._seq_pos.append(self._length-1 - ss)
                    ss += 1 
                else: 
                    self._seq_pos.append(-1)

    def _print_seq_part(self,start,ll):
        print self._seq[start:start+ll] 

    def _write_name(self, f):
        f.write("%s %s %s %s"  % (self._name, self._length, self._s, self._e))

    def _write_seq(self, f):
        f.write(self._seq) 

