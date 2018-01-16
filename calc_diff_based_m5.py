#!/usr/bin/env python

import os
import sys
import string
from calc_switch import Sequence
import calc_switch  
from calc_switch import Alignment


align = calc_switch.read_blasr_m5(sys.argv[1])
align.PreProcess()
align.left.generateSeqPos()
align.right.generateSeqPos()
align.printAlign() 
align.printAlignDiff(5)
align.getAlignDiffStat(5)
