#!/usr/bin/env python3


class Contig:
    def __init__(self, name, seq):
        self.Name = name
        self.Seq = seq
        self.Len = len(seq) 
     
    def print(self):
        print (self.Name, self.Len)

def read_Contig(filename):
    contigs = []
    f = open(filename, "r")
    #for line in filter(lambda x: len(x) > 0, map(string, f)):
    for line in f:
        if line.startswith(">"):
            name = line.split()[0][1:]
        else:
            seq = line
            contig = Contig(name,seq)  
            contig.print()
            contigs.append(contig) 
    return contigs
