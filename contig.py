#!/usr/bin/env python3


class Contig:
    def __init__(self, name, seq):
        self.Name = name
        self.Seq = seq
        self.Len = len(seq) 
     
    def print(self):
        print (self.Name, self.Len)

def read_Contig(filename):
    contigs = {}
    f = open(filename, "r")
    #for line in filter(lambda x: len(x) > 0, map(string, f)):
    seq = ""
    for line in f:
        if line.startswith(">"):
            if len(seq) > 0:
                contig = Contig(name,seq)  
                contig.print()
                contigs[name] = contig
            name = line.split()[0][1:]
            seq = "" 
        else:
            seq = seq + line.strip()
            
    if len(seq) > 0:
        contig = Contig(name,seq)  
        contig.print()
        contigs[name] = contig
    f.close()        
    return contigs
