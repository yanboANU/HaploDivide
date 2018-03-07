import sys

def read_mutation_list(filename):
    f = open(filename, "r")

    pos = []  
    for line in f:  
        words = line.split()
        if len(words) == 1 or len(words) == 3:
            pos.append(int(words[0]))
    
    return pos


def read_columns(filename):
    f = open(filename, "r")
    columns = {}
    for line in f:
        if line.startswith("reference"):
            words = line.split() 
            pos = int(words[1]) 
            assert not in columns
            columns[ words[1] ] = []
        else: 
             
