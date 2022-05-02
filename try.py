import numpy as np
import pandas as pd

def readFa(file):
    fa = open(file, "r")
    seqs = []
    seq = ""
    for line in fa:
        if not line.startswith(">"):
            seq += line[:-1]
        if len(seq) >= 2000:
            seqs.append(seq + "$")
            seq = ""
    return seqs

def readFq(file):
    fq = open(file, "r")
    reads = []
    for line in fq:
        if not line.startswith(">"):
            reads.append(line)
    return reads

def rotate(sequence):
    seqs = np.array([[sequence]])
    for i in range(0, len(sequence)):
        sequence = sequence[1:len(sequence)] + sequence[0]
        seqs = np.insert(seqs, len(seqs), values=sequence, axis=0)
    seqs = np.insert(seqs, 1, values=range(0, (len(seqs))), axis=1)
    seqs = seqs[seqs[:, 0].argsort()]
    return seqs

def BWT(seq):
    seqs = rotate(seq[0])
    ref = np.array([seqs[:, 0][-1]])
    s = [list(seqs[:, 1])]
    for i in seq[1:len(seq)]:
        seqs = rotate(i)
        ref = np.insert(ref, len(ref), values=[seqs[:, 0][-1]], axis=0)
        s.append(seqs[:, 1])
    return ref, s

def get_seed(read, seed: int, gap: int):
    FL, RL = [], []
    read = read[:-1]
    rr = read[::-1].lower()
    rr = rr.replace("a", "T")
    rr = rr.replace("c", "G")
    rr = rr.replace("g", "C")
    rr = rr.replace("t", "A")
    for i in range(0, len(read), gap):
        FL.append(read[i: i + seed])
        RL.append(rr[i: i + seed])
    return FL, RL

'''def match(seed, ref):'''




'''def build(file):



def align(file, seed, gap, ):'''


fa = readFa("../data/Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa")
FL, RL = get_seed(fa[0], 50, 10)
print(FL)
print(RL)
'''print(readF("../data/Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa"))'''
