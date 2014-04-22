#!/usr/bin/env python2.7
import sys
import pdb
import argparse
from pandas import *
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def getMutatedPositions():
    random.seed(100)
    df = read_table("depthStats", header= None, index_col=1)
    mutationPos = []
    nPos = 100 # number of positions we want to sample
    mPos = [] # mutated Positions
    flags = []
    for i in range(2, 11):
        pos = random.sample(df[df[2] == i].index.tolist(), nPos)
        flag = [i] * nPos
        mPos += pos
        flags += flag
    #mDf = DataFrame(data={'mPos':mPos, 'coverage':flags}, index=mPos)
    mDf = Series(data=flags, index=mPos)
    mDf.to_csv("mutatedPositions.csv")
    return mDf

def mutation(originalNt):
  d = {
      'A': 'C',
      'C': 'G',
      'G': 'T',
      'T': 'A',
      'N': 'N'
      }
  return d[originalNt]

if __name__ == '__main__':
    o = sys.stdout
    e = sys.stderr
    parser= argparse.ArgumentParser(description="")
    parser.add_argument("file", help="fasta file to be mutated from")
    #parser.add_argument("num", type = int,
            #help="number of SNPs you want to mutate. default = 1000", 
            #default = 1000)
    args = parser.parse_args() 
    rec = SeqIO.parse(args.file, 'fasta').next()
    np.random.seed(10)
    mutations = len(rec.seq)/300
    zeroBasedPos = np.random.randint(0, len(rec.seq), mutations)
    seq = Series(list(rec.seq))
    #seq2 = seq.copy()
    seq[zeroBasedPos] = seq[zeroBasedPos].map(mutation)
    # calculate all mutated positions excluding those positions, where 'N's reside
    mutatedPosNotN = zeroBasedPos[np.array(seq[zeroBasedPos] != 'N')]
    #combined = concat([seq, seq2], axis = 1) 
    #combined.columns = ['original', 'mutated']
    #e.write("nucleotide composition before mutation\n")
    #print combined.ix[combined['original'] != combined['mutated'], 'original'].value_counts()
    #e.write("nucleotide composition after mutation\n")
    #print combined.ix[combined['original'] != combined['mutated'], 'mutated'].value_counts()
    #combined.to_csv("mutatedPositions.csv") 
    mseq = "".join(seq)
    mrec = SeqRecord(Seq(mseq), id='Chr4', description="")
    SeqIO.write(mrec, "pChr4.fasta", 'fasta')
    np.save("mutatedPos", mutatedPosNotN)
    # zeroBasedPos = np.load("mutatedPos" + ".npy")
