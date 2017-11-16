# coding: utf-8

import gensim
import csv
import numpy as np
from numpy import linalg as LA
import sys

model = gensim.models.Word2Vec.load('/home/garner1/Work/dataset/rest2vec/modelPlus')

# infile = '/home/garner1/Work/dataset/rest2vec/restseq_plus/ACGACATC/chr10.loc-docs'
infile = sys.argv[1]
with open(infile, 'r') as f:
    reader = csv.reader(f, delimiter='|')
    signal = []
    for row in reader:
        sentence_experiment = row[1].split()
        sentence_reference = row[2].split()
        sumWV_experiment = [0]*100
        sumWV_reference = [0]*100
        count_experiment = 0.0
        count_reference = 0.0
        for word_pairs in zip(sentence_experiment,sentence_reference):
            if word_pairs[0] in model.wv.vocab: 
                count_experiment += 1.0
                sumWV_experiment += model.wv[word_pairs[0]]
            if word_pairs[1] in model.wv.vocab: 
                count_reference += 1.0
                sumWV_reference += model.wv[word_pairs[1]]
        v1 = np.asarray(sumWV_experiment)/count_experiment
        v2 = np.asarray(sumWV_reference)/count_reference
        signal.append( LA.norm(v1-v2) )
    outfile = infile + '.signal'
    np.savetxt(outfile, signal, fmt='%6.6f')    

