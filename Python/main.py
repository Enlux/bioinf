#!/usr/bin/python3
from pybloomfilter import BloomFilter
from GraphTraversal import GraphTraversal

#size of the kmer
inputFilename = "read50x_ref10K_e001.fasta"
#inputFilename = "EColi-synthetic/ecoli_test_reads.fasta"
outputFilename = "test.txt"
kmerSize = 27

#init bloom filter TODO remove it before
bloom = BloomFilter(10000000, 0.1, b'filter.bloom')


#open fasta file
f = open(inputFilename)

print("Adding kmers to bloom")
#KISS approach, read it all
lines = f.readlines()
kmers = []
for line in lines:
    if line[0] != '>':
        l = line.rstrip()
        le = len(line)
        for i in range(kmerSize, le):
            #input into the bloom
            kmer = line[i - kmerSize:i]#.encode()
            kmers.append(kmer)
            #kmer = kmer.encode()
            #bloom.add(kmer)

print("Counting k-mers")
kmerCounting = {}
for kmer in kmers:
    if kmer not in kmerCounting:
        kmerCounting[kmer] = 1
    else:
        kmerCounting[kmer] += 1

print("Filtering k-mers")
kmers = []
for kmer in kmerCounting:
    if kmerCounting[kmer] >= 3:
        kmers.append(kmer)
        bloom.add(kmer.encode())

print("Len filtered k-mers:", len(kmers))

print("Quering bloom for extension nodes")
#store set of extensions where bloom answers yes
P = []
extLetters = ['T', 'C', 'A', 'G']
for kmer in kmers:
    for e in extLetters:
        extKmer = (kmer[1:] + e)
        if extKmer.encode() in bloom:
            P.append(extKmer)

P = list(set(P))

#just a check
print("Cheking length of kmers dataset vs P (from bloom)")
print(len(kmers), len(P))

print("Building cFP")
#build cfp structure, again KISS approach of algo (page 3, algorithm 1)
D = P
cFP = []

print("TODO uncomment cFP")
'''
for m in D:
    itsIn = False
    for kmer in kmers:
        if m == kmer:
            itsIn = True
            break
    if not itsIn:
        cFP.append(m)

print(cFP)
'''

result = GraphTraversal(bloom, cFP, kmers, kmerSize, outputFilename)
