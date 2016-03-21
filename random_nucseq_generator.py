#!/usr/bin/env python
'''
Ian Rambo
Oyster Rocks Methylation Project

Random nucleotide sequence generator.
Generate a FASTA file of nucleotides with a certain frequency.
Created: January 15, 2016
Last updated: January 18, 2016
'''
import random
from random import choice
import os
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
#-----GLOBAL VARIABLES-----
#Length of each nucleotide sequence
sequenceLength = 420
#Number of FASTA entries to write to output
nSeqs = 2000
#Fixed GC content
#gc = 0.6711
#Count of created sequences
count = 0
outputDir = '/Users/imrambo/Documents/BINF868/00-random_sequences'
outputFile = 'random_sequences_00.fa'
#==============================================================================
#-----FUNCTIONS-----
#==============================================================================
def ntString(seqLen, gc) :
    '''
    Create randomized nucleotide strings of a specific length and GC content.
    '''
    DNA = str()
    #GC count
    gcCount = int(seqLen * gc)
    #AT count
    atCount = seqLen - gcCount
    for count in range(gcCount) :
        DNA += choice('CG')
    for count in range(atCount) :
        DNA += choice('AT')
    DNAlist = random.sample(DNA, len(DNA))
    DNA = ''.join(DNAlist)
    return DNA
#==============================================================================
#Create a new output directory if specified path does not exist
if not os.path.exists(outputDir) :
    print('\ncreating output folder %s\n' % outputDir)
    os.makedirs(outputDir)
else :
    print('\nnew directory not created; %s already exists\n' % outputDir)
    pass

os.chdir(outputDir)

print('\nGenerating random sequences of length %d\n' % sequenceLength)
with open(outputFile, 'w') as OUT :
    for i in range(nSeqs) :
        randomNT = str()
        #random GC content
        gc = round(random.uniform(0.0, 1.0), 2)
        #Random nucleotide sequence
        randomNT = ntString(sequenceLength, gc)
        
        #Nucleotide counts
        A = randomNT.count('A')
        C = randomNT.count('C')
        G = randomNT.count('G')
        T = randomNT.count('T')
        nucSum = sum([A,C,G,T])
        
        count += 1
        
        if nucSum == sequenceLength :
            fastaEntry = '>random_%d;gc=%.2f;length=%d\n%s\n' % (count, gc,\
                                                     len(randomNT), randomNT)
            #print(fastaEntry)
            OUT.write(fastaEntry)
        else :
            break
            print('Sequence length is not the same as desired length\n')
            print('Terminating...\n')
    
print('\n%d sequences of length %d written to %s\n' % (count, sequenceLength, outputFile))
print('\n-----DONE-----\n')
