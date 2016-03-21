#!/usr/bin/env python
'''
Ian Rambo
Oyster Rocks Methylation Project

In-silico genome generator.
Generate randomized control "genomes" based off nucleotide frequencies
of a known genome.
Created: January 18, 2016
Last updated: February 12, 2016
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
sequenceLength = 550
#Number of FASTA sequences to write to output
nSeqs = 5000
#Count of created sequences
count = 0
#Number of in-silico genomes created
genCount = 0

rootDir = '/net/biohen29/data/imrambo'
genomeDir =  '%s/BacterialGenomes' % rootDir
outputDir = '%s/00-insilico_genomes' % rootDir
#outputFile = '%s/Anaeromyxobacter_dehalogenans_insilico_sequences_00.fa' % outputDir
#==============================================================================
#-----FUNCTIONS-----
#==============================================================================
def random_ntseq(seqLen, gc) :
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
    #Shuffle the nucleotides
    DNAlist = random.sample(DNA, len(DNA))
    DNA = ''.join(DNAlist)
    return(DNA)

def read_fasta(filename, fieldsep) :
    '''
    Read a FASTA file and store header and sequence for each contig as a 
    tuple containing the header a list, while the sequence
    is stored a string. 
    '''
    sequences = []
    header = None
    with open(filename, 'r') as fasta :
        line = fasta.readline()[:-1]
        while line :
            if line[0] == '>' :
                if header :
                    sequences.append((header, seq))
                header = line[1:].split(fieldsep)
                seq = ''
            else :
                seq += line
            line = fasta.readline()[:-1]
        sequences.append((header, seq))
    return(sequences)

def sliding_window(ntSeq, startpos, extent) :
    '''
    Run a sliding window on a nucleotide sequence to sample a genomic
    scaffold for GC content at certain intervals. 
    '''
    #Starting position
    position = startpos
    #How far the window extends from position in both directions
    windowExtent = extent
    #Full size of the sliding window
    #windowSize = (windowExtent * 2) + 1
    gcValues = []
    #Moving average window
    while position <= len(ntSeq) - windowExtent :
        window = ntSeq[(position - (windowExtent + 1)):(position + windowExtent)]
        #print(len(window))
        C = window.count('C')
        G = window.count('G')
        gc_percent = round(float(C + G) / len(window), 2)
        gcValues.append(gc_percent)
        position += 1
    

    #Get the minimum and maximum GC content
    minGC = min(gcValues)
    maxGC = max(gcValues)
    
    return(minGC,maxGC)
#==============================================================================
print('\n-----BEGIN-----\n')
#Create a new output directory if specified path does not exist
if not os.path.exists(outputDir) :
    print('\ncreating output folder %s\n' % outputDir)
    os.makedirs(outputDir)
else :
    print('\nnew directory not created; %s already exists\n' % outputDir)
    pass

#Specific accession numbers to extract from
targetAccessions = []
with open('/net/biohen29/data/imrambo/targetClassAccessions.txt', 'r') as tca :
    for line in tca :
        line = line.strip()
        lineL = line.split('\t')
        a = lineL[1]
        targetAccessions.append(a)
        
for subdir, dirs, files in os.walk(genomeDir) :
    for file in files :
        #Extract the accession number
        acc = re.findall(r'^(.*?)\.fna$', file)[0]
        if acc in targetAccessions :
            fastaFile = os.path.join(subdir, file)
            outputFile = '%s/%s_insilico_sequences.fna' % (outputDir, acc)
        
            #Read the genome scaffold file
            scaffold = read_fasta(fastaFile,'|')
            #Nucleotide sequence
            ntSequence = scaffold[0][1]
            
            #Nucleotide counts for scaffold sequence
            #scaffA = ntSequence.count('A')
            #scaffC = ntSequence.count('C')
            #scaffG = ntSequence.count('G')
            #scaffT = ntSequence.count('T')
            #
            #scaffoldGC = round(float(scaffC + scaffG) / len(ntSequence), 2)
            
            #nucSum = sum([A,C,G,T])
            
            #------------------------------------------------------------------------------
            print('\n.....CALCULATING.....\n')
            #Calculate GC content for 51 bp windows
            gcStats = sliding_window(ntSequence, 51, 50)
            #GC content range
            minGC = gcStats[0]
            maxGC = gcStats[1]
            print('\nGenerating random sequences of length %d\n' % sequenceLength)
            with open(outputFile, 'w') as OUT :
                for i in range(nSeqs) :
                    A = int()
                    C = int()
                    G = int()
                    T = int()
                    randomNT = str()
                    #random GC content within bounds of GC
                    gc = round(random.uniform(minGC, maxGC), 2)
                    #print(gc)
                    #Random nucleotide sequence
                    randomNT = random_ntseq(sequenceLength, gc)
                    
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
                    
        
        else :
            pass
print('minGC = %.2f, maxGC = %.2f' % (minGC, maxGC))
print('\n%d sequences of length %d written to %s\n' % (count, sequenceLength, outputFile))
print('\n-----DONE-----\n')


