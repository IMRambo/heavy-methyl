#Ian Rambo
#University of Delaware
#BINF 868
#Last updated: September 20, 2015
#==============================================================================
import os
import re
os.chdir('/Users/imrambo/Documents/MS_thesis')

def read_FASTA(filename) :
    '''
    Read a FASTA file and store header and sequence for each contig as a list
    within a list of all contigs.
    The header string is split into a list, while the sequence
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
                header = line[1:].split('|')
                seq = ''
            else :
                seq += line
            line = fasta.readline()[:-1]
        sequences.append((header, seq))
    return(sequences)
#------------------------------------------------------------------------------
'''
GOAL: find motifs 
'''
def find_motifs(pattern, genome):
    '''
    Find the number of potential methylation motifs within a genome.
    '''
    if type(genome) is list :
        for i in range(len(genome)) :
            if type(genome[i]) is str and type(pattern) is str :
                print(len(re.findall(pattern, genome[i])))
            elif type(genome[i]) is str and type(pattern) is list :
                pattern = '|'.join(pattern)
                print(len(re.findall(pattern, genome[i])))
            else :
                print('Genome sequence is not a string')
                break

    elif type(genome) is str :
        print(len(re.findall(pattern, genome)))

    elif type(genome) is not str :
        print('Genome sequence is not a string')
#==============================================================================
sequences0306 = read_FASTA('S03-06cm_pkhAnnotated.fa')
sequences1215 = read_FASTA('S12-15cm_pkhAnnotated.fa')
sequences2427 = read_FASTA('S24-27cm_pkhAnnotated.fa')

def count_motifs(sequences) :
    counts = []
    if type(sequences) is list :
        for i in range(len(sequences)) :
            seq = sequences[i][1]
            if type(seq) is str :
                ccgg = seq.count('CCGG')
                gatc = seq.count('GATC')
                gantc = len(re.findall(r'GA[GATC]TC', seq))
                counts.append([ccgg, gatc, gantc])
    return counts
a = count_motifs(sequences0306)
print(a)