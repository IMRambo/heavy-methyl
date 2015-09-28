#Function for reading FASTA files. 
#Ian Rambo
#University of Delaware
#Last updated: September 20, 2015

import os
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

seqs0306 = read_FASTA('S03-06cm_pkhAnnotated.fa')
seqs1215 = read_FASTA('S12-15cm_pkhAnnotated.fa')
seqs2427 = read_FASTA('S24-27cm_pkhAnnotated.fa')



#for s in range(len(sequences)) :
#    if any('Deltaproteobacteria' in a for a in sequences[s][0]) :
#        if any('ase' in e for e in sequences[s][0]) :
#            print(sequences[s][0], sequences[s][1])
