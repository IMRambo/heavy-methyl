#!/usr/bin/env python
'''
Extract sequences from a FASTA file.

Created: February 25, 2016
Last modified: March 3, 2016
'''
#==============================================================================
#FUNCTIONS
#==============================================================================
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

#==============================================================================
dropboxDir = '/Users/imrambo/Dropbox/OysterRocksMethylationProject'
#contigID = '%s/methylContigID-27.txt' % dropboxDir
inputFasta = '%s/methylContigSeqs-15.fa' % dropboxDir
outFile = open('%s/methylContigSeqs-targetClass-15.fa' % dropboxDir, 'w')
outFile2 = open('%s/methylContigID-targetClass-15.txt' % dropboxDir, 'w')
nContigs = 0

print('\n-----BEGIN-----\n')

targetClass = ['Actinobacteria','Alphaproteobacteria','Bacilli',\
               'Betaproteobacteria','Clostridia','Deinococci',\
               'Deltaproteobacteria','Gammaproteobacteria']

#conIDList = []
#with open(contigID, 'r') as contig :
#    for c in contig :
#        c = c.rstrip()
#        conIDList.append(c)
        

seqs = read_fasta(inputFasta, '|')
for i in range(len(seqs)) :
    header = seqs[i][0]
    sequence = seqs[i][1]
    if header[9] in targetClass :
        outFile.write('>%s\n%s\n' % ('|'.join(header), sequence))
        outFile2.write(header[0] + '\n')
        nContigs += 1
       
    else :
        pass
outFile.close()
outFile2.close()

print('wrote %d contigs' % nContigs)
print('\n-----END-----\n')