#!/usr/bin/env python
#Ian Rambo
#BINF868
#University of Delaware
#Assignment 1: methylation motifs
#Last updated: September 9, 2015
#Howdy Doody
#==============================================================================
'''

DETERMINE POTENTIAL METHYLATION SITES IN GENOMES OF INTEREST
GOALS:
Count potential sites for cytosine and adenine methylation in genomes
Determine methyltransferase genes in genomes

'''

import re
import os
from Bio import Entrez, SeqIO
#==============================================================================
#FUNCTIONS

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
#TARGET GENOMES
#Download complete genomes from NCBI RefSeq.

'''
Dehalococcoides mccartyi strain SG1 genome
Illumina Hi-Seq 2000 shotgun sequencing
Assembled in SOAPdenovo, scaffolds constructed in Opera
Wang et al. 2014
'''
#GI numbers for scaffolds 1-5
Dmccartyi_GI = ['736351273', '36350382', '736350119', '736350202', '736350057']
#Build complete D. mccartyi genome from scaffolds
Dmccartyi_genome = list()

for i in range(len(Dmccartyi_GI)) :
    Entrez.email = 'imrambo@udel.edu'
    scaffold_h = Entrez.efetch(db = 'nucleotide',
                                        id = Dmccartyi_GI[i],
                                        rettype = 'fasta')
    scaffold = SeqIO.read(scaffold_h, 'fasta')
    scaffold_h.close()
    Dmccartyi_genome.append(str(scaffold.seq))
#------------------------------------------------------------------------------
'''
Desulfocapsa sulfexigens DSM 10523 genome
454 GS FLX Titanium pyrosequencing
Assembled in Newbler v2.5.3, gaps in assembly closed by 259 Sanger reads
Finster et al. 2013
'''
Entrez.email = 'imrambo@udel.edu'
Dsulfexigens_h = Entrez.efetch(db = 'nucleotide',
id = '451945650',
rettype = 'fasta')

Dsulfexigens_fasta = SeqIO.read(Dsulfexigens_h, 'fasta')
Dsulfexigens_h.close()
Dsulfexigens_genome = str(Dsulfexigens_fasta.seq)
#------------------------------------------------------------------------------
'''
Desulfovibrio desulfuricans ND132 complete genome
Whole genome shotgun sequencing
Assembly in Newbler v.2.3 (pre-release)
'''
Entrez.email = 'imrambo@udel.edu'
Desulfovibrio_h = Entrez.efetch(db = 'nucleotide',
id = '376294792',
rettype = 'fasta')

Desulfovibrio_fasta = SeqIO.read(Desulfovibrio_h, 'fasta')
Desulfovibrio_h.close()
Desulfovibrio_desulfuricans_genome = str(Desulfovibrio_fasta.seq)
#==============================================================================
#BACTERIAL METHYLTRANSFERASES

#Fetch adenine MTase sequence files from NCBI's Entrez databases
#File handles
#Cell cycle regulated m6A MTase (CcrM), Caulobacter crescentus
Entrez.email = 'imrambo@udel.edu'
CcrM_h = Entrez.efetch(db = 'nucleotide',
                       id = '221232939',
                       rettype = 'fasta',
                       seq_start = '398637',
                       seq_stop = '399713')
#Dam m6A MTase, E.coli K12
Entrez.email = 'imrambo@udel.edu'
Dam_h = Entrez.efetch(db = 'nucleotide',
                    id = '556503834',
                    rettype = 'fasta',
                    seq_start = '3515913',
                    seq_stop = '3515077')
#Dcm m5C MTase, E.coli K12
Entrez.email = 'imrambo@udel.edu'
Dcm_h = Entrez.efetch(db = 'nucleotide',
                    id = '556503834',
                    rettype = 'fasta',
                    seq_start = '2032317',
                    seq_stop = '2030899')

CcrM = SeqIO.read(CcrM_h, 'fasta')
Dam = SeqIO.read(Dam_h, 'fasta')
Dcm = SeqIO.read(Dcm_h, 'fasta')

CcrM_h.close()
Dam_h.close()
Dcm_h.close()

MTDict = {'CcrM':CcrM.seq,
          'Dam':Dam.seq,
          'Dcm':Dcm.seq}
#==============================================================================
#Count potential methylation motifs in target genomes
print('CCGG counts, Dehalococcoides mccartyi')
find_motifs('CCGG', Dmccartyi_genome)
print('GANTC counts, Dehalococcoides mccartyi')
find_motifs('GA[TCAG]TC', Dmccartyi_genome)
print('GATC counts, Dehalococcoides mccartyi')
find_motifs('GATC', Dmccartyi_genome)

print('CCGG, Desulfocapsa sulfexigens genome')
find_motifs('CCGG', Dsulfexigens_genome)
print('GANTC, Desulfocapsa sulfexigens genome')
find_motifs('GA[TCAG]TC', Dsulfexigens_genome)
print('GATC, Desulfocapsa sulfexigens genome')
find_motifs('GATC', Dsulfexigens_genome)
#==============================================================================
myMotifs = ['CCGG', 'GA[TCAG]TC', 'GATC']
find_motifs(myMotifs, Dmccartyi_genome)
