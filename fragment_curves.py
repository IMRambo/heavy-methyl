#!/usr/bin/env python
#Ian Rambo
#Created: January 9, 2016
#Last modified: February 28, 2016
import os
import re
import datetime
from time import strftime

'''
GOAL: run DNA curvature model script for gene fragment sequences.
DNA curvature script by Christoph Gohlke
www.com.com
'''
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
#FUNCTIONS
#==============================================================================
def read_fasta(filename, fieldsep) :
    '''
    Read a FASTA file and store header and sequence for each contig as a 
    tuple containing the header as list, while the sequence
    is stored as a string. 
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
#GLOBAL VARIABLES
rootDir = '/net/biohen29/data/imrambo/'
#rootDir = '/Users/imrambo/Documents/BINF868'
inputDir = '%s/00-ExtractData' % rootDir
outputRoot = '%s/01-remainingCurveData' % rootDir
todayDate = strftime('%d-%m-%Y')
#File with script status reports
statusFile = open('%s/curve_status-%s.txt' % (rootDir,todayDate), 'w')
#File containing target accession numbers
accFile = '%s/remainingCurveAccessions.txt' % rootDir
#Desired sequence length of input fragments
seqLen = 500
#Maximum sequence length - don't change this
maxLen = 510
#Count of genomes processed
genCount = 0
#==============================================================================
##-----MAIN-----##
#print('\n-----BEGIN-----\n')
statusFile.write('-----BEGIN @ %s-----\n\n' % strftime('%d-%m-%Y %H:%M:%S'))

if not os.path.exists(outputRoot) :
    statusFile.write('\ncreated output folder %s\n' % outputRoot)
    os.makedirs(outputRoot)
else :
    pass

#for subdir, dirs, files in os.walk(inputDir) :
#    for file in files :
#        count = 0
#        #Get the accession number
#        acc = re.findall(r'^(.*?)_CDS_\d+-\d+_extract\.fna$', file)[0]
#        if acc in accList :
#            #Create new directory
#            outputDir = '%s/%s' % (outputRoot, acc)
#            if not os.path.exists(outputDir) :
#                os.makedirs(outputDir)
#            else :
#                pass
#            fastaFile = os.path.join(subdir, file)
#            #Extract the gene fragment FASTA entries
#            fragments = read_fasta(fastaFile, '|')
#            
#            #Get curvature for sequences
#            for i in range(len(fragments)) :
#                head = fragments[i][0]
#                sequence = fragments[i][1].upper()
#                headL = head[0].split(';')
#                cds = headL[0]
#                accession = re.findall(r'^(\w+_\d+)\.\d$', headL[1])[0]
#                if len(sequence) <= maxLen :
#                    #Nucleotide gene fragment
#                    frag = sequence[:seqLen]
#                    #Calculates the global structure of a DNA molecule from its
#                    #nucleotide sequence according to the dinucleotide wedge model
#                    curveCommand = 'python %s/dnacurve.py "%s" --model="trifonov" --csv="%s/%s_%s.csv" --name="%s"' % (rootDir,frag,outputDir,accession,cds,cds)
#                    os.system(curveCommand)
#                    count += 1
#                else :
#                    print('Sequence is too long.')
#            currentTime = strftime("%Y-%m-%d %H:%M:%S")
#            statusFile.write('Processed %s at %s - %d sequences in %s\n' % (accession, currentTime, count, file))       
#            genCount += 1
#        else :
#            pass

#==============================================================================
#GET CURVATURE DATA FOR A SELECT SET OF ORGANISMS
accList = []
  
with open(accFile, 'r') as acc :
    for a in acc :
        a = a.rstrip()
        #aList = a.split('\t')
        #accession = aList[1]
        accList.append(a)
        
for a in accList[:1] :
    count = 0
    #Get the accession number
    #acc = re.findall(r'^(.*?)_CDS_\d+-\d+_extract\.fna$', file)[0]

    outputDir = '%s/%s' % (outputRoot, a)
    #Create the output folder for the accession if not already created
    if not os.path.exists(outputDir) :
        os.makedirs(outputDir)
    else :
        pass
    fastaFile = '%s/%s_CDS_200-350_extract.fna' % (inputDir, a)
    #Extract the gene fragment FASTA entries
    fragments = read_fasta(fastaFile, '|')
    
    for i in range(len(fragments)) :
        head = fragments[i][0]
        sequence = fragments[i][1].upper()
        #Nucleotide gene fragment
        frag = sequence[:seqLen]
        headL = head[0].split(';')
        cds = headL[0]
        #accession = re.findall(r'^(\w+_\d+)\.\d$', headL[1])[0]
        
        #Calculates the global structure of a DNA molecule from its
        #nucleotide sequence according to the dinucleotide wedge model
        curveCommand = 'python %s/dnacurve.py "%s" --model="trifonov" --csv="%s/%s_%s.csv" --name="%s"' % (rootDir,frag,outputDir,a,cds,cds)
        os.system(curveCommand)
        count += 1
    currentTime = strftime('%d-%m-%Y %H:%M:%S')
    statusFile.write('Processed %s at %s - %d sequences\n' % (a, currentTime, count))
    genCount += 1
        

#Finishing statements
statusFile.write('\nProcessed %d genomes' % genCount)
statusFile.write('\n-----DONE @ %s-----' % strftime('%d-%m-%Y %H:%M:%S'))       
statusFile.close()

#print('\n\n-------DONE-------')