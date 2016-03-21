#!/usr/bin/env python
#Ian Rambo
'''
Created: January 30, 2016
Last modified: January 31, 2016

GOAL: parse HMMSEARCH output to deliver the best hit. This script will
select hits based on the best e-value. 
'''
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
import os
import re
#==============================================================================
#-----GLOBAL VARIABLES-----
rootDir = '/Users/imrambo/R/ORMP/data'
inputDir = '%s/01-targetGenomes_PfamHMM_select' % rootDir
outDir = '%s/02-targetGenomes_PfamHMM_bestHit' % rootDir
fileCount = 0
#==============================================================================
#-----FUNCTIONS-----
def hmm_besthit(hmmFile) :
    '''
    Loop through a modified HMMSEARCH results file and select the best hit
    for each cds per accession.
    '''
    hmmDict = dict()
    with open(hmmFile, 'r') as hmm :
        for h in hmm :
            if h.startswith('cds') :
                h = h.rstrip()
                hFields = h.split()
                seqID = hFields[0]
                pfam = hFields[1]
                seqIDFields = seqID.split(';')
                cds = seqIDFields[0]
                accession = re.findall(r'((\w+_\d+)\.\d+)', seqIDFields[1])[0][1]
                evalue = float(hFields[2])
                #If the accession and CDS exist already, choose the hit
                #with the best e-value
                if cds in hmmDict :
                    if float(evalue) < float(hmmDict[cds]['evalue']) :
                        hmmDict[cds] = dict([('accession',accession),('cds',cds),\
                            ('evalue', str(evalue)), ('pfamily', pfam)])
                        
                
                else :
                    hmmDict[cds] = dict([('accession',accession),('cds',cds),\
                        ('evalue', str(evalue)), ('pfamily', pfam)])
        return(hmmDict)            
#==============================================================================
#-----MAIN-----
print('-----BEGIN-----')
print('...finding best hits...')

if not os.path.exists(outDir) :
    print('\ncreating output folder %s\n' % outDir)
    os.makedirs(outDir)
else :
    print('\n%s already exists\n' % outDir)
    pass

#Get the best hits for each accession
for subdir, dirs, files in os.walk(inputDir) :
    for hmmFile in files :
        hmmPath = os.path.join(subdir, hmmFile)
        hmmDict = hmm_besthit(hmmPath)
        fileHandle = hmmFile.replace('hmmsearch.txt', 'bestHit')
        #Open the output file to contain best hit results
        with open('%s/%s.txt' % (outDir, fileHandle), 'w') as bh :   
            for c in hmmDict :
                bh.write('%s\t%s\t%s\n' % (hmmDict[c]['accession'], hmmDict[c]['cds'], hmmDict[c]['pfamily']))
    
        fileCount += 1   
            
        
print('Created %d new files' % fileCount)
print('-----DONE-----')  