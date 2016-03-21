#!/usr/bin/env python
import re
import os
import sys

'''
Ian Rambo
STEP 1 OF PIPELINE

GOAL: extract DNA coding sequences (CDS) for a genome. Output is a fasta file
with an entry for each CDS.

Extract a full CDS or a fragment of a CDS that includes an upstream promoter.
CDS fragment extends x bases upstream from TSS
(promoter) and x bases downstream (gene).
Calculates the reverse complement for minus strand sequences.

Creation date: November 19, 2015
Last updated: January 20, 2015
'''

'''
HOW TO USE THIS SCRIPT
This script was designed to process genomes downloaded from NCBI. As such, FASTA
and GFF files must have matching identifiers (e.g. NC_021711.fna, NC_021711.gff).

I. Two inputs are required: a directory with genome scaffold FASTA files, and a
directory with GFF files for those genomes.
II. Output files are written to a desired output directory. If the output directory
does not exist, the script will create it.
III. There are two options for this script: extract entire regions of interest
(CDS, gene, mRNA, etc.) or extract gene fragments with upstream and downstream
regions of specified length.

To change which region is extracted, change the value of the variable "region"
within Global Variables.

Usage examples:
To extract a gene fragment with a 100 bp upstream region and 320 bp downstream
region:
python MotifPipeline/gff_extractor.py "/net/biohen29/data/imrambo/BacterialGenomes/" "/net/biohen29/data/imrambo/all.gff/" "/net/biohen29/data/imrambo/00-ExtractData" 100 320

To extract entire CDS sequences:
python MotifPipeline/gff_extractor.py "/net/biohen29/data/imrambo/BacterialGenomes/" "/net/biohen29/data/imrambo/all.gff/" "/net/biohen29/data/imrambo/00-ExtractData"

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

def rev_comp(sequence) :
    '''
    Get the reverse compliment of a DNA sequence.
    '''
    compDict = {'A':'T','T':'A','C':'G','G':'C'}
    compSeq = str()
    for base in sequence :
        if base in compDict :
            compBase = compDict[base]
            compSeq += compBase
        else :
        #Add an unknown base
            compSeq += base
    
    compSeq = compSeq[::-1]
    
    return compSeq

def find_file(filename, path) :
    for root, dirs, files in os.walk(path) :
        if filename in files :
            return os.path.join(root, filename)
            
#==============================================================================
#---------COMMAND LINE ARGUMENTS-----------------
if len(sys.argv) > 1 :
    #Root input directory containing FASTA files
    #FASTA filenames must be in the format <NCBI accession>.fna
    fastaRoot = sys.argv[1]
    #Root input directory containing GFF files
    gffRoot = sys.argv[2]
    #Root output folder
    outDir = sys.argv[3]
    #If the output folder does not exist, create it
    if not os.path.exists(outDir) :
        print('\ncreating output folder %s\n' % outDir)
        os.makedirs(outDir)
    else :
        pass

if len(sys.argv) > 4 :
    #How many bp upstream from TSS
    upstream = int(sys.argv[4])
    #How many bp downstream from TSS
    downstream = int(sys.argv[5])
#==============================================================================
#---------GLOBAL VARIABLES-----------------------
#Count the number of genomes processed
genomes = 0
#Count of valid gene fragments written to output
count = 0
#Total count of CDS in genome
cdsTotal = 0
#Count of CDS sequences shorter than desired fragment length
short = 0
#Count of overlaps
overlap = 0
#Count of fragments excluded due to invalid start codons
invalSC = 0
#Count of sequences with N fill for promoter region
withN = 0
#Which region to extract from - gene, coding sequence (CDS), mRNA, etc.
region = 'CDS'
#==============================================================================
#---------MAIN-----------------------------------
print('\n...extracting %s sequences...\n' % region)

#Specific accession numbers to extract from
targetAccessions = []
with open('/net/biohen29/data/imrambo/targetClassAccessions.txt', 'r') as tca :
    for line in tca :
        line = line.strip()
        lineL = line.split('\t')
        a = lineL[1]
        targetAccessions.append(a)
        
if len(sys.argv) > 4 and all(isinstance(x, int) for x in [upstream, downstream]) :
    print('fragments extend %d bp upstream and %d bp downstream from TSS' % (upstream, downstream))
else :
    print('extracting entire %s sequences\n' % region)

for subdir, dirs, files in os.walk(fastaRoot) :
    for file in files :
        #Extract the accession number
        acc = re.findall(r'^(.*?)\.fna$', file)[0]
        if acc in targetAccessions :
            fastaFile = os.path.join(subdir, file)
            #Get the accession number
            #acc = re.findall(r'(.*?)\.fna', file)[0]
            #Extract the genome scaffold FASTA entry
            scaffold = read_fasta(fastaFile, '|')
            #scaffold FASTA header contents
            HEAD = scaffold[0][0]
            #scaffold sequence
            SEQ = scaffold[0][1].upper()
            gffName = re.sub(r'\.fna$', '.gff', file)
            #Find the respective GFF file
            gffFile = find_file(gffName, gffRoot)
            #Display the accession number of the genome being processed
            print('processing %s\n' % acc)
            #Determine output file format based on extraction type
            if len(sys.argv) > 4 and all(isinstance(x, int) for x in [upstream, downstream]) :
                #Output FASTA file containing extracted CDS or gene fragments
                outFile = '%s_%s_%d-%d_extract.fna' % (acc, region, upstream, downstream)
            else :
                outFile = '%s_%s_extract.fna' % (acc, region)
                
            genomes += 1
            
            with open(gffFile, 'r') as gff, open('%s/%s' % (outDir, outFile), 'w') as OUT :
                #List containing CDS start positions
                cdsStart = []
                #List containing CDS end positions
                cdsEnd = []
                for entry in gff :
                    if not entry.startswith('#') :
                        entry = entry.rstrip()
                        entryL = entry.split('\t')
                        name = entryL[2]
                        #Select for the desired region
                        if name == region :
                            cdsTotal += 1
                            accession = entryL[0]
                            #Length of the intergenic region, if applicable
                            #intergenic = int()
                            start = int(entryL[3]) #CDS start
                            end = int(entryL[4]) #CDS end
                            cdsStart.append(start)
                            cdsEnd.append(end)
                            #Length of CDS
                            cdsLen = end - start
                                      
                            #if len(cdsEnd) > 1 :
                            ##Calculate length of intergenic region
                            #    intergenic = start - cdsEnd[-2]
                            #else :
                            #    pass
                            
                            #DNA strand sense
                            sense = entryL[6]
                            #CDS attributes
                            attr = entryL[8]
                            attributes = attr.split(';')
                            ID = re.sub(r'ID\=', '', attributes[0])  
                            if 'product=' in attr :
                                product = re.findall(r';product\=(.*?);', attr)[0]
                            else :
                                pass
                            if 'protein_id' in attr :
                                proteinID = re.findall(r';protein_id\=(.*?);', attr)[0]
                            else :
                                pass
                            if len(sys.argv) > 4 and all(isinstance(x, int) for x in [upstream, downstream]) :
                                '''
                                Extract a gene fragment containing a hypothetical promoter region
                                with bases extending upstream and downstream from the
                                transcription start site. 
                                '''
                                #Is region long enough to include upstream promoter
                                #and downstream gene fragment (and potential +3 stop codon)?
                                if start >= upstream and cdsLen >= downstream + 3 :
                                    if len(cdsStart) >= 2 : #Extract second CDS or greater
                                        '''
                                        Test that the hypothetical promoter region does
                                        not lie within previous CDS fragment, and that
                                        the end positition of a CDS is not the same as the
                                        start of the subsequent one.
                                        '''
                                        if start - upstream > cdsStart[-2] + downstream and cdsEnd[-2] != cdsStart[-1] :
                                            #Extract upstream promoter region and downstream gene of CDS
                                            sequence = SEQ[(start - upstream):(start + downstream)]
                                                         
                                            if sense == '+' :
                                                #Find the start codon
                                                startCodon = sequence[:3]
                                                #Test that the start codon is valid, and is not a stop codon
                                                if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                                    OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                                    sense,proteinID,product,sequence))
                                                    count += 1
                                                
                                                else :
                                                    invalSC += 1
                                                        
                                            elif sense == '-' :
                                                #If minus strand, give the reverse compliment
                                                rcSequence = rev_comp(sequence)
                                                #Find the start codon
                                                startCodon = rcSequence[:3]
                                                #Make sure the start codon is valid, and is not a stop codon
                                                if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                                    OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                                    sense,proteinID,product,rcSequence))
                                                    count += 1
                                                else :
                                                    invalSC += 1
                                                    
                                            else :
                                                print('Strand sense not found for %s' % ID)
                                        else :
                                            overlap += 1
                                    else :
                                        #Extract the first CDS
                                        if start >= upstream :
                                            sequence = SEQ[(start - upstream):(start + downstream)]
                                            if sense == '+' :
                                                #Find the start codon
                                                startCodon = sequence[:3]
                                                #Test that the start codon is valid, and is not a stop codon
                                                if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                                    OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                                    sense,proteinID,product,sequence))
                                                    count += 1
                                                
                                                else :
                                                    invalSC += 1
                                                
                                            elif sense == '-' :
                                                #If minus strand, give the reverse compliment
                                                rcSequence = rev_comp(sequence)
                                                #Find the start codon
                                                startCodon = rcSequence[:3]
                                                #Make sure the start codon is valid, and is not a stop codon
                                                if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                                    OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                                    sense,proteinID,product,rcSequence))
                                                    count += 1
                                                
                                                else :
                                                    invalSC += 1
                                            else :
                                                  print('Strand sense not found for %s' % ID)
                                        else :
                                        #If the first CDS position is less than the upstream length,
                                        #fill the upstream positions with N
                                            nfill = 'N' * (upstream - start)
                                            sequence = nfill + SEQ[:(start + downstream)]
                                            withN += 1
                                        
                                            if sense == '+' :
                                                #Find the start codon
                                                startCodon = sequence[:3]
                                                #Test that the start codon is valid, and is not a stop codon
                                                if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                                    OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                                    sense,proteinID,product,sequence))
                                                    count += 1
                                                
                                                else :
                                                    invalSC += 1
                                                
                                            elif sense == '-' :
                                                #If minus strand, give the reverse compliment
                                                rcSequence = rev_comp(sequence)
                                                #Find the start codon
                                                startCodon = rcSequence[:3]
                                                #Make sure the start codon is valid, and is not a stop codon
                                                if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                                    OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                                    sense,proteinID,product,rcSequence))
                                                    count += 1
                                                else :
                                                    invalSC += 1
                                            else :
                                                  print('Strand sense not found for %s' % ID)
                                else :
                                    short += 1
                            else :
                                '''
                                Extract the entire sequence.
                                Determine the reverse compliment if applicable.
                                '''
                                sequence = SEQ[start:end]                      
                                if sense == '+' :
                                    #Find the start codon
                                    startCodon = sequence[:3]
                                    #Test that the start codon is valid, and is not a stop codon
                                    if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                        OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                        sense,proteinID,product,sequence))
                                        count += 1
                                    
                                    else :
                                        invalSC += 1
                                
                                elif sense == '-' :
                                    #If minus strand, give the reverse compliment
                                    rcSequence = rev_comp(sequence)
                                    #Find the start codon
                                    startCodon = rcSequence[:3]
                                    #Make sure the start codon is valid, and is not a stop codon
                                    if re.match(r'[ACTG]{3}', startCodon) and startCodon not in ['TAA','TAG','TGA'] :
                                        OUT.write('>%s;%s;%d;%d;%s;%s;%s\n%s\n' % (ID,accession,start,end,\
                                        sense,proteinID,product,rcSequence))
                                        count += 1
                                    else :
                                        invalSC += 1
                                        
                                else :
                                    print('Strand sense not found for %s' % ID)
                                                    
                                            
                                            
                                            
        else :
            pass
                      
            
##Print the output stats
print('\nProcessed %d genomes\n' % genomes)
print('\n%d %s sequences of %d total written\n' % (count,region,cdsTotal))
print('%d %s sequence(s) shorter than required length\n' % (short,region))
print('%d fragment(s) excluded due to overlap with an existing fragment\n' % overlap)
print("%d fragment(s) with N's as placeholders in promoter\n" % withN)
print('%d sequence(s) excluded due to invalid start codons\n' % invalSC)

print('-----DONE-----\n')

#####    #####
#-----END-----
#####    #####