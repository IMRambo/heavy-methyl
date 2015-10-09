#!/usr/bin/env python
#Ian Rambo
#Oyster Rocks Methylation Project
#University of Delaware
#Last updated: October 5, 2015
#==============================================================================
#COG header annotation
'''
GOAL: add RAMMCAP COG results to FASTA headers.
'''
#==============================================================================
import os
import re
#==============================================================================
def function_codes(filename) :
    '''
    Obtain function descriptions for COG function codes.
    '''
    fCodes = dict()
    with open(filename, 'r') as func :
        for f in func :
            f = f.rstrip()
            if not f.startswith('#') :
                fL = f.split('\t')
                fCodes[fL[0]] = fL[1]
        return fCodes  
    func.close()
#------------------------------------------------------------------------------
def cognames(namefile, funcfile) :
    '''
    Create dictionary of COG names and functions.
    Use COG name file and COG function key file.
    '''
    COG = dict()
    fCodes = function_codes(funcfile)
    
    with open(namefile, 'r') as names :
        for record in names :
            record = record.rstrip()
            if not record.startswith('#') :
                recordL = record.split('\t')
                
                if len(recordL[1]) > 1 :
                    function = []
                    for code in recordL[1] :
                        if code in fCodes :
                            function.append(fCodes[code])
                            COG[recordL[0]] = dict([('fcode', recordL[1]),
                                ('function', ','.join(function)),
                                ('name', recordL[2])])
                
                elif len(recordL[1]) == 1 :
                    COG[recordL[0]] = dict([('fcode', recordL[1]),
                        ('function', fCodes[recordL[1]]),
                        ('name', recordL[2])])
                else :
                    pass
            else :
                pass
        return COG
    names.close()
#------------------------------------------------------------------------------
def parse_cog_results(cogresultsfile):
    '''
    Create a dictionary of RAMMCAP COG results containing queries, COG names,
    function codes, enzyme names, function descriptions, evalues, and scores.
    '''
    output_dict = dict()
    COG = cognames('./proteinAnnotation/cognames2003-2014.txt',
               './proteinAnnotation/cogfunction2003-2014.txt')
    with open(cogresultsfile, 'r') as cog :
        #Skip the header line
        cog.readline()
        #next(cog)
        for result in cog :
            result = result.rstrip()
            #Remove open reading frame number
            result = re.sub('\.\d+\tCOG', '\tCOG', result)
            resultL = result.split('\t')
            if resultL[0] in output_dict :
                #If contig is already seen in results dictionary, choose result with
                #the lowest e-value and highest score
                if float(resultL[2]) < float(output_dict[resultL[0]]['evalue']) and\
                float(resultL[4]) > float(output_dict[resultL[0]]['score']) and\
                resultL[1] in COG :
                    output_dict[resultL[0]] = dict([('query',resultL[0]),('hit',resultL[1]),
                        ('fcode',COG[resultL[1]]['fcode']),('name',COG[resultL[1]]['name']),
                        ('function',COG[resultL[1]]['function']),('evalue',str(resultL[2])),
                        ('score',str(resultL[4]))])
                    
            #If contig seen for the first time, add to output dictionary
            else :
                if resultL[1] in COG :
                    output_dict[resultL[0]] = dict([('query',resultL[0]), ('hit',resultL[1]),
                        ('fcode',COG[resultL[1]]['fcode']), ('name',COG[resultL[1]]['name']),
                        ('function',COG[resultL[1]]['function']),('evalue',str(resultL[2])),
                        ('score',str(resultL[4]))])
                    
            
    cog.close()
    return output_dict
#------------------------------------------------------------------------------
def read_FASTA(filename) :
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
                header = line[1:].split('|')
                seq = ''
            else :
                seq += line
            line = fasta.readline()[:-1]
        sequences.append((header, seq))
    fasta.close()
    return(sequences)
#------------------------------------------------------------------------------
def fasta_write_test(fasta) :
    '''
    Test if FASTA output file wrote correctly. 
    '''
    with open(fasta, 'r') as fa :
        lines = fa.readlines()
        if lines[0].startswith('>') and lines[0].endswith('\n') and not '>' in lines[-1] :
            print('FASTA file %s written successfully.' % fasta)
        else :
            print('FASTA file %s not written successfully.' % fasta)
    fa.close()
    return
#==============================================================================
os.chdir('/Users/imrambo/Documents/MS_thesis/')

#Read PhymmBL/Kraken/FOAM annotated FASTA
sequences = read_FASTA('S24-27cm_pkhAnnotated_01.fa')

#RAMMCAP COG results test file - header and first ten results
#testCOG = './proteinAnnotation/S03-06cm_COG/S03-06cm_COG_parse_output.1'
#testCOG = './proteinAnnotation/S12-15cm_COG/S12-15cm_COG_parse_output.1'
testCOG = './proteinAnnotation/S24-27cm_COG/S24-27cm_COG_parse_output.1'

cogResDict = parse_cog_results(testCOG)

#Output FASTA with COG annotation added to header
outputHandle = 'S24-27cm_pkhcAnnotated_01.fa'
outfile = open(outputHandle, 'w')

#Append COG results to sequence headers
for i in range(0, len(sequences)) :
    header = sequences[i][0]
    seqid = sequences[i][0][0]
    sequence = sequences[i][1]
    if seqid in cogResDict :
        annotation = [cogResDict[seqid]['hit'], cogResDict[seqid]['evalue'],
                      cogResDict[seqid]['score'], cogResDict[seqid]['fcode'],
                      cogResDict[seqid]['name'], cogResDict[seqid]['function']]
        header.extend(annotation)
        outfile.write('>' + '|'.join(header) + '\n' + sequence + '\n')
    
    else :
        annotation = ['NA']*6
        header.extend(annotation)
        outfile.write('>' + '|'.join(header) + '\n' + sequence + '\n')

outfile.close()

fasta_write_test(outputHandle)