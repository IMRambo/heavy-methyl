#!/usr/bin/env python
'''
Ian Rambo
Created: February 4, 2016
Last modified: February 4, 2016

GOAL: Create new methylation score files for methylation state comparison.
'''

#==============================================================================
rootDir = "/Users/imrambo/R/ORMP/data"
metFile1 = "%s/S06-08-QuantMetMeasureData_copy.txt" % rootDir
metFile2 = "%s/S15-08-QuantMetMeasureData_copy.txt" % rootDir
outFile = "%s/MetCCGG-0615.txt" % rootDir
#==============================================================================
#outFile.write()

#Dictionary containing information for the first methylation score file
met1Dict = dict()
with open(metFile1, 'r') as met1 :
    met1.readline()
    for m in met1 :
        m = m.rstrip()
        mVars = m.split('\t')
        contig = mVars[0]
        metPos = mVars[33]
        #Header containing taxonomic and protein annotation
        taxGenHeader = '|'.join(mVars[1:33])
        #Methylated score metric
        GPmet = float(mVars[34])
        #Non-methylated score metric
        GPumt = float(mVars[35])
        #Methylation score
        GPscore = float(mVars[36])
        localD = float(mVars[37])
        solo = int(mVars[38])
        metDesc = str()
        if GPmet < 5.65 and GPumt < 2.00 :
            metDesc = 'UNK'
        elif GPmet > 5.65 and GPumt < 2.00 :
            metDesc = 'MET'
        elif GPmet < 5.65 and GPumt > 2.00 :
            metDesc = 'UMT'
        elif GPmet > 5.65 and GPumt > 2.00 :
            metDesc = 'MIX'
        else :
            pass
        
        met1Dict[(contig, metPos)] = dict([('contig', contig), ('header', taxGenHeader),\
            ('position', metPos),('GPmet',str(GPmet)),('GPumt', str(GPumt)),('GPscore', str(GPscore)),\
        ('localD', str(localD)),('solo', solo),('metDesc', str(metDesc))])

    

with open(metFile2, 'r') as met2 :
    met2.readline()
    for m in met2 :
        m = m.rstrip()
        mVars = m.split('\t')
        contig = mVars[0]
        metPos = mVars[33]
        #Header containing taxonomic and protein annotation
        taxGenHeader = '|'.join(mVars[1:33])
        GPmet = float(mVars[34])
        GPumt = float(mVars[35])
        GPscore = float(mVars[36])
        localD = float(mVars[37])
        solo = mVars[38]
        metDesc = str()
        if GPmet < 5.65 and GPumt < 2.00 :
            metDesc = 'UNK'
        elif GPmet >= 5.65 and GPumt < 2.00 :
            metDesc = 'MET'
        elif GPmet < 5.65 and GPumt > 2.00 :
            metDesc = 'UMT'
        elif GPmet >= 5.65 and GPumt > 2.00 :
            metDesc = 'MIX'
        else :
            pass
        
        if (contig, metPos) in met1Dict :
            shift = '%s:%s' % (met1Dict[(contig, metPos)]['metDesc'], metDesc)
            print(shift)
            
        
            
            

        
    
    
    
