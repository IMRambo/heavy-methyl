#!/usr/bin/env python
'''
Ian Rambo
Created: January 24, 2016
Last modified: January 25, 2016

GOAL: transpose motif cluster data for analysis in R.
'''
import os
import re
#==============================================================================
#-----GLOBAL VARIABLES-----
rootDir = '/net/biohen29/data/imrambo'
inputDir = '%s/01-ClusterData' % rootDir
outputDir = '%s/02-CompiledClusterData' % rootDir
geneCount = 0

if not os.path.exists(outputDir) :
    print('\ncreating output folder %s\n' % outputDir)
    os.makedirs(outputDir)
else :
    print('\n%s already exists\n!' % outputDir)
    pass

#-----OUTPUT FILES-----
outATmotif = open('%s/ATmotif_compiled.txt' % outputDir, 'w')
outCpGnorm = open('%s/CpGnorm_compiled.txt' % outputDir, 'w')
outGATCmotif = open('%s/GATCmotif_compiled.txt' % outputDir, 'w')
outGCmotif = open('%s/GCmotif_compiled.txt' % outputDir, 'w')
outCpGOE = open('%s/CpG-OE_compiled.txt' % outputDir, 'w')
outCpGmotif = open('%s/CpGmotif_compiled.txt' % outputDir, 'w')
outGATCexpect = open('%s/GATCexpected_compiled.txt' % outputDir, 'w')
outGCcontent = open('%s/GCcontent_compiled.txt' % outputDir, 'w')
outGpCnorm = open('%s/GpCnorm_compiled.txt' % outputDir, 'w')
#==============================================================================
print('\n-----BEGIN-----\n')

outfiles = [outATmotif,outCpGnorm,outGATCmotif,outGCmotif,outCpGOE,outCpGmotif,outGATCexpect,outGCcontent,outGpCnorm]
for f in outfiles :
    f.write('cds\taccession\tntPosition\tscore\n')

for subdir, dirs, files in os.walk(inputDir) :
    for file in files :
        scoreFile = os.path.join(subdir, file)
        print('transposing %s ...\n' % scoreFile)
        metric = file.replace('.txt', '')
        with open(scoreFile, 'r') as sf :
            for line in sf :
                line = line.rstrip()
                lineL = line.split('\t')
                header = lineL[0].split(';')
                cds = header[0].replace('>', '')
                accession = re.findall(r'^(.*?)\.\d+$', header[1])[0]
                scores = lineL[1:]
                scores = ['NA' if i == '-' else i for i in scores]
                if metric == "ATcontent" :       
                    for i in range(len(scores)) :
                        outATmotif.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))
                elif metric == "CpGnorm" :
                    for i in range(len(scores)) :
                        outCpGnorm.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))
                elif metric == "GATCmotif" :
                    for i in range(len(scores)) :
                        outGATCmotif.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))
                elif metric == "GCmotif" :
                    for i in range(len(scores)) :
                        outGCmotif.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))
                elif metric == "observed-expected" :
                    for i in range(len(scores)) :
                        outCpGOE.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))
                elif metric == "CpGmotif" :
                    for i in range(len(scores)) :
                        outCpGmotif.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))    
                elif metric == "GATCexpected" :
                    for i in range(len(scores)) :
                        outGATCexpect.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))    
                elif metric == "GCcontent" :
                    for i in range(len(scores)) :
                        outGCcontent.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))    
                elif metric == "GpCnorm" :
                    for i in range(len(scores)) :
                        outGpCnorm.write('%s\t%s\t%d\t%s\n' % (cds,accession,i+1,scores[i]))    
                else :
                    print('Cannot find score metric for ' % file)
                
                geneCount += 1
                
for f in outfiles :
    f.close()

print('transposed %d genes\n' % geneCount)
print('\n-----DONE-----\n')