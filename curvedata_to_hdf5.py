#!/usr/bin/env python

'''
Ian Rambo
Created: February 11, 2016
Last modified: February 29, 2016
GOAL: Create an HDF5 file containing curvature data.

'''
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
import numpy as np
import os
import re
import h5py
import csv
import datetime
from time import strftime
#==============================================================================
#-----GLOBAL VARIABLES-----
#I/O directories
rootDir = '/net/biohen29/data/imrambo'
outDir = '%s/02-CurveDataHDF5' % rootDir
#Input directory with subdirectories for each accession containing csv files
curveDir = '%s/01-targetCurveData' % rootDir

#Input files
lineageFile = '%s/taxonomic_lineages_genomes.txt' % rootDir
targetAccessionFile = '%s/targetClassAccessions.txt' % rootDir

#Output files
hdf5_out = '%s/targetCurveData.hdf5' % outDir
progressFile = '%s/curveHDF5ProgressFile-%s.txt' % (outDir, strftime('%d-%m-%Y'))

#Length of sequence fragments
seqlen = 500
#Genome count
nGenomes = 0
#Gene count
nGenes = 0
#==============================================================================
#-----FUNCTIONS-----
def printname(name) :
    '''
    Print the names in an HDF5 file.
    '''
    print name
#==============================================================================
#-----MAIN CODE-----
print('-----BEGIN-----')

if not os.path.exists(outDir) :
    print('\ncreating output folder %s\n' % outDir)
    os.makedirs(outDir)
else :
    print('\nDirectory %s already exists, doing nothing...\n' % outDir)
    pass


#The progress file will store the date and time when each entry was completed.
progressFile = open(progressFile, 'w')

#Obtain accession numbers
#targetAccessions = []
#with open(targetAccessionFile, 'r') as ta :
#    for line in ta :
#        line = line.rstrip()
#        lineL = line.split('\t')
#        accession = lineL[1]
#        targetAccessions.append(accession)

linPathList = []
with open(lineageFile, 'r') as genLin :
    for lin in genLin :
        lin = lin.rstrip()
        linPath = '/' + lin
        linPathList.append(linPath)
    

#Create a compound data type
#dt = np.dtype([('ntPosition', np.int32),('curvature', np.float64),('bendAngle', np.float64),\
            #('curvatureAngle', np.float64)])
#dt = np.dtype([('ntPosition', np.int32),('curvature', np.float64)])
dt = np.dtype([('curvature', np.float64)])

print('...generating datasets...\n')
with h5py.File(hdf5_out, 'w') as h5file :
    '''
    Create a hierarchy based on taxonomic rank. A dataset is created for
    each cds per accession. 
    '''
    for subdir, dirs, files in os.walk(curveDir) :
        if re.search(r'NC_\d+', subdir) :
            accession = subdir.split('/')[6]
            path = [l for l in linPathList if accession in l][0]

            grp = h5file.create_group(path)
            
        
            '''
            Obtain curvature data for each accession and add to respective HDF5 dataset.
            '''
            for curveCSV in files :
                curveCSVPath = os.path.join(subdir, curveCSV)
            
                cds = re.findall(r'NC_\d+_(cds\d+)\.csv', curveCSVPath.split('/')[-1])[0]
                
                with open(curveCSVPath, 'r') as c :
                    progressFile.write('Processing %s of %s\n' % (cds, accession))
                    c.readline()
                    #ntPosData = []
                    curvatureData = []
                    #bendAngleData = []
                    #curveAngleData = []
                    readCSV = csv.reader(c, delimiter = ",")  
                    for row in readCSV :
                        #Column with nucleotide position
                        #ntPos = int(row[1])
                        #Column with curvature data
                        curvature = float(row[2])
                        #Column with bend angle
                        #bendAngle = float(row[3])
                        #Column with curvature angle
                        #curveAngle = float(row[4])
                        #ntPosData.append(ntPos)
                        curvatureData.append(curvature)
                        #bendAngleData.append(bendAngle)
                        #curveAngleData.append(curveAngle)
                    #If the csv file has less bases than desired, fill 
                    #if len(ntPosData) < seqlen :
                    #    remainingPos = list(range(len(ntPosData), seqlen+1))
                    #    ntPosData.extend(remainingPos)
                    #else :
                    #    pass
                    if len(curvatureData) < seqlen :
                        lenDiff = seqlen - len(curvatureData)
                        remainingCurve = [0] * lenDiff
                        curvatureData.extend(remainingCurve)
                    else :
                        pass
                    #Convert data lists to numpy arrays
                    #ntPosArray = np.array(ntPosData, dtype = np.int32)
                    curvatureArray = np.array(curvatureData, dtype = np.float64)
                    #bendAngleArray = np.array(curvatureData, dtype = np.float64)
                    #curveAngleArray = np.array(curvatureData, dtype = np.float64)
                  
                    dPath = '%s/%s' % (path, cds)
                    grp.create_dataset(cds, (seqlen, ), dtype = dt,\
                                                             chunks = (125,), compression='gzip',\
                                compression_opts=6)
                    dset = grp[dPath]
                    #dset[..., 'ntPosition'] = ntPosArray
                    dset[..., 'curvature'] = curvatureArray
                    #dset[..., 'bendAngle'] = bendAngleArray
                    #dset[..., 'curvatureAngle'] = curveAngleArray

                        #If the dataset has not been made yet, create it and add
                        #nucleotide positions and values
                    #else :
                    #    newDset = grp.create_dataset(cds, (500, ), dtype = dt, chunks = (125,))
                    #    newDset[..., 'ntPosition'] = ntPos
                nGenes += 1       
                        
            nGenomes += 1
            #Time the genome processing was completed
            currentTime = strftime("%Y-%m-%d %H:%M:%S")
            progressFile.write('processed %s at %s - genome %d\n' % (accession,currentTime,nGenomes))
#--------------
#                            
#
progressFile.write('\nprocessed %d genes from %d genomes\n' % (nGenes, nGenomes))                                          
progressFile.write('-----DONE @ %s-----\n' % strftime("%Y-%m-%d %H:%M:%S"))
progressFile.close()
print('-----DONE-----\n')
#-----END-----
