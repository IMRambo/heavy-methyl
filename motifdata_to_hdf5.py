
'''
Ian Rambo
Created: January 27, 2016
Last updated: January 29, 2016
GOAL: Create an HDF5 file containing sequence motif frequency and curvature data.

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
rootDir = '/Users/imrambo/R/ORMP/data'
outDir = '%s/02-MotifDataHDF5' % rootDir
motifDir = '%s/01-Test_ClusterData' % rootDir

#Input files
lineageFile = '%s/taxonomic_lineages_genomes.txt' % rootDir
targetAccessionFile = '%s/targetClassAccessions.txt' % rootDir

#Output files
hdf5_out = '%s/test.hdf5' % outDir
progressFile = '%s/progressFile-00.txt' % outDir  

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
    print('\n%s already exists, doing nothing...\n' % outDir)
    pass

progressFile = open(progressFile, 'w')
#Obtain accession numbers
targetAccessions = []
with open(targetAccessionFile, 'r') as ta :
    for line in ta :
        line = line.rstrip()
        lineL = line.split('\t')
        accession = lineL[1]
        targetAccessions.append(accession)
        
#Create a compound data type
dt = np.dtype([('ntPosition', np.int32),('ATmotif', np.float64),('CpGmotif', np.float64),\
            ('CpGnorm', np.float64),('CpGoe', np.float64),\
            ('GATCexpected', np.float64),('GATCmotif', np.float64),\
            ('GCcontent', np.float64),('GCmotif', np.float64),\
            ('GpCnorm', np.float64)])

print('...generating datasets...\n')
with h5py.File(hdf5_out, 'w') as h5file, open(lineageFile, 'r') as genpath :
    '''
    Create a hierarchy based on taxonomic rank. A dataset is created for
    each cds per accession. 
    '''
    for path in genpath :
        #Taxonomic lineage - Phylum to Accession (stand-in for species)
        path = path.rstrip()
        pathL = path.split('/')
        accession = pathL[-1]
        if accession in targetAccessions :
            #hierarchical taxonomic lineage for accession
            groupName = '/' + path
            
            #Directory containing motif datasets
            mfileDir = '%s/%s_CDS_200-350_extract.fna' % (motifDir, accession)
            
            for subdir, dirs, files in os.walk(mfileDir) :
                #Create a nested group for the accession
                grp = h5file.create_group(groupName)
                '''
                Obtain motif frequency data for each CDS of each accession
                and add to its respective HDF5 dataset.
                '''
                for mfile in files :
                    mfilePath = os.path.join(subdir, mfile)
                    if accession in mfilePath :
                        #print('processing motif data for %s...' % accession)
                        #Which motif is being looked at?
                        metric = mfile.replace('.txt', '')
                        with open(mfilePath, 'r') as mf :
                            #This is where you made the group before
                            for line in mf :
                                nGenes += 1
                                line = line.rstrip()
                                lineL = line.split('\t')
                                seqHead = lineL[0].split(';')
                                #cds number as string
                                cdsStr = seqHead[0].replace('>', '')
                                #cds number as integer
                                cdsInt = int(seqHead[0].replace('>cds', ''))
                                #flag if the base pair was 'N'
                                flag = '-'
                                #motif frequency values
                                values = lineL[1:]
                                #convert values to a NumPy array of floats and replace flag values
                                #with a very small number
                                floatValues = np.array([float(1e-63) if v == flag else float(v) for v in values], dtype = np.float64)
                                #nucleotide positions
                                ntPos = np.arange(1, 501, dtype = np.int32)
                                #Is a dataset for this CDS already in the group? If so, add the values for
                                #this score metric
                                if cdsStr in grp :
                                    dPath = '%s/%s' % (groupName, cdsStr)
                                    dset = grp[dPath]
                                    dset[..., metric] = floatValues
                                #If the dataset has not been made yet, create it and add
                                #nucleotide positions and values
                                else :
                                    newDset = grp.create_dataset(cdsStr, (500, ), dtype = dt,\
                                                                 chunks = (125,), compression='gzip',\
                                    compression_opts=6)
                                    newDset[..., 'ntPosition'] = ntPos
                                    newDset[..., '%s' % metric] = floatValues
                                  
                          
                    else :
                        print('Directory %s not found' % mfilePath)
                        
            nGenomes += 1
            #Time the genome processing was completed
            currentTime = strftime("%Y-%m-%d %H:%M:%S")
            progressFile.write('processed %s at %s ... genome %d of %d\n' % (accession,currentTime,\
                                                                      nGenomes,len(targetAccessions)))
    #--------------------------------------------------------------------------
            #'''
            #Obtain curvature data for each accession and add to respective HDF5 dataset.
            #'''
            #for subdir, dirs, files in os.walk('%s/%s' % (curveDir, accession)) :
            #    for curveCSV in files :
            #        #accession = re.findall(r'^(\w{2}_\d+)_cds\d+\.csv$', curveCSV)[0]
            #        curveCSVPath = os.path.join(subdir, curveCSV)
            #        cds = re.findall(r'_(cds\d+)\.csv', curveCSV)[0]
            #        
            #        #print('processing curvature data for %s...' % accession) 
            #        with open(curveCSVPath, 'r') as c :
            #            c.readline()
            #            curveData = []
            #            readCSV = csv.reader(c)
            #            for row in readCSV :
            #                #Column with curvature data
            #                cd = float(row[2])
            #                curveData.append(cd)
            #            if len(curveData) == seqlen :
            #                cdArray = np.array(curveData)
            #                ntPos = np.arange(1, 501, dtype = np.int32)
            #                if cds in grp :
            #                    dPath = '%s/%s' % (groupName, cds)
            #                    dset = grp[dPath]
            #                    dset[..., 'curvature'] = cdArray
            #                #If the dataset has not been made yet, create it and add
            #                #nucleotide positions and values
            #                else :
            #                    newDset = grp.create_dataset(cds, (500, ), dtype = dt, chunks = (125,))
            #                    newDset[..., 'ntPosition'] = ntPos
            #                    newDset[..., 'curvature'] = cdArray
            #                    
            #                
            #            else :
            #                print('Curvature length not the same as sequence length')   
                            

progressFile.write('\nprocessed %d genes from %d genomes\n' % (nGenes, nGenomes))                                          
progressFile.write('-----DONE-----\n')
progressFile.close()
print('-----DONE-----\n')
#-----END-----
