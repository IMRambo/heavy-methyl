#!/Users/imrambo/Documents/BINF868/env/bin python
#PacBio SMRT: E. coli 20kb Size Selected Library with P6 C4 chemistry
#Ian Rambo
#BINF 868
#Last updated: September 26, 2015
#==============================================================================
'''
GOAL: learn how to work with PacBio SMRT data in HDF5 format. 
'''
import h5py
import numpy as np
import os
import pbcore
from pbcore.io import BasH5Reader
#==============================================================================
def printname(name) :
    '''
    A simple print function for use with iterating over an entire HDF5 file
    via the Group methods visit() and visititems().
    '''
    print(name)
#==============================================================================
os.chdir('/Users/imrambo/Documents/BINF868/EcoliK12_SMRT')
filepath = './Analysis_Results/m140930_121059_sherri_c100688052550000001823139503241542_s1_p0.1.bax.h5'
'''
Ecoli_hdf5_1 = h5py.File(filepath, 'r')

Ecoli_hdf5_1.name

#Ecoli_hdf5_1.visit(printname)
keys = Ecoli_hdf5_1.keys()
print(keys)

basecallZMW = Ecoli_hdf5_1['PulseData/BaseCalls/ZMWMetrics/']
print(basecallZMW.name)
print(basecallZMW['BaseRate'].dtype)
print(basecallZMW['BaseRate'].shape)

baseRate = basecallZMW['BaseRate']
readScore = basecallZMW['ReadScore']
print(baseRate[460:500])
print(readScore[0:30])
print(np.mean(readScore))

Ecoli_hdf5_1.close()
'''

#basH5 file
filepath = './Analysis_Results/m140930_121059_sherri_c100688052550000001823139503241542_s1_p0.bas.h5'
Ecoli_bas = h5py.File(filepath, 'r')
Ecoli_bas.visit(printname)


ZMW = Ecoli_bas['PulseData/BaseCalls']
print(ZMW.name)
print(ZMW['HoleStatus'].dtype)
print(ZMW['HoleStatus'].shape)
ecbas = pbcore.io.BasH5Reader(Ecoli_bas)

Ecoli_bas.close()

    
