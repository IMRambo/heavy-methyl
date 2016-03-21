#!/usr/bin/env python
'''
Create an HDF5 file containing curvature and motif data.
'''
import numpy as np
import os
import re
import h5py

with h5py.File('testfile.hdf5', 'w') as f :
    bigData = np.ones((100,1000))
    f["my dataset"] = bigData
    dset = f["my dataset"]
    #print(dset)
    #print(dset.shape)
    out = dset[...]
    #print(out)
    
    dset[1:4,1] = 2.0
   
    big_out = np.empty((100,1000), dtype = np.float64)
    dset.read_direct(big_out)
    #print(dset[...])
    dset = f.create_dataset('filled', (2,2), dtype = np.int32, fillvalue = 42)
    #print(dset[...])
    
    dset = f.create_dataset('resizable', (2,2), maxshape = (2,2))
    #print(dset.shape)
    dset = f.create_dataset('unlimited', (2,2), maxshape = (2,None))
    
    dset = f.create_dataset('chunked', (100,480,640),
                            dtype = 'i1', chunks = (1,64,64))
    #Check the chunk shape
    #print(dset.chunks)


with h5py.File('groups.hdf5', 'w') as g :
    tPhylum = g.create_group('phylum')
    #print(tPhylum)
    #print(tPhylum.name)
    tClass = tPhylum.create_group('class')
    #print(tClass.name)
    tGenus = g.create_group('/phylum/class/order/family/genus')
    #print(tGenus.name)
    
    g['dataset1'] = 1.0
    g['dataset2'] = 2.0
    g['dataset3'] = 3.0
    tClass['dataset4'] = 4.0
    tClass['dataset5'] = 5.0
    
    dset1 = g['dataset1']
    dset4 = g['/phylum/class/dataset4']
    #print(dset1[...])
    #print(dset4[...])
    #print(len(g))
    #print(len(g['/phylum/class']))
    
    groups = [x for x in g['/phylum/class']]
    print(groups)
    
    groups = [(x, y) for x, y in g.iteritems()]
    print(groups)
    

    if 'dataset4' in g['/phylum/class'] :
        print 'yup'
    else :
        print 'nope'
        
    #visitor iteration
gname = 'subgroup_1'
f = h5py.File('visit_test.hdf5', 'w')
f.create_dataset('top_dataset', data = 1.0)
f.create_group( 'top_group_1' )
f.create_group( 'top_group_1/%s' % gname )
f.create_dataset('top_group_1/subgroup_1/sub_dataset_1', data = 1.0)

f.create_group( 'top_group_2' )
f.create_dataset('top_group_2/sub_dataset_2', data = 1.0)

def printname(name) :
    print name
    
#f.visit(printname)


