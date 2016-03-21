#!/usr/bin/env python
rootdir = '/net/biohen29/data/imrambo'

completed = []
outFile = open('%s/remainingCurveAccessions.txt' % rootdir, 'w')
with open('%s/completedCurveAccessions.txt' % rootdir, 'r') as cca :
    for c in cca :
        c = c.rstrip()
        completed.append(c)
        
        
with open('%s/targetClassAccessions.txt' % rootdir, 'r') as ta :
    for t in ta :
        t = t.rstrip()
        tList = t.split('\t')
        acc = tList[1]
        if acc in completed :
            pass
        else :
            outFile.write('%s\n' % acc)
            

outFile.close()