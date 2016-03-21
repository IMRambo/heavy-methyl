#!/usr/bin/env python
'''
Parse transposase Pfam accessions and descriptions from pfam.xfam.org search.

Created: March 6, 2016
Last modified: March 6, 2016
Ian Rambo
'''
import re
import urllib.request
import urllib.parse
#==============================================================================
#GLOBAL VARIABLES
rootDir = '/Users/imrambo/Dropbox/OysterRocksMethylationProject'

print('\n-----BEGIN-----\n')
query = 'transposase'

urlString = 'http://pfam.xfam.org/search/keyword?query=%s&submit=Submit' % query

url = urllib.request.urlopen(urlString)
urlSource = str(url.read())
pfam = re.findall(r'family/PF\d+\">(PF\d+)<', urlSource)
pfamDesc = re.findall(r'<td class\=\"desc\">(.*?)</td>', urlSource)

outfile = open('%s/transposase_pfams.txt' % rootDir, 'w')
if len(pfam) == len(pfamDesc) :
    print('length matches')
    pfamDict = {}
    for p in pfamDesc :
        if re.search(r'[Tt]ransposase', p) :
            pfamDict[pfam[p]] = pfamDesc[p]
        else :
            pass
    for p in pfam :
        if p in pfamDict :
            outfile.write(p + '\t' + pfamDict[p] + '\n')
        else :
            pass



outfile.close()

print('\n-----DONE-----\n')