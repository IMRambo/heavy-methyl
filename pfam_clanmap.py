#!/usr/bin/env python
import urllib.request
import urllib.parse
import re

'''
Ian Rambo
Date created: January 22, 2016
Last modified: January 22, 2016

GOAL: parse Pfam clan entries from pfam.xfam.org to create a
mapping file of Pfam clans and their member families.
'''
#==============================================================================
#-----GLOBAL VARIABLES-----
outDir = '/Users/imrambo/R/ORMP/data'
#Pfam clan mapping file
clanMapFile = '%s/PfamClanMapping.txt' % outDir
#Number of clans 
clanCount = 0
#==============================================================================
print('-----BEGIN-----\n')

allClansURL = urllib.request.urlopen('http://pfam.xfam.org/clan/browse')
urlSource = allClansURL.read()

#Clan accession numbers
clanAccessions = re.findall(r'http://pfam.xfam.org/clan/(CL\d+)">', str(urlSource))

with open(clanMapFile, 'w') as clanMap :
    clanMap.write('PfamClan\tPfamFamily\n')
    for clan in clanAccessions :
        print('Finding members of %s' % clan)
        clanURL = ''
        clanMembers = []
        clanURL = 'http://pfam.xfam.org/clan/%s' % clan
        clanURLReq = urllib.request.urlopen(clanURL)
        clanURLSource = str(clanURLReq.read()).strip()
        clanMembers = re.findall(r'<a title="PF\d+"\\n\s+href=".*?">\\n\s+(.*?)</a>', clanURLSource)
        #print(clanMembers)
        for member in clanMembers :
            clanMap.write('%s\t%s\n' % (clan, member))
        
        clanCount += 1

clanMap.close()
print('Parsed data for %d clans\n' % clanCount)
print('-----DONE-----\n')