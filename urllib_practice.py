#!/usr/bin/env python

#Ian Rambo
#BINF868
#Last updated: September 10, 2015
#note: must use python3 command via terminal to run python 3
#urllib requires python 3


'''
GOAL:
Use urllib to get data from the web.
'''
#Tutorial section: Python 3 Programming Tutorial #35 - urllib module
#YouTube - sentdex
import urllib.request
import urllib.parse

#Make request to the URL
#x = urllib.request.urlopen('https://www.google.com')

#print(x.read())
'''
#Let's get data from NCBI
#Tutorial example has changed
url = 'http://www.ncbi.nlm.nih.gov/gquery/'
values = {'term':'Desulfobacterales'}

data = urllib.parse.urlencode(values)
data = data.encode('utf-8')
#Request the URL with the variables
req = urllib.request.Request(url, data)
#Response - visit the URL
resp = urllib.request.urlopen(req)
respData = resp.read()

print(respData)
'''
'''
try :
    x = urllib.request.urlopen('https://www.google.com/search?q=test')

    print(x.read())

except Exception as e :
    print(str(e))
#403 code: Forbidden. Google sees that it's a program.

try :
    url = 'https://www.google.com/search?q=test'
    #Information on you
    headers = {}
    #User-Agent - type of browser. This looks like Internet Explorer.
    #Allows you to bypass filters for robots.
    headers['User-Agent'] = 'Mozilla/5.0 (X11; Linux i686) AppleWebKit/537.17\
    (KHTML, like Gecko) Chrome/24.0.1312.27 Safari/537.17'
    req = urllib.request.Request(url, headers = headers)
    #Response
    resp = urllib.request.urlopen(req)
    respData = resp.read()

    saveFile = open('withHeaders.txt', 'w')
    saveFile.write(str(respData))
    saveFile.close()

except Exception as e :
    print(str(e))
'''
#==============================================================================
#Tutorial section: Python 3 Programming Tutorial #36 - Regular Expressions w/re
#YouTube - sentdex
#This will apply to the next tutorial- using regex to parse websites
'''
Identifiers:
\d any number
\D anything but a number
\s space
\S anything but a space
\w any character
\W anything but a character
. = any character EXCEPT newline
\b whitespace around words
\. a period

Modifiers:
{1, 3} we're expecting 1-3
+ match 1 or more
? match 0 or 1
* match 0 or more
$ match end of string
^ match beginning of string
| either or
[] range or "variance" e.g. [A-Za-z] first uppercase, second lower
{x} expect "x" amount

White Space Characters:
\n newline *
\s space *
\t tab *
\e escape
\f form feed
\r return

Don't forget!:
If you want to use these:
. + * ? [ ] $ ^ () {} | \
You must escape them with \
'''
'''
import re

exampleString = 'Jessica is 15 years old, and Daniel is 27 years old.\
Edward is 97, and his grandfather, Oscar, is 102.'

#Pull the names and ages from string.

ages = re.findall(r'\d{1,3}', exampleString)
#Find one capital letter, and 0 or more lower case letters
names = re.findall(r'[A-Z][a-z]*', exampleString)

print(ages)
print(names)
ageDict = dict()
x = 0
for eachName in names :
    ageDict[eachName] = ages[x]
    x += 1
'''
#==============================================================================
'''
#Tutorial section: Python 3 Programming Tutorial #37
#Parse websites with re and urllib
#YouTube - sentdex
import urllib.request
import urllib.parse
import re

url = 'http://pythonprogramming.net/parse-website-using-regular-expressions-urllib/'
# data = urllib.parse.urlencode(values)
# data = data.encode('utf-8')
req = urllib.request.Request(url)
resp = urllib.request.urlopen(req)
respData = resp.read()

#print(respData)
#Find everything between paragraph tags
paragraphs = re.findall(r'<p>(.*?)</p>', str(respData))

#for eachP in paragraphs :
    #print(eachP)
'''
#==============================================================================
'''
#Parse FASTA from NCBI. Nucleotide database, E. coli K-12 Dam MTase
url = 'http://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3'
values = {'report':'fasta','from':'3515077',
          'to':'3515913','strand':'true'}
data = urllib.parse.urlencode(values)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
resp = urllib.request.urlopen(req)
respData = resp.read().decode('utf-8')

print(respData)
'''
'''
url = urllib.request.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?')
values = {'dbfrom':'gene','db':'nuccore',
          'id':'947893','term':'srcdb+refseq[prop]'
          }

data = urllib.parse.urlencode(values)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
resp = urllib.request.urlopen(req)
respData = resp.read()
print(respData)
'''

'''
#Let's get data from NCBI
#Tutorial example has changed
url = 'http://www.ncbi.nlm.nih.gov/gquery/'
values = {'term':'Desulfobacterales'}

data = urllib.parse.urlencode(values)
data = data.encode('utf-8')
#Request the URL with the variables
req = urllib.request.Request(url, data)
#Response - visit the URL
resp = urllib.request.urlopen(req)
respData = resp.read()
'''
import re
import os

#Search NCBI for bacterial methyltransferase gene ID, Accession number, strand, and loci.
url = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/gene/?term=bacteria[Orgn]+AND+dam[Gene]+AND+alive[prop]')
urlSource = url.read()
#Gene UID
geneID = re.findall(r'>ID:\s*(\d+)<', str(urlSource))
#GenBank Accession numbers
accession = re.findall(r'<td>\D*(NC_\d+\.?\d+)\s*', str(urlSource))
#Gene location information; start locus, stop locus, strand
location = re.findall(r'\s*\((\d+\.\.\d+\,?\s*\w*)\)', str(urlSource))
#Gene name - finding more than 20 hits. Why?
#geneName = re.findall(r'style="background-color\:">(.*?)</span>', str(urlSource))
#geneName = [x.lower() for x in geneName]
#==============
strand = []
for c in location :
    if 'complement' in c :
        strand.extend('2')
    else :
        strand.extend('1')

locusStart = []
for s in location :
    locusStart.extend(re.findall(r'^\d+', s))

locusEnd = []
for e in location :
    locusEnd.extend(re.findall(r'\.\.(\d+)\,?', e))
    
assert len(geneID) == len(accession) == len(strand) == len(locusStart) == len(locusEnd) 

#List of GI numbers for dam MTase nucleotide records
GI = []
for i in range(0, len(accession)) :
    print('Parsing GI number for accession %s' % accession[i])
    url = 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % accession[i]
    values = {'report':'fasta','from':locusStart[i],
              'to':locusEnd[i],'strand':'true'}
    data = urllib.parse.urlencode(values)
    data = data.encode('utf-8')
    #Request the URL with the variables
    req = urllib.request.Request(url, data)
    #Response - visit the URL
    resp = urllib.request.urlopen(req)
    respData = resp.read()
    #Return first instance of regex match for GI number
    uid = re.findall(r'uid\=(\d+)\&\w*\;', str(respData))[0]
    GI.append(uid)

#Fetch nucleotide sequences for dam MTases
batch_size = 4
os.chdir('/Users/imrambo/Documents/BINF868/')

outputHandle = open('dam.fa', 'w')

for i in range(0, len(GI)) :
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    values = {'db':'nuccore','id':GI[i],
              'seq_start':locusStart[i], 'seq_stop':locusEnd[i],
              'rettype':'fasta','retmax':str(batch_size)}
    data = urllib.parse.urlencode(values)
    data = data.encode('utf-8')
    #Request the URL with the variables
    req = urllib.request.Request(url, data)
    #Response - visit the URL
    resp = urllib.request.urlopen(req)
    respData = resp.read().decode('utf-8')
    respData = str(respData).replace('\n\n', '\n')
    outputHandle.write(respData)