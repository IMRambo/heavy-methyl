#!/usr/bin/env python
#Ian Rambo
#BINF868
#University of Delaware
#Last updated: September 18, 2015
#==============================================================================
'''
GOAL:
Download nucleotide sequences for bacterial DNA methyltransferase genes.
METHOD:
Parse NCBI databases using urllib. 
'''
#==============================================================================
import re
import os
import urllib.request
import urllib.parse

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
#==============================================================================


##==============================================================================
##Count potential methylation motifs in target genomes
#print('CCGG counts, Dehalococcoides mccartyi')
#find_motifs('CCGG', Dmccartyi_genome)
#print('GANTC counts, Dehalococcoides mccartyi')
#find_motifs('GA[TCAG]TC', Dmccartyi_genome)
#print('GATC counts, Dehalococcoides mccartyi')
#find_motifs('GATC', Dmccartyi_genome)
#
#print('CCGG, Desulfocapsa sulfexigens genome')
#find_motifs('CCGG', Dsulfexigens_genome)
#print('GANTC, Desulfocapsa sulfexigens genome')
#find_motifs('GA[TCAG]TC', Dsulfexigens_genome)
#print('GATC, Desulfocapsa sulfexigens genome')
#find_motifs('GATC', Dsulfexigens_genome)
##==============================================================================
#myMotifs = ['CCGG', 'GA[TCAG]TC', 'GATC']
#find_motifs(myMotifs, Dmccartyi_genome)
