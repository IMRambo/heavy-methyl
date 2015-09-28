#!/usr/bin/env python
#Ian Rambo
#BINF868
#University of Delaware
#Last updated: September 23, 2015
#==============================================================================
'''
GOAL:
Download nucleotide sequences for bacterial DNA methyltransferase genes.
METHOD:
Parse HTML source code from NCBI using regular expressions. 
'''
#==============================================================================
import re
import os
import urllib.request
import urllib.parse
#==============================================================================
#Functions
#------------------------------------------------------------------------------
def get_gene_gi_HTML(accession, start, end) :
    '''
    Parse gene nucleotide GI numbers from NCBI HTML using gene Accession numbers. 
    '''
    geneGI = []
    assert len(accession) == len(start) == len(end)
    for i in range(0, len(accession)) :
        print('Parsing GI number for accession %s' % accession[i])
        url = 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % accession[i]
        values = {'report':'fasta','from':start[i],
                  'to':end[i],'strand':'true'}
        data = urllib.parse.urlencode(values)
        data = data.encode('utf-8')
        #Request the URL with the variables
        req = urllib.request.Request(url, data)
        #Response - visit the URL
        resp = urllib.request.urlopen(req)
        respData = resp.read()
        #Return first instance of regex match for GI number
        uid = re.findall(r'uid\=(\d+)\&\w*\;', str(respData))[0]
        geneGI.append(uid)
    return(geneGI)
#------------------------------------------------------------------------------
def efetch_gene_nucleotide(GI, start, end, batch_size, filename) :
    '''
    Download gene nucleotide sequences in FASTA format using the NCBI
    efetch eutil. 
    '''
    print('\nDownloading gene sequences...\n')
    outputHandle = open('%s' % filename, 'w')
    assert len(GI) == len(start) == len(end)
    for i in range(0, len(GI)) :
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        values = {'db':'nuccore','id':GI[i],
                  'seq_start':start[i], 'seq_stop':end[i],
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
#------------------------------------------------------------------------------
def efetch_genome_nucleotide(GI, batch_size, filename) :
    '''
    Download genomic nucleotide sequences in FASTA format using the NCBI
    efetch eutil. 
    '''
    print('\nDownloading genome sequences...\n')
    outputHandle = open('%s' % filename, 'w')
    for i in range(0, len(GI)) :
            url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
            values = {'db':'nuccore','id':GI[i],
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
'''
#Search NCBI for bacterial methyltransferase
#gene IDs, Accession numbers, strands, and loci.
gene = 'dam'
url = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/gene/?term=bacteria[Orgn]+AND+' + gene + '[Gene]+AND+alive[prop]')
HTML = url.read()
#Gene UID
geneID = re.findall(r'>ID:\s*(\d+)<', str(HTML))
#GenBank Accession numbers
accession = re.findall(r'<td>\D*(NC_\d+\.?\d+)\s*', str(HTML))
#Gene location information; start locus, stop locus, strand
location = re.findall(r'\s*\((\d+\.\.\d+\,?\s*\w*)\)', str(HTML))
#Gene name - finding more than 20 hits. Why?
#geneName = re.findall(r'style="background-color\:">(.*?)</span>', str(urlSource))
#geneName = [x.lower() for x in geneName]
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------    
#List of GI numbers for dam MTase nucleotide records
mtaseGI = get_gene_gi_HTML(accession, locusStart, locusEnd)

os.chdir('/Users/imrambo/Documents/BINF868/')
#Fetch nucleotide sequences for dam MTases
efetch_gene_nucleotide(mtaseGI, locusStart, locusEnd, 5, 'dam.fa')
'''
#==============================================================================
#Bacterial target genomes
#Define taxon of interest for NCBI search
taxon = 'Desulfuromonadales'

url = 'http://www.ncbi.nlm.nih.gov/nuccore/?term=' + taxon + '[Orgn]+AND+%28biomol_genomic[PROP]+AND+refseq[filter]%29'
#Search results for taxon
URL = urllib.request.urlopen(url)
HTML = URL.read()
#Parse nucleotide GI numbers
genomeGI = re.findall(r'GI:</dt>\s*<dd>(\d+)</dd>', str(HTML))

print('\nParsing %s genome GIs...\n' % taxon)
efetch_genome_nucleotide(genomeGI, 5, '%s_genomes.fa' % taxon)




#==============================================================================
#Parse txt file from NCBI gene search
accession = []
start = []
stop = []
with open('mtase_gene_results_test.txt', 'r') as gres :
    for record in gres :
        record = record.rstrip()
        recordl = record.split()
        accession.append(recordl[0])
        start.append(recordl[1])
        stop.append(recordl[2])
        
gres.close()

print(accession)
        