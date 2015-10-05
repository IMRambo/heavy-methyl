#!/usr/bin/env python
#Ian Rambo
#BINF868
#University of Delaware
#Last updated: October 4, 2015
#==============================================================================
'''
GOAL:
Download nucleotide sequences for bacterial DNA methyltransferase genes.
METHOD:
Download nucleotide sequences from NCBI ENTREZ with urllib.
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
def get_gene_gi(accession, start, end) :
    '''
    Parse gene nucleotide GI numbers from NCBI using gene Accession numbers. 
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
    #Entrez.email = "imrambo@udel.edu"
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
        Headers = {}
        #Firefox 40.1
        Headers['User-Agent'] = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0) Gecko/20100101 Firefox/40.1'
        #Request the URL with the variables
        req = urllib.request.Request(url, data, headers = Headers)
        #Response - visit the URL
        resp = urllib.request.urlopen(req)
        respData = resp.read().decode('utf-8')
        respData = str(respData).replace('\n\n', '\n')
        outputHandle.write(respData)
    outputHandle.close()
    return
#------------------------------------------------------------------------------
def efetch_genome_nucleotide(GI, batch_size, filename) :
    '''
    Download genomic nucleotide sequences in FASTA format using the NCBI
    efetch eutil via urllib.
    '''
    #Import Entrez to tell NCBI who you are with Entrez.email
    Entrez.email = 'imrambo@udel.edu'
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
    outputHandle.close()
    return
#------------------------------------------------------------------------------
def fasta_write_test(fasta) :
    '''
    Test if FASTA output file wrote correctly. 
    '''
    with open(fasta, 'r') as fa :
        lines = fa.readlines()
        if lines[0].startswith('>') and lines[0].endswith('\n') and not '>' in lines[-1] :
            print('FASTA file written successfully.')
        else :
            print('FASTA file not written successfully.')
    fa.close()
    return
#==============================================================================
#START
#==============================================================================
#BACTERIAL DNA METHYLTRANSFERASES
#==============================================================================
'''
#Parse txt file from NCBI gene search
os.chdir('/Users/imrambo/Documents/BINF868/')
accession = []
start = []
stop = []
#txt file containing gene search results, downloaded from NCBI
geneInfo = 'dcm_genome_result.txt'


nonTerms = ['hypothetical protein', 'RNA methyltransferase', 'putative']
with open(geneInfo, 'r') as gres :
    for record in gres :
        record = record.rstrip()
        recordL = record.split('\t')
        if not any(x in recordL[7] for x in nonTerms) and len(recordL) == 16 :
            accession.append(recordL[11])
            start.append(recordL[12])
            stop.append(recordL[13])
           
        
gres.close()

batch_size = 500
mtgi = get_gene_gi(accession, start, stop)

outputHandle = open('dcm_genome.fa', 'w')
url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
assert len(mtgi) == len(start) == len(stop)
print('Downloading nucleotide sequences...')
for i in range(0, len(mtgi)) :
    values = {'db':'nuccore','id':mtgi[i],
              'seq_start':start[i], 'seq_stop':stop[i],
              'rettype':'fasta','retmax':str(batch_size)}
    data = urllib.parse.urlencode(values)
    data = data.encode('utf-8')
    Headers = {}
    #Firefox 40.1
    #Headers['User-Agent'] = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0) Gecko/20100101 Firefox/40.1'
    #Request the URL with the variables
    req = urllib.request.Request(url, data)#, headers = Headers)
    #Response - visit the URL
    resp = urllib.request.urlopen(req)
    respData = resp.read().decode('utf-8')
    respData = str(respData).replace('\n\n', '\n')
    outputHandle.write(respData)

outputHandle.close()

print('Download successful.')
'''
#==============================================================================
#HTML regex parsing method. Caveat - I don't know how to download
#more than 20 records at a time. 
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
mtaseGI = get_gene_gi(accession, locusStart, locusEnd)

os.chdir('/Users/imrambo/Documents/BINF868/')
#Fetch nucleotide sequences for dam MTases
efetch_gene_nucleotide(mtaseGI, locusStart, locusEnd, 5, 'dam.fa')
'''
#==============================================================================
#TARGET GENOMES
#==============================================================================
os.chdir('/Users/imrambo/Documents/BINF868')



