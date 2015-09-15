#!/usr/bin/env python
#Ian Rambo
#BINF868
#University of Delaware
#Assignment 1: target genomes, MTases, methylation motifs
#Last updated: September 15, 2015
#==============================================================================
'''
DETERMINE POTENTIAL METHYLATION SITES IN GENOMES OF INTEREST
GOALS:
Download nucleotide sequences for RefSeq target genomes
Download nucleotide sequences for bacterial methyltransferase genes
Count potential sites for cytosine and adenine methylation in genomes
'''
import urllib.request
import urllib.parse
import re
import os
os.chdir('/Users/imrambo/Documents/BINF868')
#==============================================================================
#FUNCTIONS
def find_motifs(pattern, genome):
    '''
    Find the number of potential methylation motifs within a genome.
    '''
    if type(genome) is list :
        for i in range(len(genome)) :
            if type(genome[i]) is str and type(pattern) is str :
                print(len(re.findall(pattern, genome[i])))
            elif type(genome[i]) is str and type(pattern) is list :
                pattern = '|'.join(pattern)
                print(len(re.findall(pattern, genome[i])))
            else :
                print('Genome sequence is not a string')
                break

    elif type(genome) is str :
        print(len(re.findall(pattern, genome)))

    elif type(genome) is not str :
        print('Genome sequence is not a string')

def efetch_genome_nucleotide(gi, batch_size, outputHandle) :
    '''
    Download multiple genome nucleotide sequences in FASTA format from NCBI
    using the Entrez efetch utility. Write FASTA sequences to output file.
    '''
    print('Fetching sequences...')
    for i in range(0, len(gi)) :
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        values = {'db':'nuccore','id':gi[i],
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
        
def efetch_gene_nucleotide(gi, start, stop, name, batch_size, outputHandle) :
    '''
    Download multiple gene nucleotide sequences in FASTA format from NCBI
    using the Entrez efetch utility. Write FASTA sequences to an output
    file.
    '''
    assert len(gi) == len(start) == len(stop) == len(name)
    print('Fetching gene nucleotide sequences from NCBI...')
    for i in range(0, len(gi)) :
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        values = {'db':'nuccore','id':gi[i],
                  'seq_start':start[i], 'seq_stop':stop[i],
                  'rettype':'fasta','retmax':str(batch_size)}
        data = urllib.parse.urlencode(values)
        data = data.encode('utf-8')
        #Request the URL with the variables
        req = urllib.request.Request(url, data)
        #Response - visit the URL
        resp = urllib.request.urlopen(req)
        respData = resp.read().decode('utf-8')
        respData = str(respData).replace('\n\n', '\n')
        #outputHandle.write(respData)
        #Add gene name to FASTA header
        if respData.startswith('>') :
            respData = re.sub(r'^>\w+\|\d+\:\d+\-\d+.*\n$',
                              r'^>\w+\|\d+\:\d+\-\d+.*\|%s\n$'% geneName[i], respData)
            outputHandle.write(respData)
        else :
            outputHandle.write(respData)
#==============================================================================
#TARGET GENOME NUCLEOTIDE SEQUENCES
#==============================================================================
'''
Parse NCBI for bacterial target RefSeq genomes.
'''
#Deltaproteobacteria
#order Desulfovibrionales
taxon = 'desulfovibrionales'
url = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/nuccore/?term=desulfovibrionales[Orgn]+AND+srcdb_refseq[prop]+AND+biomol_genomic[prop]+AND+%28biomol_genomic[prop]+AND+refseq[filter]%29+AND+srcdb_refseq[PROP]')
urlSource = url.read()
genomeGI = re.findall(r'>Select\s?item\s?(\d+)</label>', str(urlSource))

genomeOutHandle = open('desulfovibrionales.fa', 'w')

efetch_genome_nucleotide(genomeGI, 5, genomeOutHandle)  
genomeOutHandle.close()
#==============================================================================
#BACTERIAL METHYLTRANSFERASE NUCLEOTIDE SEQUENCES
#==============================================================================
'''
Parse NCBI to retrieve data for bacterial methyltransferase genes.
'''
#ALL BACTERIAL METHYLTRANSFERASES STORED ON NCBI GENE DATABASE
#url = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/gene/?term="Bacteria"[Organism]+AND+methyltransferase[All+Fields]+AND+("source_genomic"[properties]+AND+"genetype+protein+coding"[Properties]+AND+alive[prop])')
#Search for dcm methyltransferase
url = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/gene/?term=bacteria[Orgn]+AND+dcm[Gene]+AND+alive[prop]')
urlSource = url.read()
#Gene UID
geneID = re.findall(r'>ID:\s*(\d+)<', str(urlSource))
#Gene name
geneName = re.findall(r'link_uid\=\d+\">(.*?)</a>', str(urlSource))
#GenBank Accession numbers
accession = re.findall(r'<td>\D*(NC_\d+\.?\d+)\s*', str(urlSource))
#Gene location information: start locus, stop locus, strand
location = re.findall(r'\s*\((\d+\.\.\d+\,?\s*\w*)\)', str(urlSource))
#Gene description
#geneDesc = re.findall(r'<td>.*\s?\[<em>\w+</em>\]</td>', str(urlSource))

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

#Parse GI numbers and place into list
GI = []
assert len(accession) == len(locusStart) == len(locusEnd)
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

geneOutHandle = open('dcm.fa', 'w')

efetch_gene_nucleotide(GI, locusStart, locusEnd, geneName, 5, geneOutHandle)

geneOutHandle.close()

#==============================================================================
#SCRAPYARD
##==============================================================================
##==============================================================================
#
##TARGET GENOMES
##Download complete genomes from NCBI RefSeq.
#
#'''
#Dehalococcoides mccartyi strain SG1 genome
#Illumina Hi-Seq 2000 shotgun sequencing
#Assembled in SOAPdenovo, scaffolds constructed in Opera
#Wang et al. 2014
#'''
##GI numbers for scaffolds 1-5
#Dmccartyi_GI = ['736351273', '36350382', '736350119', '736350202', '736350057']
##Build complete D. mccartyi genome from scaffolds
#Dmccartyi_genome = list()
#
#for i in range(len(Dmccartyi_GI)) :
#    Entrez.email = 'imrambo@udel.edu'
#    scaffold_h = Entrez.efetch(db = 'nucleotide',
#                                        id = Dmccartyi_GI[i],
#                                        rettype = 'fasta')
#    scaffold = SeqIO.read(scaffold_h, 'fasta')
#    scaffold_h.close()
#    Dmccartyi_genome.append(str(scaffold.seq))
##------------------------------------------------------------------------------
#'''
#Desulfocapsa sulfexigens DSM 10523 genome
#454 GS FLX Titanium pyrosequencing
#Assembled in Newbler v2.5.3, gaps in assembly closed by 259 Sanger reads
#Finster et al. 2013
#'''
#Entrez.email = 'imrambo@udel.edu'
#Dsulfexigens_h = Entrez.efetch(db = 'nucleotide',
#id = '451945650',
#rettype = 'fasta')
#
#Dsulfexigens_fasta = SeqIO.read(Dsulfexigens_h, 'fasta')
#Dsulfexigens_h.close()
#Dsulfexigens_genome = str(Dsulfexigens_fasta.seq)
##------------------------------------------------------------------------------
#'''
#Desulfovibrio desulfuricans ND132 complete genome
#Whole genome shotgun sequencing
#Assembly in Newbler v.2.3 (pre-release)
#'''
#Entrez.email = 'imrambo@udel.edu'
#Desulfovibrio_h = Entrez.efetch(db = 'nucleotide',
#id = '376294792',
#rettype = 'fasta')
#
#Desulfovibrio_fasta = SeqIO.read(Desulfovibrio_h, 'fasta')
#Desulfovibrio_h.close()
#Desulfovibrio_desulfuricans_genome = str(Desulfovibrio_fasta.seq)
#=====================================================================
##BACTERIAL METHYLTRANSFERASES
#
##Fetch adenine MTase sequence files from NCBI's Entrez databases
##File handles
##Cell cycle regulated m6A MTase (CcrM), Caulobacter crescentus
#Entrez.email = 'imrambo@udel.edu'
#CcrM_h = Entrez.efetch(db = 'nucleotide',
#                       id = '221232939',
#                       rettype = 'fasta',
#                       seq_start = '398637',
#                       seq_stop = '399713')
##Dam m6A MTase, E.coli K12
#Entrez.email = 'imrambo@udel.edu'
#Dam_h = Entrez.efetch(db = 'nucleotide',
#                    id = '556503834',
#                    rettype = 'fasta',
#                    seq_start = '3515913',
#                    seq_stop = '3515077')
##Dcm m5C MTase, E.coli K12
#Entrez.email = 'imrambo@udel.edu'
#Dcm_h = Entrez.efetch(db = 'nucleotide',
#                    id = '556503834',
#                    rettype = 'fasta',
#                    seq_start = '2032317',
#                    seq_stop = '2030899')
#
#CcrM = SeqIO.read(CcrM_h, 'fasta')
#Dam = SeqIO.read(Dam_h, 'fasta')
#Dcm = SeqIO.read(Dcm_h, 'fasta')
#
#CcrM_h.close()
#Dam_h.close()
#Dcm_h.close()
#
#MTDict = {'CcrM':CcrM.seq,
#          'Dam':Dam.seq,
#          'Dcm':Dcm.seq}
#==============================================================================

#==============================================================================
#Count potential methylation motifs in target genomes
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