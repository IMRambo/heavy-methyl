
# Test: download multiple genome and gene nucleotide sequences from NCBI.
# Ian Rambo
# September 10, 2015

'''
GOAL: download UIDs from NCBI Gene database. Biopython 9.15.1
Use Entrez.elink to get nucleotide sequences from Gene db. Biopython 9.7
Create output FASTA containing nucleotide sequences for all current
dam methyltransferases. Biopython 9.15.1
'''
from Bio import SeqIO, Entrez

#Sulfur and sulfate reducing deltaproteobacteria orders
dpbOrders = ['Desulfovibrionales','Desulfobacterales',
'Desulfurellales','Desulfuromonadales']

dpbGI = list()
dpbGI_count = list()
dpbWebEnv = list()
dpbQueryKey = list()

#Hello, NCBI
Entrez.email = 'imrambo@udel.edu'
for i in range(len(dpbOrders)) :
    #Search for deltaproteobacteria
    dpbSearchHandle = Entrez.esearch(db = 'nucleotide',
    term = '%s[Orgn] AND srcdb_refseq[prop] AND biomol_genomic[prop] AND\
    (biomol_genomic[prop] AND refseq[filter])' % dpbOrders[i],
    retmax = 10000, usehistory = 'y')
    dpbSearchResults = Entrez.read(dpbSearchHandle)
    dpbSearchHandle.close()
    dpbGI.append(dpbSearchResults['IdList'])
    dpbGI_count.append(dpbSearchResults['Count'])
    dpbWebEnv.append(dpbSearchResults['WebEnv'])
    dpbQueryKey.append(dpbSearchResults['QueryKey'])
    
#Convert elements in list dpbGI_count from strings to integers
dpbGI_count = [ int(i) for i in dpbGI_count ]
print(sum(dpbGI_count))

#Make sure that record counts and number of GI entries are equal
for i in range(0, len(dpbGI)) :
    assert dpbGI_count[i] == len(dpbGI[i])

#assert len(dpbWebEnv) == len(dpbQueryKey)

import os
os.chdir('/Users/imrambo/Documents/BINF868/seq_data/')
outHandle = open('deltaproteobacteriaGenomes_RefSeq.fa', 'w')
batch_size = 4
for i in range(0, len(dpbGI)) :
    for start in range(0, sum(dpbGI_count), batch_size) :
        end = min(dpbGI_count[i], start + batch_size)
	print('Downloading genomic sequence %i to %i' % (start + 1, end))
	Entrez.email = 'imrambo@udel.edu'
	dpbFetchHandle = Entrez.efetch(db = 'nucleotide', rettype = 'fasta',
	retmode = 'text', retstart = start, retmax = batch_size,
	id = ','.join(dpbGI[i]))
	dpbFetchRecord = dpbFetchHandle.read()
	dpbFetchHandle.close()
	dpbFetchRecord = dpbFetchRecord.rstrip()
	outHandle.write(dpbFetchRecord)

#outHandle.close()

#print(type(dpbWebEnv[0]), len(dpbWebEnv))
# print(str(dpbWebEnv))
# print(type(dpbWebEnv))
# myList = [str(dpbWebEnv), 'hallo']
# print(myList)
#Instead of using a list of GI numbers, the WebEnv session cookie and
#query key will be used.
#dpbWebEnv = str(dpbSearchResults['WebEnv'])
# print(type(dpbWebEnv))
# print(dpbWebEnv)
#dpbQueryKey = dpbSearchResults['QueryKey']

#print(dpbWebEnv)





# Entrez.email = 'imrambo@udel.edu'
# handle = Entrez.efetch(db = 'gene', id = '1232567', rettype = 'fasta',
# retmode = 'text')
# record = SeqIO.read(handle, 'genbank')
# print(str(record.seq))

# Entrez.email = 'imrambo@udel.edu'
# damSearchHandle = Entrez.esearch(db = 'Gene',
# #Search for bacterial dam methyltransferases, Current records (alive)
# term = 'bacteria[Orgn] AND dam[Gene] AND alive[prop]',
# retmax = 100, usehistory = 'y')
#
# damSearchResults = Entrez.read(damSearchHandle)
# damSearchHandle.close()
#
# Entrez.email = 'imrambo@udel.edu'
# handle = Entrez.elink(dbfrom = 'gene',
# id = '947893', linkname = 'gene_nuccore')
# record = Entrez.read(handle)
# handle.close()
# print(record[0]['LinkSetDb'][0]['Link'][0])
# print(record[0]['LinkSetDb'][0]['Link'][1])
#
# handle = Entrez.efetch(db = 'gene', id = '1232567',
# rettype = 'fasta', retmode = 'text')
# record = SeqIO.read(handle, 'fasta')

# UIDList = damSearchResults['IdList']
# UIDCount = int(damSearchResults['Count'])
# print(UIDCount)
# assert UIDCount == len(UIDList)

#Instead of using a list of UID numbers, the WebEnv session cookie and
#query key can be used. Biopython tutorial 9.15.1
# damWebEnv = damSearchResults['WebEnv']
# damQueryKey = damSearchResults['QueryKey']

# record = Entrez.read(Entrez.elink(dbfrom = 'Gene', db = 'nucleotide',
# id = '947893'))
#
# for linksetdb in record[0]["LinkSetDb"]:
#     print(linksetdb["DbTo"], linksetdb["LinkName"], len(linksetdb["Link"]))
#
# print(record[0]['LinkSetDb'][0]['Link'][0])
# print(record[0]['LinkSetDb'][0]['Link'][1])
# batch_size = 2
#
# import os
# os.chdir('/Users/imrambo/Documents/BINF868')
# damFASTA = open('damMTase.fa', 'w')



# damSearchRecord = Entrez.read(damSearchHandle)
# #dam_handle.close()
# print(damSearchRecord["IdList"])
# print(len(damSearchRecord["IdList"]))


#dpbWebEnv = dpbSearchResults['WebEnv']
#dpbQueryKey = dpbSearchResults['QueryKey']
# for i in range(len(dpbWebEnv)) :
#     print(dpbWebEnv[i])
#     print(type(dpbWebEnv[i]))
# for i in range(len(dpbQueryKey)) :
#     print(type(dpbQueryKey[i]))
#     print(dpbQueryKey[i])
# print(len(dpbQueryKey))
# print(dpbQueryKey)
