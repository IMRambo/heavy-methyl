#Ian Rambo
#Last updated: October 8, 2015
#==============================================================================
import os
import re
import urllib.request
import urllib.parse
#==============================================================================
'''
GOAL: Obtain genome nucleotide sequences from NCBI in FASTA format using the
NCBI eutils and urllib.
'''
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/
#==============================================================================
def fasta_write_test(fasta) :
    '''
    Test if FASTA output file wrote correctly. 
    '''
    with open(fasta, 'r') as fa :
        lines = fa.readlines()
        if lines[0].startswith('>') and lines[0].endswith('\n') and not '>' in lines[-1] :
            print('FASTA file %s written successfully.' % fasta)
        else :
            print('FASTA file %s not written successfully.' % fasta)
    fa.close()
    return
#==============================================================================
os.chdir('/Users/imrambo/Documents/BINF868')

taxFile = 'targetGenusCleanSort.txt'
targetTaxa = []
with open(taxFile, 'r') as taxa :
    for t in taxa :
        t = t.rstrip()
        targetTaxa.append(t)
        #searchTerms.extend(t)
    #return searchTerms
taxa.close()

#testTaxa = targetTaxa[1:5]      
#
#print(testTaxa)
#------------------------------------------------------------------------------
#Esearch
#------------------------------------------------------------------------------
#Enter target taxon - script will work for any taxonomic level
taxon = 'Zymomonas'
#Search for complete genome sequences
esearch_term = '"bacteria"[orgn] AND "%s"[orgn] AND "complete genome"[title] AND "biomol_genomic"[prop]' % taxon

try :
    xml = 'esearch_xml/%s_esearch.xml' % taxon
    outputHandle = open(xml, 'w')  
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    values = {'db':'nuccore','term':esearch_term,
              'rettype':'uilist','retmode':'uilist',
              'retmax':'9000','usehistory':'y'}
    data = urllib.parse.urlencode(values)
    data = data.encode('utf-8')
    Headers = {}
    #Firefox 40.1 User Agent
    Headers['User-Agent'] = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0)\
                             Gecko/20100101 Firefox/40.1'
    #Request the URL with the variables
    req = urllib.request.Request(url, data, headers = Headers)
    #Response - visit the URL
    resp = urllib.request.urlopen(req)
    respData = resp.read().decode('utf-8')
    respData = str(respData)
    outputHandle.write(str(respData))
    outputHandle.close()
    
except Exception as e :
    print(str(e))
#------------------------------------------------------------------------------
ID = []
queryKey = []
webEnv = []
for filename in os.listdir("./esearch_xml") :
    if filename.startswith(taxon) :
        print(filename)
        with open('./esearch_xml/%s' % filename, 'r') as xml :
            for line in xml :
                line = line.rstrip()
                qk = re.findall(r'<QueryKey>(\d+)</QueryKey>', line)
                queryKey.extend(qk)
                we = re.findall(r'<WebEnv>(.*?)</WebEnv>', line)
                webEnv.extend(we)
                gid = re.findall(r'<Id>(\d+)</Id>', line)
                ID.extend(gid)
        xml.close()
print(len(ID))
#print(queryKey)
#print(webEnv)

#------------------------------------------------------------------------------
#Efetch
#------------------------------------------------------------------------------
os.chdir('/Users/imrambo/Documents/BINF868/TargetGenomes')
print('Efetch')
targetGenomes = '%sGenomes.fa' % taxon
outputHandle = open(targetGenomes, 'w')

count = int(len(ID))
batchSize = 40
for start in range(0, count, batchSize) :
    end = min(count, start+batchSize)
    print('Preparing to download record %i to %i' % (start+1, end))
    try :
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        values = {'db':'nuccore','query_key':str(queryKey[0]),
                  'WebEnv':str(webEnv[0]),'rettype':'fasta',
                  'retmode':'fasta','retstart':str(start),
                  'retmax':str(batchSize)}
                     
        data = urllib.parse.urlencode(values)
        data = data.encode('utf-8')
        Headers = {}
        #Firefox 40.1 User Agent
        Headers['User-Agent'] = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0)\
                                 Gecko/20100101 Firefox/40.1'
        #Request the URL with the variables
        #print('Requesting URL...')
        req = urllib.request.Request(url, data, headers = Headers)
        #Response - visit the URL
        resp = urllib.request.urlopen(req)
        respData = resp.read().decode('utf-8')
        respData = str(respData)
        #print('Downloading FASTA sequences...')
        outputHandle.write(str(respData))
           
            
    except Exception as e :
        print(str(e))
#print("Download complete")

outputHandle.close()
print('Running filter_empty_seqs.awk ... please be patient')
filterGenomes = '%sGenomesFilter.fa' % taxon
awk = 'awk -f ../filter_empty_seqs.awk %s > %s' % (targetGenomes, filterGenomes)
os.system(awk)
fasta_write_test(filterGenomes)

#==============================================================================
#Code to automate the process with a list of genomes, and Elink from genomes
#to nucleotide. Turned out to be a pain in the ass.
#==============================================================================
#batchSize = 100000
#print('Esearch')
#
#for i in range(0, len(testTaxa)) :
#    esearch_term = '"Bacteria"[Organism] AND "%s"[Organism]' % testTaxa[i]
#    #print(esearch_term)
#    try :
#        xml = 'esearch_xml/%s_esearch.xml' % testTaxa[i]
#        outputHandle = open(xml, 'w')  
#        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
#        values = {'db':'genome','term':esearch_term,
#                  'rettype':'uilist','retmode':'uilist'}
#        data = urllib.parse.urlencode(values)
#        data = data.encode('utf-8')
#        Headers = {}
#        #Firefox 40.1 User Agent
#        Headers['User-Agent'] = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0)\
#                                 Gecko/20100101 Firefox/40.1'
#        #Request the URL with the variables
#        req = urllib.request.Request(url, data, headers = Headers)
#        #Response - visit the URL
#        resp = urllib.request.urlopen(req)
#        respData = resp.read().decode('utf-8')
#        respData = str(respData).replace('\n\n', '\n')
#        outputHandle.write(str(respData))
#        outputHandle.close()
#        
#    except Exception as e :
#        print(str(e))
##------------------------------------------------------------------------------
#ID = []
#for t in testTaxa :
#    for filename in os.listdir("./esearch_xml") :
#        if filename.startswith(t) :
#            with open('./esearch_xml/%s' % filename, 'r') as xml :
#                for line in xml :
#                    line = line.rstrip()
#                    gid = re.findall(r'<Id>(\d+)</Id>', line)
#                    ID.extend(gid)
#            xml.close()
#print(len(ID))
#
##------------------------------------------------------------------------------
##Elink
##------------------------------------------------------------------------------
#print('Elink')
#try :
#    xml = './elink_xml/%s_elink.xml' % taxon
#    outputHandle = open(xml, 'w')  
#    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
#    values = {'dbfrom':'nuccore','db':'nuccore',
#              'cmd':'neighbor_history','id':','.join(ID),
#              'rettype':'uilist','retmode':'uilist'}
#              
#    data = urllib.parse.urlencode(values)
#    data = data.encode('utf-8')
#    Headers = {}
#    #Firefox 40.1 User Agent
#    Headers['User-Agent'] = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0)\
#                             Gecko/20100101 Firefox/40.1'
#    #Request the URL with the variables
#    req = urllib.request.Request(url, data, headers = Headers)
#    #Response - visit the URL
#    resp = urllib.request.urlopen(req)
#    respData = resp.read().decode('utf-8')
#    respData = str(respData).replace('\n\n', '\n')
#    outputHandle.write(str(respData))
#    outputHandle.close()
#    
#except Exception as e :
#    print(str(e))
#
#queryKey = []
#webEnv = []
#
#with open('./elink_xml/%s_elink.xml' % taxon, 'r') as xml :
#    for line in xml :
#        line = line.strip()
#        if '<QueryKey>' in line :
#            qk = re.findall(r'<QueryKey>(\d+)</QueryKey>', line)
#            queryKey.extend(qk)
#        elif '<WebEnv>' in line :
#            we = re.findall(r'<WebEnv>(.*?)</WebEnv>', line)
#            webEnv.extend(we)
#        else :
#            pass
#        
#        
#xml.close()
#print(webEnv)
#print(queryKey)