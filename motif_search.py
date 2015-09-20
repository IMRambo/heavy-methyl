#Ian Rambo
#University of Delaware
#BINF 868

#==============================================================================
'''
GOAL: find motifs 
'''
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
        
        

os.chdir('/Users/imrambo/Documents/BINF868')

sequences = []
descr = None

with open('Ypestis_coding.fa', 'r') as yp :
    line = yp
    

    
            