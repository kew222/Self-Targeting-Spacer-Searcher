#This script will search for homologs of a protein and generate a list of nucleotide accession numbers 
#where they reside

#This script was written in 2017 by Kyle Watters in the Doudna Lab at UC Berkeley, Berkeley, CA
#Copyright (c) 2017 Kyle Watters. All rights reserved.

from __future__ import division
#from Bio.Seq import Seq
import getopt
import sys
import time
import httplib
from Bio import Entrez
from Bio.Blast import NCBIWWW,NCBIXML
from urllib2 import HTTPError  # for Python 2
Entrez.email = "watters@berkeley.edu"

help_message = '''
homology_locater_lite.py [-f file with multiple proteins | -s single protein name]

Options:
    -h, --help    Raise this help message
    -f, --file    Give a file with a list of protein accession numbers
    -s, --single  Give a single name of a protein accession number to search

*Note: only returns proteins found in Prokaryotes.        
'''
class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    def __init__(self):
        pass

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "f:hs:", ["file=", "help", "single="])
        except getopt.error, msg:
            raise Usage(msg)
        
        file_name = ''
        protein_to_search = ''
        
        for option, value in opts:
            if option in ("-f","--file"):
                file_name = value  
            if option in ("-s", "--single"):
                protein_to_search = value
            if option in ("-h", "--help"):
                raise Usage(help_message)  
                                                                    
        if len(args) != 0:
            raise Usage(help_message)    
                                    
        return args,file_name,protein_to_search

def truncate(f, n):
    #Truncates/pads a float f to n decimal places without rounding
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

def protein_blast(protein_to_search,E_value_limit=0.0001):

    attempt_num = 0
    while True:
        try:
            blastp2 = NCBIWWW.qblast("blastp", 'nr', protein_to_search, entrez_query='Prokaryote', expect=E_value_limit, hitlist_size=100000)#, format_type="Text")
            blast_record = NCBIXML.read(blastp2)
            break
        except httplib.IncompleteRead:  #If get an incomplete read, retry the request up to 3 times
            if attempt_num == 10:
                print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                return
            else:
                print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
            attempt_num += 1  
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt_num)
                attempt_num += 1
                time.sleep(1)
            else:
                raise 
        except RuntimeError:
            #NCBI probably closed the connection early, happens with poor internet connections
            if attempt_num == 3:
                print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                return
            else:
                print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
            attempt_num += 1
        except Exception as e:
                print('Unknown error: {0}'.format(e))
                attempt_num += 1
        time.sleep(5)
    return blast_record

def get_Nuc_uIDs(nr_blast_results,num_limit=1000000):
    
    acr_nr_genomes = []
    for result in nr_blast_results:
        protein_id = result[0]
        attempt_num = 0
        while True:
            try:
                handle = Entrez.elink(dbfrom='protein', db='nucleotide', id=protein_id, retmax=num_limit)
                record = Entrez.read(handle)
                handle.close()
                break
            except httplib.IncompleteRead:  #If get an incomplete read, retry the request up to 4 times
                if attempt_num == 10:
                    print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                    return
                else:
                    print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
                attempt_num += 1  
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt_num)
                    attempt_num += 1
                    time.sleep(1)
                else:
                    raise 
            except RuntimeError:
                #NCBI probably closed the connection early, happens with poor internet connections
                if attempt_num == 5:
                    print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                    return
                else:
                    print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
                attempt_num += 1
            except Exception as e:
                print('Unknown error: {0}'.format(e))
                attempt_num += 1
            time.sleep(5)
        
        holder = []
        for linked in record:
            try:
                for link in linked["LinkSetDb"][0]["Link"]:
                    holder.append(link['Id'])   
            except IndexError:
                print("No nucleotide links from genome ID {0}. Skipping...".format(protein_id))
                break
        acr_nr_genomes.append(holder)  
    
    return acr_nr_genomes

def fetch_nuc_accessions(acr_nr_genomes,nr_blast_results):
    
    protein_num = 0; complete_results = []; Accs_all = []
    for refgenomes in acr_nr_genomes:
        attempt_num = 0
        while True:
            try:
                handle2 = Entrez.efetch(db='nucleotide', rettype="acc", id=refgenomes)   #Get Acc#s of those found from search
                Accs = []
                for Id in handle2:
                    if Id.strip() != '':
                        Accs.append(Id.strip())       
                handle2.close()  
                break
            except httplib.IncompleteRead:  #If get an incomplete read, retry the request up to 3 times
                if attempt_num == 10:
                    print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                    return
                else:
                    print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
                attempt_num += 1  
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt_num)
                    attempt_num += 1
                    time.sleep(1)
                else:
                    raise 
            except RuntimeError:
                #NCBI probably closed the connection early, happens with poor internet connections
                if attempt_num == 3:
                    print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                    return
                else:
                    print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
                attempt_num += 1
            except Exception as e:
                print('Unknown error: {0}'.format(e))
                attempt_num += 1
                
        #Format of stored results: [protein acc, %i, query coverage, [list of found genomes] ] 
        complete_results.append([nr_blast_results[protein_num],Accs])
        Accs_all += Accs
    
        protein_num += 1
                  
    return complete_results, Accs_all
    
def get_assem_uIDs(acr_nr_genomes,num_limit=100000):
    
    #First, flatten the nucleotide lists
    all_nucleotides = []
    for genomes in acr_nr_genomes:
        all_nucleotides += genomes
    
    chunk_size = 100; Assem_uIDs = []
    for chunk in range(0,len(all_nucleotides),chunk_size):
        nucleotide_chunk = all_nucleotides[chunk:chunk+chunk_size]  
        attempt_num = 0
        while True:
            try:
                handle = Entrez.elink(dbfrom='nucleotide', db='assembly', id=nucleotide_chunk, retmax=num_limit)
                record = Entrez.read(handle)
                handle.close()
                break
            except httplib.IncompleteRead:  #If get an incomplete read, retry the request up to 3 times
                if attempt_num == 3:
                    print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                    return
                else:
                    print("httplib.IncompleteRead error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
                attempt_num += 1  
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt_num)
                    attempt_num += 1
                    time.sleep(1)
                else:
                    raise 
            except RuntimeError:
                #NCBI probably closed the connection early, happens with poor internet connections
                if attempt_num == 3:
                    print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Reached limit of {0} failed attempts.".format(attempt_num))
                    return
                else:
                    print("Runtime error at Entrez step linking assembly numbers to nucleotide database. Attempt #{0}. Retrying...".format(attempt_num))
                attempt_num += 1
            except Exception as e:
                print('Unknown error: {0}'.format(e))
                attempt_num += 1
                
        holder = []; nucleotide_num = 0
        for linked in record:
            try:
                for link in linked["LinkSetDb"][0]["Link"]:
                    holder.append(link['Id']) 
                    break  #Only take the first entry, most recent
            except IndexError:
                print("No assembly links from nucleotide uID {0}. Skipping...".format(nucleotide_chunk[nucleotide_num]))
                holder.append("")
            nucleotide_num += 1
            
        Assem_uIDs += holder
    
    #build a disctionary
    conversion_dict = {}; index = 0
    for el in all_nucleotides:
        conversion_dict[el] = Assem_uIDs[index]
        index += 1
              
    return conversion_dict
                                             

def find_homologs(protein_to_search):
    E_value_limit = 0.001; num_limit = 100000
    
    #Check the NR database first and give priority. Will check INSDC
    print("Blasting against NCBI nr protein database...")
    
    blast_record = protein_blast(protein_to_search,E_value_limit)
      
    nr_blast_results = []
    #Go through each record, get the %identity, query coverage, and name, use to find nucleotide records
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            protein_ident = str(alignment.title.split('|')[3])
            if "." in protein_ident:  #this will happen with pdb entries, etc. things that aren't versioned
                percent_ident = truncate(hsp.identities / hsp.align_length * 100,1)
                #NOTE: query coverage calculated represents the coverage over only the best single alignment (see break below)
                query_coverage = truncate(len(hsp.query.replace('-','')) / blast_record.query_length * 100,1)
                nr_blast_results.append([protein_ident,str(percent_ident)+'%',str(query_coverage)+'%'])  
            break  #only consider the best alignment
            
    #Take the protein results and link them to the nucleotide database
    acr_nr_genomes = get_Nuc_uIDs(nr_blast_results)
    
    #Get the nucleotide accession numbers to print out
    complete_results,Accs_all = fetch_nuc_accessions(acr_nr_genomes,nr_blast_results)
    
    #Make a dictionary between nucleotide uIDs and accessions
    collapsed_acr_nr_genomes = []
    for x in acr_nr_genomes:
        collapsed_acr_nr_genomes += x
    Nuc_dict = dict(zip(Accs_all,collapsed_acr_nr_genomes))
    
    #Make a dictionary with the assembly uIDs
    convert_N_to_A = get_assem_uIDs(acr_nr_genomes)
    
    #Store the results      
    with open("Protein_homologs_found.txt", "a") as save_file: 
        save_file.write("Searched protein: {0}\n".format(protein_to_search))
        save_file.write("Genome Acc\tAssembly uID\tProtein Acc\t%Identity\tQuery Coverage\n".format(protein_to_search))
        for result in complete_results:
            for genome in result[1]:
                assembly_uID = convert_N_to_A[Nuc_dict[genome]]
                save_file.write(genome + '\t' + assembly_uID + '\t'  +  '\t'.join(result[0]) + '\n')
        save_file.write('\n')

    return complete_results, convert_N_to_A, Nuc_dict

def main(argv=None):
    
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args,file_name,protein_to_search = params.parse_options(argv)
            
            if file_name != '':
                with open(file_name,'rU') as file1:
                    lines = file1.readlines()
                proteins_to_search = [line.strip() for line in lines]
            else:
                proteins_to_search = [protein_to_search]
                
            for protein in proteins_to_search:
                print("Looking up homologs for {0}...".format(protein))
                complete_results,convert_N_to_A,Nuc_dict = find_homologs(protein)
                
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())

