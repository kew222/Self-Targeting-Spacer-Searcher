#Get a list of Accession #s from a list of Nucleotide Accession #s

from __future__ import division
from Bio import Entrez
from user_email import email_address
Entrez.email = email_address
import httplib
import getopt, sys

help_message = '''
python nucleotide_acc_to_assembly.py <file with list of Nuc Acc #s>  '''

def get_version():
    return "0.0.1"
    
class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    
    def __init__(self):
        pass
    
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "vh",["help","version"])
        
        except getopt.error, msg:
            raise Usage(msg)
         
        for option, value in opts:
            if option in ("-v", "--version"):
                print "python nucleotide_acc_to_assembly.py v%s" % (get_version())
                exit(0)
            
        if len(args) != 1:
            raise Usage(help_message)    
                                                                                                                             
        return args

def get_Accs(IDs,database="assembly"):
    
    attempt = 1
    while attempt < 4:
        Accs = []
        handle = Entrez.efetch(db=database, rettype="acc", id=IDs)   #Get Acc#s of those found from search
        try:    
            handle = Entrez.efetch(db=database, rettype="acc", id=IDs)   #Get Acc#s of those found from search
            for Id in handle:
                if Id.strip() != "":
                    Accs.append(Id.strip())          
            handle.close()
            break
        except:
            if attempt < 4:
                attempt += 1
            else:
                raise
        
    return Accs

def link_to_assembly(nucleotide_list,database1='nucleotide',database2='assembly',assemblies=[],num_limit=100000):
    
    attempt_num = 1
    while True:
        try:
            handle4 = Entrez.elink(dbfrom=database1, db=database2, id=nucleotide_list, retmax=num_limit)
            record4 = Entrez.read(handle4)
            handle4.close()
            break
        except httplib.IncompleteRead:  #If get an incomplete read, retry the request up to 3 times
            if attempt_num == 3:
                print("httplib.IncompleteRead error at Entrez step linking {0} to {1} database. Reached limit of {2} failed attempts.".format(database1,database2,attempt_num))
                return
            else:
                print("httplib.IncompleteRead error at Entrez step linking {0} to {1} database. Attempt #{2}. Retrying...".format(database1,database2,attempt_num))
            attempt_num += 1
        except IndexError as e: 
            print("IndexError: ", e)
            if attempt_num == 3:
                print("Encountered problems linking {0} to {1} database. Reached limit of {2} failed attempts.".format(database1,database2,attempt_num))
                return         
        except RuntimeError:
                #NCBI probably closed the connection early, happens with poor internet connections
                if attempt_num == 3:
                    print("Runtime error at Entrez step linking {0} to {1} database. Reached limit of {2} failed attempts.".format(database1,database2,attempt_num))
                    return
                else:
                    print("Runtime error at Entrez step linking {0} to {1} database. Attempt #{2}. Retrying...".format(database1,database2,attempt_num))
                attempt_num += 1
        except Exception as e:
            print('Unknown Error: {0}'.format(e))
            attempt_num += 1  
    
    nucleotide_num = 0
    for linked in record4:
        try:
            for link in linked["LinkSetDb"][0]["Link"]:
                assemblies.append(link['Id']) 
                break  
        except IndexError:
            print("No {0} links from {1} ID {2}. Skipping...".format(database2,database1,nucleotide_list[nucleotide_num]))
            assemblies.append("")
        nucleotide_num += 1 
    
    return assemblies   

def NCBI_search(search,database,num_limit=100000,tag="",exclude_term=""):
    
    attempt_num = 1
    while True:
        try:
            search_term = search + tag + exclude_term   #otherwise, you get isolated viruses
            handle = Entrez.esearch(db=database,term=search_term, retmax=num_limit)
            record = Entrez.read(handle)
            handle.close()
            break
        except:  #If get an incomplete read, retry the request up to 3 times
            if attempt_num == 5:
                print("httplib.IncompleteRead error at Entrez genome search step. Reached limit of {0} failed attempts.".format(attempt_num))
                return
            else:
                print("httplib.IncompleteRead error at Entrez genome search step. #{0}. Retrying...".format(attempt_num))
            attempt_num += 1
    genomes = record["IdList"]
    return genomes
      

                        
def main(argv=None):
    
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            
            with open(args[0], 'rU') as file1:
                lines = file1.readlines()
            Nuc_Accs = [x.strip() for x in lines]
        else:
            if type(argv) != list:
                print('Must pass a list of nucleotide acessions to nucleotide_acc_to_assembly_acc')
            Nuc_Accs = argv
        
        #Shrink down the size of the list for speed:
        compressed_Nuc_Accs = []
        for acc in Nuc_Accs:
            if acc not in compressed_Nuc_Accs:
                compressed_Nuc_Accs.append(acc)
        
        #Convert with chunks at a time
        chunk_size = 100; Assem_uIDs_all = []; activity_ticker = 1; Nuc_uIDs = []
        for chunk in range(0,len(compressed_Nuc_Accs),chunk_size):
            print("Getting Accessions for batch {0} of {1}....".format(activity_ticker,-(-len(compressed_Nuc_Accs)//chunk_size)))
            starting_pos = chunk; recursive_chunk = len(compressed_Nuc_Accs[chunk:chunk+chunk_size])
            while True:  #breaks down the chunks further if some of the search terms return nothing to figure out which are the problem, as they are not indexed
                nucleotide_chunk = compressed_Nuc_Accs[starting_pos:starting_pos+recursive_chunk]  
                Nuc_uIDs_temp = NCBI_search(", ".join(nucleotide_chunk),"nucleotide")
                if len(Nuc_uIDs_temp) == recursive_chunk:
                    Nuc_uIDs += Nuc_uIDs_temp
                    starting_pos = starting_pos+recursive_chunk
                    recursive_chunk = chunk_size - (starting_pos - chunk)
                elif len(Nuc_uIDs_temp) == 0:
                    Nuc_uIDs += [""] * recursive_chunk  #add blanks for each position that's coming up as no good, skip over
                    print('problems with {0}'.format(nucleotide_chunk))
                    starting_pos = starting_pos+recursive_chunk
                    recursive_chunk = chunk_size - (starting_pos - chunk)
                else:
                    print('found a bad nucleotide accession number, splitting up chunk...trying {0} to {1}...'.format(starting_pos,starting_pos+recursive_chunk))
                    recursive_chunk = recursive_chunk // 2
                if len(Nuc_uIDs) == activity_ticker * chunk_size:
                    break
                elif chunk+chunk_size > len(compressed_Nuc_Accs[chunk:chunk+chunk_size]) and len(Nuc_uIDs) == len(compressed_Nuc_Accs):
                    break
            activity_ticker += 1
            
        Nuc_dict = dict(zip(compressed_Nuc_Accs,Nuc_uIDs))
        compressed_Nuc_uIDs = []
        for Nuc in Nuc_uIDs:
            if Nuc != "":
                compressed_Nuc_uIDs.append(Nuc)
        
        chunk_size = 100; Assem_uIDs_all = []
        for chunk in range(0,len(compressed_Nuc_uIDs),chunk_size):
            nucleotide_chunk = compressed_Nuc_uIDs[chunk:chunk+chunk_size]  
            Assem_uIDs = link_to_assembly(nucleotide_chunk,'nucleotide','assembly',[])
            Assem_uIDs_all += Assem_uIDs
            
        Assem_dict = dict(zip(compressed_Nuc_uIDs,Assem_uIDs_all))    
            
        #build a dictionary:
        convert_N_to_A = {}
        for acc in compressed_Nuc_Accs:
            key = Nuc_dict[acc]
            if key != "":
                convert_N_to_A[acc] = Assem_dict[Nuc_dict[acc]]
            else:
                convert_N_to_A[acc] = "-"       
            
        #Convert the original uncompressed list of Nucleotide accessions to assemblies:
        Assem_uIDs = []
        for acc in Nuc_Accs:
            Assem_uIDs.append(convert_N_to_A[acc])

        if __name__ == "__main__":
            with open('Assembly_Accs.txt','w') as file1:
                for Acc in Assem_uIDs:
                    file1.write(Acc + '\n')
        else:
            return Assem_uIDs
                       
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())

             