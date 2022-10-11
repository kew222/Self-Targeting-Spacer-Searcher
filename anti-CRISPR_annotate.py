#Takes STSS output and annotates which genomes have anti-CRISPRs in them

from homology_locator_blastp import find_homologs
from all_Acrs_known import known_Acrs  #dictionary with protein: common name
import getopt, sys

help_message = '''
anti-CRISPR_annote.py is basically a wrapper for homology_locator_lite.py which takes a protein list input and an STSS output (with a leading
annotation for the assembly uID), looks up the proteins and adds which are present in a final column.

Usage: python anti-CRISPR_annote.py [-f file with multiple proteins | -s single protein name] -i <STSS input file>

Options:
    -h, --help    Raise this help message
    -f, --file    Give a file with a list of protein accession numbers (defaults to internal list updated in 2022)
    -s, --single  Give a single name of a protein accession number to search
    -i, --STSS    

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
            opts, args = getopt.getopt(argv[1:], "f:hs:i:", ["file=", "help", "single=", "STSS="])
        except getopt.error as msg:
            raise Usage(msg)
        
        file_name = ''
        protein_to_search = ''
        STSS_file = ''
        
        for option, value in opts:
            if option in ("-f","--file"):
                file_name = value  
            if option in ("-i","--STSS"):
                STSS_file = value
            if option in ("-s", "--single"):
                protein_to_search = value
            if option in ("-h", "--help"):
                raise Usage(help_message)  
                                                                    
        if len(args) != 0 or STSS_file == '':
            raise Usage(help_message)    
                                    
        return args,file_name,protein_to_search,STSS_file
                           
                
def main(argv=None):
    
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args,file_name,protein_to_search, STSS_file = params.parse_options(argv)
            
            #Load in the STSS data
            with open(STSS_file, 'r') as file1:
                STSS_data = file1.readlines()
            
            STSS_data = [x.strip()+'\t' for x in STSS_data]                              
            Assem_uIDs = [x.split('\t')[0] for x in STSS_data]
            
            #Run homolog_locator_lite.py to find all of the homologs
            if file_name != '':
                with open(file_name,'r') as file1:
                    lines = file1.readlines()
                proteins_to_search = [line.strip() for line in lines]
            elif protein_to_search != "":
                proteins_to_search = [protein_to_search]
            else:
                proteins_to_search = [key for key,value in known_Acrs.iteritems()] #search list in all_Acrs.py by default
                         
            for protein in proteins_to_search:
                print("Looking up homologs for {0}...".format(protein))
                complete_results,convert_N_to_A, Nuc_dict = find_homologs(protein)
                
                #annotate the STSS results
                for result in complete_results:
                    for genome in result[1]:
                        assembly_uID = convert_N_to_A[Nuc_dict[genome]]
                        protein_found = result[0][0]
                        if assembly_uID in Assem_uIDs:
                            indexes = [index for index, value in enumerate(Assem_uIDs) if value == assembly_uID]
                            for index in indexes:
                                old_data = STSS_data[index]
                                if old_data[-1] != "\t":
                                    old_data += ', '  #add a comma for clarity
                                STSS_data = STSS_data[:index] + [old_data + "{0} ({1})".format(known_Acrs[protein],protein_found)] + STSS_data[index+1:]
                
            #Add a indicator for lines that had no anti-CRISPR
            index = 0
            for line in STSS_data:
                old_data = STSS_data[index]
                if old_data[-1] == "\t":
                    old_data += 'None'  
                    STSS_data = STSS_data[:index] + [old_data] + STSS_data[index+1:]
                index += 1
                
            with open('anti_CRISPRs_annotations.txt','w') as file1:
                for el in STSS_data:
                    file1.write(el + '\n')
            
                
    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
        return 2

if __name__ == "__main__":
    sys.exit(main())                                                              
                                                                                                          
                                                                                                                                                                                                  
                               