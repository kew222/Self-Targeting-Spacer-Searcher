#!/usr/bin/env python

#Wrapper for STSS.py. This code takes a list of bacterial genomes and searches through each group (in group mode) or a list of individual 
#genomes (in update mode) using the STSS.

#This script was written in 2016 by Kyle Watters
#Copyright (c) 2016 Kyle Watters. All rights reserved.

from STSS import self_target_search,NCBI_search,link_assembly_to_nucleotide,link_nucleotide_to_assembly,download_genomes,get_Accs
import getopt
import sys
import os

bin_path = os.path.dirname(os.path.realpath(__file__)) + "/bin/"
HMM_dir = os.path.dirname(os.path.realpath(__file__)) + "/HMMs/"

help_message = '''
STSS_auto_update.py uses STSS.py to take a search string and grab all of the genomes from NCBI with that tag. 
It can also accept a list of nucleotide genomes to filter out genomes that were already run.
This script is really meant to use both, if no genomes to exclude, STSS.py is recommended to use instead.
Remember quotes if the search term is multiple words, do not include tags like [organism] - this is automatic.

Usage:
   STSS_auto_update.py [options] <search string>

Options
-h, --help                      Opens help message
-v, --version                   Displays version number

-e, --exclude <filename>        Accepts a list of genomes to exclude from the update function                

-o, --prefix <string>           Prefix for filenames in output (ex. prefix of 'bagel' gives: bagel_Spacers...islands.txt)
-f, --force-redownload          Forces redownloading of genomes that were already downloaded
-l, --limit <N>                 Limit Entrez search to the first N results found (default: 1000)
--CDD                           Use the Conserved Domain Database to identify Cas proteins (default is to use HMMs)
--complete-only                 Only return complete genomes from NCBI
-E, --E-value  <N>              Upper limit of E-value to accept BLAST results as protein family (default: 1e-3)
--percent-reject <N>            Percentage (of 100%) to use a cutoff from the average spacer length to reject validity of a proposed locus (default 40)

--all-islands                   Include all unknown proteins found within a predicted MGE
--outside-islands               Include proteins from predicted MGEs when the spacer is outside a predicted MGE island
                                (automatically includes --all-islands option)
-s, --spacers <N>               Number of spacers needed to declare a CRISPR locus (default: 3)
--pad-locus <N>                 Include a buffer around a potential locus to prevent missed repeats from
                                appearing as hits (default: 100)
-d, --Cas-gene-distance <N>     Window around an array to search for Cas proteins to determine CRISPR subtype
                                (default: 20000 - input 0 to search whole genome)                                 
--skip-PHASTER                  Skip PHASTER analysis (currently can't upload search files)                                
'''


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
            opts, args = getopt.getopt(argv[1:], "hvE:l:s:fe:o:d:",
                                       ["limit=",
                                       "help",
                                       "all-islands",
                                       "outside-islands",
                                       "E-value=",
                                       "spacers=",
                                       "NCBI",
                                       "skip-family-search",
                                       "families-limit=",
                                       "cluster-search",
                                       "pad-locus=",
                                       "force-redownload",
                                       "percent-reject=",
                                       "complete-only",
                                       "skip-PHASTER",
                                       "rerun-PHASTER=",
                                       "align-families",
                                       "exclude=",
                                       "version"
                                       "CDD",
                                       "prefix="
                                       "Cas-gene-distance="])
        except getopt.error, msg:
            raise Usage(msg)
        
        update_group_file = ''
        search = ''
        provided_dir = ''
        Accs_input = ''
        default_limit = 1000000
        num_limit = 1000000
        E_value_limit = 1e-3
        all_islands = False
        in_islands_only = True
        repeats = 4  #zero indexed (actually 3)
        skip_family_search = False
        families_limit = 300
        pad_locus = 40
        skip_family_create = True
        complete_only = False
        skip_PHASTER = False
        rerun_PHASTER = False
        skip_alignment = True
        PHASTER_input_file = ''
        percent_reject = 40
        redownload = False
        exclude_file = ""
        CDD = False
        Cas_gene_distance = 20000
        prefix = ''
        
        for option, value in opts:
            if option in ("-v", "--version"):
                print "anti_CRISPR_miner.py v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)  
            if option in ("-l","--limit"):
                num_limit = int(value)   
            if option in ("-E", "--E-value"):
                E_value_limit = float(value)     
            if option == "--all-islands":
                all_islands = True
            if option == "--outside-islands":
                in_islands_only = False  
                all_islands = True
            if option in ("-s","--spacers"):
                repeats = int(value) + 1 
            if option == "--pad-locus":
                pad_locus = int(value) 
            if option in ('c',"--cluster-search"):
                skip_family_create = False 
            if option == "--complete-only":
                complete_only = True   
            if option == "--skip-PHASTER":
                skip_PHASTER = True  
            if option in ("-f","--force-redownload"):
                redownload = True   
            if option in ("-e","--exclude"):
                exclude_file = value      
            if option == "--CDD":
                CDD = True
            if option in ('-d',"--Cas-gene-distance"):
                Cas_gene_distance = int(value) 
            if option in ("-o","--prefix"):
                prefix = value + "_"     
                                                                                                                                                                    
        if len(args) != 1:
            print("Wrong number of arguments, received {0}, should be exactly 1: the search string.".format(len(args)))
            raise Usage(help_message)   
            
        return args,num_limit,E_value_limit,provided_dir,search,all_islands,in_islands_only,repeats,skip_family_search,families_limit,pad_locus,skip_family_create,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,rerun_PHASTER,PHASTER_input_file,skip_alignment,exclude_file,update_group_file,Accs_input,CDD,Cas_gene_distance,prefix
    
    def check(self):
        pass

def import_list(update_group_file):
    
    with open(update_group_file, 'r') as file1:
        groups_to_search = [x.strip() for x in file1.readlines()]
    
    return groups_to_search

def rescue_list(groups_remaining):  #Export the remaining to search so you can continue after a break (intended or unintended)
    with open("groups_remaining_export.log", 'w') as file2:  
        for group in groups_remaining:
            file2.write(group + "\n")

                    
def main(argv=None):
    
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args,num_limit,E_value_limit,provided_dir,search,all_islands,in_islands_only,repeats,skip_family_search,families_limit,pad_locus,skip_family_create,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,rerun_PHASTER,PHASTER_input_file,skip_alignment,exclude_file,update_group_file,Accs_input,CDD,Cas_gene_distance,prefix = params.parse_options(argv)
            params.check()
        
        search_string = args[0]
        orig_dir = os.getcwd()+'/'
        
        exclude_genomes = []; excluded_assemblies = []
        if exclude_file != "":
            print("Loading genomes to exclude...")
            with open(exclude_file, 'rU') as file1:
                exclude_genomes = [x.strip() for x in file1.readlines()]
            
            #Find the nucleotide uIDs for the excluded genomes
            query_chunk = 90
            num_chunks = len(range(0, len(exclude_genomes), query_chunk))
            i = 1
            for piece in range(0, len(exclude_genomes), query_chunk):
                print("Looking up uIDs for excluded genomes part {0} of {1}...".format(i,num_chunks))
                if piece + query_chunk < len(exclude_genomes):
                    search_list = exclude_genomes[piece:piece+query_chunk]
                else:                
                    search_list = exclude_genomes[piece:len(exclude_genomes)]
                found_nucleotide = NCBI_search(", ".join(search_list),"nucleotide",num_limit=100000,tag="",exclude_term="")
                found_assemblies = link_nucleotide_to_assembly(found_nucleotide,assemblies=[],num_limit=100000)
                excluded_assemblies += found_assemblies
                i+= 1
                                
        #Query NCBI for a complete list of genomes from the provided search string
        print("Finding list of assemblies from search of '{0}'...".format(search_string))
        all_assemblies = NCBI_search(search_string,"assembly",num_limit=1000000,tag="[organism]",exclude_term=" NOT phage NOT virus")
         
        print('Found {0} assemblies from search.'.format(len(all_assemblies)))    
        
        if all_assemblies == []:
            print("No genomes found needing to be run. Exiting...")
            exit()
            
        #filter out the excluded genomes
        filtered_assemblies = []
        for ID in all_assemblies:
            if ID not in excluded_assemblies:
                filtered_assemblies.append(ID)
        excluded_assemblies = []  #clear memory
        
        #Print out a list of the assemblies to be downloaded and excluded
        with open("{0}assembly_uids_to_be_downloaded.txt".format(prefix), "w") as file1:
            for assembly in filtered_assemblies:
                file1.write(assembly + "\n")       
        with open("{0}excluded_assembly_uids_found.txt".format(prefix), "w") as file1:
            for assembly in filtered_assemblies:
                file1.write(assembly + "\n")       
    
        #Remove duplicate assembly numbers
        filtered_assemblies = list(set(filtered_assemblies))                    
        print("{0} of {0} genomes found need to be downloaded.".format(len(filtered_assemblies),len(all_assemblies)))
        
        #Link assemblies to nuccore
        print("Getting nucleotide uIDs for assemblies to update....")
        found_complete = 0; found_WGS = 0; total = 0; complete_IDs = []; WGS_IDs = []; wgs_master_GIs = []; num_genomes = 0
        found_complete,found_WGS,total,complete_IDs,WGS_IDs,wgs_master_GIs,num_genomes = link_assembly_to_nucleotide(filtered_assemblies,num_limit,complete_only,num_genomes,complete_IDs,WGS_IDs,wgs_master_GIs,False)   
        
        #Make a new directory to run in
        if not os.path.exists("{0}Genomes_Update".format(prefix)):
            os.mkdir("{0}Genomes_Update".format(prefix))
        try:
            os.chdir("{0}Genomes_Update".format(prefix))
            current_dir=orig_dir+"{0}Genomes_Update/".format(prefix)
        except:
            print("Cannot get into {0}. Exiting...".format("{0}Genomes_Update/".format(prefix)))
            sys.exit()
        
        #Perform a filter on the nucleotide records to remove potential duplicates
        #Get accession numbers from UIDs
        batch_size = 100 ; Accs = []
        all_IDs = complete_IDs+wgs_master_GIs
        for start in range(0, len(all_IDs), batch_size):
            end = min(len(all_IDs), start+batch_size)
            IDs = all_IDs[start:end]
            Accs += get_Accs(IDs)
        Acc_convert_to_GI = dict(zip(Accs,all_IDs))
        
        #Check list against those already run and remove duplicates
        GIs_to_remove = []
        for Acc in Accs:
            if Acc in exclude_genomes:
                GIs_to_remove.append(Acc_convert_to_GI[Acc])    
        for ID in complete_IDs:
            if ID in GIs_to_remove:
                index = complete_IDs.index(ID)
                complete_IDs.pop(index)
        for ID in wgs_master_GIs:
            if ID in GIs_to_remove:
                index = wgs_master_GIs.index(ID)
                wgs_master_GIs.pop(index) 
                WGS_IDs.pop(index)           
        
        #Update the sums
        found_complete = len(complete_IDs)
        found_WGS = len(wgs_master_GIs)
        total = found_complete + found_WGS            
        
        #Print out a list of the genomes to downloaded
        with open("{0}Genomes_to_download.txt".format(prefix), "w") as file1:
            for genome in (complete_IDs+wgs_master_GIs):
                file1.write(genome + "\n")       
        
        num_limit = 10000; search = ""; redownload=False; provided_dir=""; fastanames={}  #variables required for download, set to defaults
        fastanames,Acc_convert_to_GI = download_genomes(total,num_limit,num_genomes,found_complete,search,redownload,provided_dir,current_dir,found_WGS,complete_IDs,WGS_IDs,wgs_master_GIs,fastanames,False,prefix)
        
        #Use the provided genomes pathway to run self-targeting search
        provided_dir = "{0}downloaded_genomes/".format(prefix)
        try:
            protein_list = self_target_search(provided_dir,search,num_limit,E_value_limit,repeats,pad_locus,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,current_dir,bin_path,Cas_gene_distance,HMM_dir,prefix,CDD,False)
            protein_list = [] #Currently not used, clear memory
        except SystemExit:
            pass #Presumedly, this is raised when the anti_CRISPR_miner function finds that further analysis in that group isn't necessary
        os.chdir(orig_dir)                                                                                                       
                                                                                          
                                     
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())
