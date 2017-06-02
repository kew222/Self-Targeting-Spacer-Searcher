#!/usr/bin/env python

#Wrapper for STSS.py. This code takes a list of bacterial genomes and searches through each group (in group mode) or a list of individual 
#genomes (in update mode) using the STSS.

#This script was written in 2016 by Kyle Watters
#Copyright (c) 2016 Kyle Watters. All rights reserved.

from anti_CRISPR_miner import anti_CRISPR_search,link_nucleotide_to_bioproject,link_assembly_to_nucleotide,gather_assemblies_from_bioproject_IDs,download_genomes,get_Accs
import getopt
import sys
import subprocess
from datetime import datetime
import os

help_message = '''
STSS_wrapper.py uses anti_CRISPR_miner.py to search for all of the genomes in a list provided
by the user. It can function by group or just run a list of genomes. It does not support protein clustering.

Usage:
   STSS.py [options] [--group <groups_file> | --update <IDs_file> ]

Options
-h, --help                      Opens help message
-v, --version                   Displays version number
-g, --group  <groups_file>      Search by group. Makes a directory for each group provided and searches it on NCBI
-u, --update <IDs_file>         Takes a list of genomes, searches for all those that were analyzed, and runs those that weren't
                                (intended to be used to update self-targeting list as genomes are uploaded to NCBI;
                                 Currently, only supports use of bioproject IDs, available on NCBI ftp server)
--auto-update                   Automatically pulls the updates from NCBI (don't need assemblies file)
-e, --exclude <filename>        Accepts a list of genomes to exclude from the update function                

-f, --force-redownload          Forces redownloading of genomes that were already downloaded
-l, --limit <N>                 Limit Entrez search to the first N results found (default: 1000)
--complete-only                 Only return complete genomes from NCBI
-E, --E-value  <N>              Upper limit of E-value to accept BLAST results as protein family (default: 1e-3)
--percent-reject <N>            Percentage (of 100%) to use a cutoff from the average spacer length to reject validity of a proposed locus (default 40)
--all-islands                   Include all unknown proteins found within a predicted MGE
--outside-islands               Include proteins from predicted MGEs when the spacer is outside a predicted MGE island
                                (automatically includes --all-islands option)
-s, --spacers <N>               Number of spacers needed to declare a CRISPR locus (default: 3)
--pad-locus <N>                 Include a buffer around a potential locus to prevent missed repeats from
                                appearing as hits (default: 100)
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
            opts, args = getopt.getopt(argv[1:], "hvg:u:E:l:s:fe:",
                                       ["limit=",
                                       "help",
                                       "group=",
                                       "update=",
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
                                       "auto-update",
                                       "skip-PHASTER",
                                       "rerun-PHASTER=",
                                       "align-families",
                                       "exclude=",
                                       "version"])
        except getopt.error, msg:
            raise Usage(msg)
        
        group = False
        update = False
        update_group_file = ''
        search = ''
        provided_dir = ''
        Accs_input = ''
        default_limit = 100000
        num_limit = 0
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
        auto_update = False
        
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
            if option in ("-g","--group"):
                group = True
                update_group_file = value  
            if option in ("-u","--update"):
                update = True 
                update_group_file = value
            if option == "--auto-update":
                auto_update = True  
                update = True 
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
                                                                                                                                                                    
        if len(args) != 0:
            if not auto_update:
                print("Wrong number of arguments, received {0}, should be none.".format(len(args)))
                raise Usage(help_message)   
            
        if group and update:
            print("Please select either --group or --update")  
        elif not group and not update:
            print("Please select either --group or --update") 
        
        return args,num_limit,E_value_limit,provided_dir,search,all_islands,in_islands_only,repeats,skip_family_search,families_limit,pad_locus,skip_family_create,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,rerun_PHASTER,PHASTER_input_file,skip_alignment,group,update,exclude_file,update_group_file,Accs_input,auto_update
    
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

def update_genomes_list(filename="NCBI_prokaryotes.txt"):
    
    #Get list from their ftp server
    ftp_cmd = "curl ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt -o {0}".format(filename)
    update_genomes = subprocess.Popen(ftp_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output,error = update_genomes.communicate() 
    #if error != '':
    #    print([error])
    #    sys.exit()
    with open(filename, "rU") as filetemp:
        lines = filetemp.readlines()
    bioprojectIDs = []
    for line in lines:
        splitit = line.split("\t")
        if splitit[3].isdigit():
            bioprojectIDs.append(splitit[3])
    return bioprojectIDs
                    
def main(argv=None):
    
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args,num_limit,E_value_limit,provided_dir,search,all_islands,in_islands_only,repeats,skip_family_search,families_limit,pad_locus,skip_family_create,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,rerun_PHASTER,PHASTER_input_file,skip_alignment,group,update,exclude_file,update_group_file,Accs_input,auto_update = params.parse_options(argv)
            params.check()
        
        orig_dir = os.getcwd()+'/'
        
        if group:
            groups_to_search = import_list(update_group_file)
            for group_to_search in groups_to_search:
                groups_remaining = groups_to_search[groups_to_search.index(group_to_search):]
                rescue_list(groups_remaining)
                dir_name = group_to_search.replace(" ","_")
                if not os.path.exists(dir_name):
                    os.mkdir(dir_name)
                try:
                    os.chdir(dir_name)
                    print("\nCurrently in {0} directory.".format(dir_name))
                except:
                    print("Cannot get into {0}. Exiting...".format(dir_name))
                    sys.exit()
                provided_dir = ''  #Force searching with the wrapper
                print("Searching '{0}'...\n{1}\n".format(group_to_search,str(datetime.now())))
                
                #Identify genomes that contain self-targeting spacers     
                current_dir=orig_dir+dir_name+"/"
                assemblies = []
                try:
                    protein_list = anti_CRISPR_search(provided_dir,group_to_search,num_limit,E_value_limit,all_islands,in_islands_only,repeats,skip_family_search,families_limit,pad_locus,skip_family_create,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,current_dir)
                    protein_list = [] #Currently not used, clear memory
                except SystemExit:
                    pass #Presumedly, this is raised when the anti_CRISPR_miner function finds that further analysis in that group isn't necessary
                os.chdir(orig_dir)
        
        elif update:
            exclude_genomes = []
            if exclude_file != "":
                print("Loading genomes to exclude")
                with open(exclude_file, 'rU') as file1:
                    exclude_genomes = [x.strip() for x in file1.readlines()]
            
            #Get a list of bioprojects associated with these nucleotide entries
     #       bioproject_chunk_size = 50
            excluded_bioprojects = []
            #Loop over chunks in case they are large
     #       for chunk in range(0,len(exclude_genomes),bioproject_chunk_size):
     ##           if chunk + bioproject_chunk_size < len(exclude_genomes):
      #              bioprojects_chunk = exclude_genomes[chunk:chunk+bioproject_chunk_size]
      #          else: 
      #              bioprojects_chunk = exclude_genomes[chunk:len(exclude_genomes)]
      #          excluded_bioprojects += link_nucleotide_to_bioproject(bioprojects_chunk,[])   #results in a list of assembly uids to exclude
      #      exclude_genomes = []
            
      #      #Print out a list of the bioprojects to exclude
      #      with open("Bioproject_uids_to_exclude.txt", "w") as file1:
      #          for bioproject in excluded_bioprojects:
      #              file1.write(bioproject + "\n")       
            
            if auto_update:
                #Query NCBI for a list of genomes to update
                print("Fetching all prokaryotic genomes from NCBI...")
                bioprojectIDs = update_genomes_list()
            else:
                #Import the assemblies
                print("Importing Bioproject uIDs to run...")
                bioprojectIDs = import_list(update_group_file)
            
            text1 = "Finding assemblies for bioprojects..."
            text2 = " and filtering out previously run genomes..."
            if excluded_bioprojects == []:
                print(text1)
            else:
                print(text1[:-3]+text2)
            
            #Remove excluded_assemblies form all_assemblies
            filtered_bioprojects = []
            for ID in bioprojectIDs:
                if ID not in excluded_bioprojects:
                    filtered_bioprojects.append(ID)
            excluded_bioprojects = []  #clear memory
            
            #Link bioproject IDs to assemblies, but breakup into groups of 100
            query_chunk = 50
            assemblies = []
            for piece in range(0, len(filtered_bioprojects), query_chunk):
                if piece + query_chunk < len(filtered_bioprojects):
                    IDs_to_send = filtered_bioprojects[piece:piece+query_chunk]
                else:
                    IDs_to_send = filtered_bioprojects[piece:len(bioprojectIDs)]
                assemblies += gather_assemblies_from_bioproject_IDs(IDs_to_send)   
            
            #Print out a list of the assemblies to downloaded
            with open("Assembly_uids_found.txt", "w") as file1:
                for assembly in assemblies:
                    file1.write(assembly + "\n")       
            
            #Remove duplicate assembly numbers
            filtered_assemblies = list(set(assemblies))                    
                                                                           
            #Link assemblies to nuccore
            found_complete = 0; found_WGS = 0; total = 0; complete_IDs = []; WGS_IDs = []; wgs_master_GIs = []; num_genomes = 0; num_limit=100000; complete_only=False
            found_complete,found_WGS,total,complete_IDs,WGS_IDs,wgs_master_GIs,num_genomes = link_assembly_to_nucleotide(filtered_assemblies,num_limit,complete_only,num_genomes,complete_IDs,WGS_IDs,wgs_master_GIs)  
            
            #Make a new directory to run in
            if not os.path.exists("Genomes_Update"):
                os.mkdir("Genomes_Update")
            try:
                os.chdir("Genomes_Update")
                current_dir=orig_dir+"Genomes_Update/"
            except:
                print("Cannot get into {0}. Exiting...".format("Genomes_Update/"))
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
            with open("Genomes_to_download.txt", "w") as file1:
                for genome in (complete_IDs+wgs_master_GIs):
                    file1.write(genome + "\n")       
            
            num_limit = 10000; search = ""; redownload=False; provided_dir=""; fastanames={}  #variables required for download, set to defaults
            fastanames,Acc_convert_to_GI = download_genomes(total,num_limit,num_genomes,found_complete,search,redownload,provided_dir,current_dir,found_WGS,complete_IDs,WGS_IDs,wgs_master_GIs,fastanames)
            
            #Use the provided genomes pathway to run anti_CRISPR search
            provided_dir = "downloaded_genomes/"
            try:
                protein_list = anti_CRISPR_search(provided_dir,search,num_limit,E_value_limit,all_islands,in_islands_only,repeats,skip_family_search,families_limit,pad_locus,skip_family_create,complete_only,skip_PHASTER,percent_reject,default_limit,redownload,current_dir)
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
