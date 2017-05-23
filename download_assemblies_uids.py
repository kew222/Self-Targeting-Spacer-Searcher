from anti_CRISPR_miner import link_assembly_to_nucleotide,download_genomes
import os, sys

orig_dir = os.getcwd()+'/'

#import assemblies UIds
with open("sorted_assemblies.txt", "rU") as file1:
    lines = file1.readlines()

assemblies = [x.strip() for x in lines]

filtered_assemblies = list(set(assemblies))                     
                                                                           
#Link assemblies to nuccore
found_complete = 0; found_WGS = 0; total = 0; complete_IDs = []; WGS_IDs = []; wgs_master_GIs = []; num_genomes = 0; num_limit=150000; complete_only=False
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

num_limit = 15000; search = ""; redownload=False; provided_dir=""; fastanames={}  #variables required for download, set to defaults
fastanames,Acc_convert_to_GI = download_genomes(total,num_limit,num_genomes,found_complete,search,redownload,provided_dir,current_dir,found_WGS,complete_IDs,WGS_IDs,wgs_master_GIs,fastanames)
