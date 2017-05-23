
import glob
import os
from anti_CRISPR_miner import import_data, locus_re_annotator,output_results


orig_dir = os.getcwd()+'/'

with open("genomes_to_fix.txt", "r") as file1:
    lines = file1.readlines()
    lines = [line.strip() for line in lines]

for species in lines:
    os.chdir(orig_dir + species + "/")
    print("\nCurrently in {0} directory.".format(species))
    spacer_files = glob.glob(orig_dir + species + "/Spacers_*")   #gets the list of spacers in that species
    
    for spacer_file in spacer_files:
        print("Analyzing {0} in {1}...".format(spacer_file.split("/")[-1],species))
        imported_data = import_data(spacer_file)
        re_analyzed_data = locus_re_annotator(imported_data)
        output_results(re_analyzed_data,{},{},spacer_file)   #Overwrites the old file with new data           

os.chdir(orig_dir)