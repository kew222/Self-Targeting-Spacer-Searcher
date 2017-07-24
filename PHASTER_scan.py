#!/usr/bin/env python
# -*- coding: utf-8 -*-

#PHASTER_scan.py Checks through all of the contigs of a genome and looks for phage islands with PHASTER.

#This script was written in 2017by Kyle Watters
#Copyright (c) 2017 Kyle Watters. All rights reserved.

#Dependencies include:  BLAST (local), biopython, CRT CRISPR finder tool, Clustal Omega (and thus argtable2)

import getopt
import sys, os
from STSS import link_nucleotide_to_assembly, link_assembly_to_nucleotide, query_PHASTER, get_Accs
from Bio import Entrez

Entrez.email = "watters@berkeley.edu"

help_message = '''
PHASTER_scan.py takes a nucleotide accession number, looks up the contigs from the assembly and searches each contig for prophages with PHASTER

Usage:
   PHASTER_scan.py [options] <Accession number or text file with list>

Options
-h, --help                      Opens help message
-v, --version                   Displays version number
-l, list                        Input is a text file list (default: False)

'''

def get_version():
    return "0.1.0"

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    
    def __init__(self):
        pass
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvl",
                                       ["help",
                                       "version",
                                       "list"])
        except getopt.error, msg:
            raise Usage(msg)
        
        list_input = False
        for option, value in opts:
            if option in ("--list", "-l"):
                list_input = True
                                                                           
        if len(args) != 1:
            raise Usage(help_message)    
                                    
        return args, list_input

def parse_PHASTER(lines,Acc_to_search,results=[]):
    
    #expression1 = re.compile("(Totally\s+)(\d+)(\s+intact prophage regions have been identified)")
    record = False
    for line in lines:
        if not record:
            #a = expression1.match(line)
            #if a is not None:
            #    if int(a.group(2)) == 0: #there are no prophages
            #        break
            if line.strip()[:20] == "-"*20:
                record = True
        else:
            if line != "" and line != "\n":
                data = line.strip().split()
                results.append([Acc_to_search] + [data[0]] + [data[1]] + [data[6]] + [data[4]])  
    return results
             
def output_results(all_Acc_results):
    
    #Export results to a separate file 
    with open("PHASTER_scan_results.txt","a") as file1:
        for key,values in all_Acc_results.iteritems():
            file1.write("New search: {0}\n".format(key))
            file1.write("Contig\tProphage #\tRegion Length\t# Proteins\tLower Range\tUpper Range\n")
            for result in values:
                file1.write("\t".join(result[:4] + result[4].split("-")) + "\n")
            file1.write("\n\n")                                                       
                                                                                                                                
def main(argv=None):
    
    current_dir = os.getcwd()+'/'
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args, list_input = params.parse_options(argv)
        
        all_Acc_results = {}
        if not list_input:
            Acc_nums = [args[0]]
        else:
            with open(args[0], "rU") as file1:
                lines = file1.readlines()
            Acc_nums = []
            for line in lines:
                if line.strip() != '':
                    Acc_nums.append(line.strip())
        
        for Acc_num in Acc_nums:
            #First link to assembly
            print("Finding assembly for {0}...".format(Acc_num))
            assemblies = link_nucleotide_to_assembly([Acc_num])
            
            #Link back to nucleotide
            print("Finding contig Accs {0}'s assembly...".format(Acc_num))
            IDs = link_assembly_to_nucleotide(assemblies)
            
            print("Searching for prophages with PHASTER for contigs in {0}'s parent assembly...".format(Acc_num))
            for ID_set in IDs:
                #Get Accession number list
                Accs = get_Accs(ID_set)
                results = []
                #Run Phaster searches iteratively
                for Acc_to_search in Accs:
                    PHASTER_file = current_dir+"PHASTER_analysis/" + Acc_to_search.split(".")[0] + ".txt"
                    try:
                        with open(PHASTER_file,'rU') as file1:
                            lines = file1.readlines()
                        skip_entry = False
                    except IOError:
                        lines,skip_entry = query_PHASTER(Acc_to_search,PHASTER_file,current_dir,post=False)
                    if not skip_entry:
                        results = parse_PHASTER(lines,Acc_to_search,results)
                    else:
                        results = []
            all_Acc_results[Acc_num] = results
        
        output_results(all_Acc_results)                                    
                            
                                                        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())
