#!/usr/bin/env python

#Compiles all unique results from all spacers for import into Excel

#This script was written in 2016 by Kyle Watters
#Copyright (c) 2016 Kyle Watters. All rights reserved.

import glob
import sys
import getopt

help_message = '''
Spacer_data_compiler.py creates a compiled list of spacers, removing duplicates. It can accept files as input to
exclude previously recorded spacers (to prevent copy-paste duplicates)

Usage:
    Spacer_data_compiler.py [ <file1> <file2> ... <filen> ]
    
'''

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    
    def __init__(self):
        pass
    
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "h",
                                       ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)  
                                
        return args
        
    def check(self):
        pass

def main(argv=None):

    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()
        
        exclude_list = []
        for exclude_file in args:
            with open(exclude_file, 'r') as exclude:
                lines = exclude.readlines()
            for line in lines:
                if line[:16] != "Target Accession#" and line not in exclude_list:  #Exclude the header line
                    exclude_list.append(line)
        
        print("Please wait...compiling self-targets list")
        files = glob.glob("*/Spacers_*.*")   #gets all Spacer results files in the current directory 
        all_results = []   
        for eachfile in files:
            with open(eachfile,'r') as readingfile:
                lines = readingfile.readlines()
            for line in lines:
                if line not in all_results+exclude_list:   #Only take new results
                    all_results.append(line)
        with open("compiled_Spacers.txt", 'w') as writingfile:
            for line in all_results:
                writingfile.write(line)
                
        #Also count the number of unique genomes scanned
        genomes = []
        print("Please wait...compiling genomes list")
        files = glob.glob("*/genomes_analyzed.txt")
        for eachfile in files:
            with open(eachfile, 'r') as readingfile:
                lines = readingfile.readlines()
            for genome in lines:
                if genome != '\n' and genome.strip() not in genomes:
                    genomes.append(genome.strip())
        print("{0} unique genomes have been searched.".format(len(genomes)))
        with open("Unique_genomes_scanned.txt", 'w') as writingfile:
            for genome in genomes:
                writingfile.write(genome + '\n')
    
    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
        return 2

if __name__ == "__main__":
    sys.exit(main())

