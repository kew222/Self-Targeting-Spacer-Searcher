# Self-Targeting Spacer Searcher (STSS)

## About STSS

STSS was written to search prokaryotic genomes for CRISPR arrays and determine if any of the spacers in the array target the organism's own genome. Each self-targeting spacer (STS) that is found is checked for mutations in the repeat, target sequence, or PAM as well as for missing Cas genes. What gene(s) is/are being targeted is also determined and whether the targeted positions occur in a prophage are checked with PHASTER (Arndt, D., Grant, J., Marcu, A., Sajed, T., Pon, A., Liang, Y., Wishart, D.S. (2016) PHASTER: a better, faster version of the PHAST phage search tool. Nucleic Acids Res., 2016 May 3.).

### Installation

#### Requirements:

- Python 2 (Python 3 may work but has not been thoroughly tested)
- Biopython
- requests (Python package)
- blastn (available from NCBI)
- Clustal Omega
- HMMER 3
- Internet connection

#### Instructions

First, clone the STSS repository with:

`git clone https:github.com/kew222/Self-Targeting-Spacer-Search-tool.git`

Install blastn from NCBI (ftp:ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

For recognition by STSS, all of the required binaries need to be visible in the bin/ directory in the repository. We recommend installing or placing all binaries in /usr/local/bin/ and creating a link to the STSS bin/ directory with (as an example):

`ln -s /usr/local/bin/blastn STSS/bin/blastn`

If not storing the binaries at /usr/local/bin is undesired, simply replace with the location of the binaries. Alternatively, the binaries can be placed directly in the bin/ directory in STSS.

Install Clustal Omega from (http:www.clustal.org/omega/#Download). Again, precompiled binaries are available for quick use, although compiling from source is also an option. After choosing the appropriate operating system, link to bin/ with (ex.):

`ln -s /usr/local/bin/clustalo STSS/bin/clustalo`

The last required non-Python binary is HMMER3 (http:hmmer.org/download.html) for which binaries are again available. Following the same method as before, move the binaries (or compile) to a centralized location and link at minimum nhmmscan and hmmscan to bin/, for example:

`ln -s /usr/local/bin/*hmm* STSS/bin/`

With all of the binaries in place for STSS, all that remains to set up the Python packages to run. This can be done by simplying running:

`python setup.py install` 

from the STSS directory. This will check for Biopython and the requests package and install them if they are missing. At this point everything should be installed.

Before running STSS, however, you will need to edit user_email.py (using any standard text editor) to input your email address, which is needed for running the NCBI tools. This will only need to be done once. STSS has no email collecting code, etc. so we will never see your email address. 


### Workflow

In order to find STSs, STSS goes through the following steps (described in more detail in: **PUBLICATION HERE**):

1. Gathers genomes (either provided or by searching NCBI)
2. Uses the CRISPR Recogition Tool (CRT; Bland, C. et al. CRISPR Recognition Tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats. BMC Bioinformatics 8, 209–8 (2007)) to identify CRISPR arrays
   - Because CRT cannot handle genomes with long degenerate stretches, these are effectively masked before CRT analysis
   - Also, CRT has a propensity to put repeats in the spacer sequences. Arrays are double-checked for similar sequences at the beginning and ends of the spacers to find these mistakenly incorporated repeats and fixes them.
3. Performs a BLAST search on the genome with each spacer found (must be outside bounds of found arrays)
4. The array containing the STS is checked for Cas gene content using HMMER3 or the Conserved Domain Database (CDD; NCBI)
5. The array consensus repeat is determined and the repeats on either side of the STS are checked
6. The direction of the array is determined with either the Cas genes and/or the repeat sequence
7. The targeted gene in the genome is determined, possibly using the CDD to find domains if not annotated
8. Any mutations between the target and the guide RNA are determined as well as the upstream/downstream neighboring sequences
9. Last, the targeted contig is checked for prophages with PHASTER to determine if the STS is in a prophage


### Usage

#### Running STSS

STSS can accept a couple different formats for searching genomes for self targeting spacers. In the simplist case, the user might want to search for STSs in all of one organism say _E. coli_. To search all _E. coli_ genomes:

`python STSS.py --search "Escherichia coli"`

Note that the quotes will be required for multi-word terms. Also, the search term is used to search the [organism] tag on NCBI.

Alternatively, the user could determine the genomes he/she wants to search and download them (fasta format!). In that case, STSS can be given a directory containing genomes to search in:

`python STSS.py --dir downloaded_genomes/`

This method is also useful is the genomes to be searched are not on NCBI. However, some of the capabilities of STSS will be limited if there isn't an NCBI Accession number to use or it can't determine it from the fasta file. This can be somewhat overcome though if there is a GenBank formatted file (.gb) provided with the same name for STSS to find.

Last, STSS can also accept a list of assemblies (NCBI Accessions) to get the genomes of interest (assuming they are listed in a file named assemblies.txt):

`python STSS.py --Accs assemblies.txt`


The results from any running method will be same. There will be two tab-delimited files with the same formats, one containing the STSs that were found to be in prophages and the other that were not. If the PHASTER analysis was skipped (see options below), only one file will be output. 

#### Options

Option |  Description
-------| ------------
-h, --help               |       Opens help message  
-v, --version              |     Displays version number  
--dir <directory>            |   Use directory of genomes (fasta-formatted) instead of searching NCBI  
--search <"NCBI search term"> |  Use NCBI nucleotide database to find genomes  
--Accs <Assembly_list_file> |    Search genomes based on a given list of assemblies (incompatible with search)  
-o, --prefix <string>        |   Prefix for filenames in output (ex. prefix of 'bagel' gives: bagel_Spacers...islands.txt)  
-f, --force-redownload      |    Forces redownloading of genomes that were already downloaded  

Use -f if using search and want to force the redownload of any genomes that already exist in the downloaded_genomes directory 

-n, --no-ask            |        Force downloading of genomes regardless of number found (default: ask) 

By default, STSS will ask the user if he/she wants to continue with the download if there are large number of files returned. There is a delay while searching NCBI, so turning the option off will prevent the need to wait to confirm the download if hard drive space is not an issue. 

-l, --limit <N>         |        Limit Entrez search to the first N results found (default: 10000) 

Lowering this value is unnessecary for small searches, but may need to be raised for large scale searches.  

--CDD                    |       Use the Conserved Domain Database to identify Cas proteins (default is to use HMMs) 

By default, HMMs are used to try to identify Cas proteins near arrays or in the genome (see -d). However, the CDD can also be used. The main difference other than the use of HMMs vs. PSWMs, is which database was updated most recently CDD can be slow to update, but depending on how often the provided HMMs are updated, it may be more useful to use CDD. Using HMMER is much faster than the CDD due to the need to use webservers for the CDD. 

--complete-only         |        Only return complete genomes from NCBI 
--rerun-loci <filename>    |     Rerun the locus annotater to recheck the nearby locus for Type and completeness from provided spacer search results file 

When --rerun-loci is used, a results file from a previous run must be given. This option will rerun all of the data collection steps (steps 4-9 in Workflow above) using the previously found STSs.  

-E, --E-value  <N>        |      Upper limit of E-value to accept BLAST results as protein family (default: 1e-4) 

Note that this E-value is used for most of the scans (HMM, CDD, BLAST, etc.), so changing it for one will change it for all. 

--percent-reject <N>       |     Percentage (of 100%) to use a cutoff from the average spacer length to reject validity of a proposed locus (default 25%). Lower values are more strigent. 

Adjust the percent-reject value will tune how much deviation is allowed from the average spacer length to try to determine array false positives. A lower percent will be more stringent and reject more arrays as false positives. This was an ad-hoc addition after we observed that most of the arrays identified by CRT that was not actually CRISPR arrays typically had a spacers with a variety of lengths, while correct arrays don't vary as much. Note, however, that Class 1 arrays do have some natural variability in length, so being too stringent can reject good arrays. 

-s, --spacers <N>       |        Number of spacers needed to declare a CRISPR locus (default: 3)  
--pad-locus <N>         |        Include a buffer around a potential locus to prevent missed repeats from appearing as hits (default: 100) 

This pad-locus value is an ad-hoc correction to prevent spacers near the edges of arrays from being picked up as false positives. The main reason for including this factor is due to how CRT determines position, which ignores separations between contigs. 

-d, --Cas-gene-distance <N>  |   Window around an array to search for Cas proteins to determine CRISPR subtype (default: 20000 - input 0 to search whole genome)  

When determining the array Type, the genes up- and downstream of identified arrays are checked to see if and what Cas genes they contain using HMMs or PSWMs. By default, STSS searches 20k bases away from the start of an array. However, there are cases where a lot of Cas genes are present for an array that could cause some genes to be missed outside this range. In these cases, it may be avisable to increase the search distance. As a second option, all of the contigs of a genome (i.e., the whole genome) can be searched if '0' is input. Be aware, however, that STSS will still try to guess a CRISPR type based on the Cas genes identified, and may be meaningless if multiple CRISPR systems exist in the genome. 

--skip-PHASTER                |  Skip PHASTER analysis (currently can't upload search files)  
-p, --rerun-PHASTER <filename> | Rerun PHASTER to recheck islands from provided Spacer search results file  

Only reruns the PHASTER analysis and does not recheck anything else in the results. Requires a results file from a previous run in the STSS output format.









