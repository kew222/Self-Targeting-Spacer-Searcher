#functions that were originally written to look for proteins that re-appear in self-targeting genomes


def mine_proteins(protein_list,Acc_to_search,proteins_found,region,hit_num=None,all_islands=False):
    
    #Now want to mine all of the unassigned/predicted proteins 
    webcall = "http://phaster.ca/jobs/{0}/detail.txt".format(Acc_to_search)
    tries = 0
    while True:
        tries += 1
        try:
            q = requests.get(webcall, timeout=20)  #do a PHASTER search in the potential hits genome
            break
        except (requests.exceptions.ReadTimeout, requests.exceptions.ConnectTimeout, requests.exceptions.ConnectionError):
            time.sleep(5)  #wait 5 seconds and retry
        if tries > 3:
            print("PHASTER server not responding to query for details on {0}. Skipping...".format(Acc_to_search))
            break
    q2 = q.text
    q3 = q2.split("\n")
    protein_list.append([proteins_found])
    right_region = False
    for protein in q3:
        if right_region == True and protein.strip() != '' and protein.strip()[:3] != '###':
            island_protein = [str(x).strip() for x in filter(None, protein.split("  "))]
            phrases_to_search = ["hypothetical", "phage-like"]
            for phrase in phrases_to_search: 
                if island_protein[1].find(phrase) > -1:
                    protein_list[hit_num].append(island_protein) 
        elif right_region == False:
            if protein[:12] == "#### region ":
                if all_islands:
                    right_region = True   #if recording all proteins from islands, skip finding if in same one as spacer
                else:                                        
                    current_region = int(protein.split("region ")[1].split(" ####")[0])
                    if current_region == region: ##That is, is the current region cycling through the one the spacer is in
                        right_region = True
        elif protein.strip() == '' or protein.strip()[:3] == "###" and right_region == True:
            right_region = False
            if not all_islands:
                break     #found the right region and all the proteins have been checked
    hit_num+=1
    return protein_list,hit_num  
    
def family_cluster(protein_list,E_value_limit=1e-3):

    #Iteratively BLAST all proteins against all other proteins, but only blast those that haven't been aligned yet
    #If makes a certain cutoff, store as a site for that protein
    protein_hits_dict = {}
    protein_hits_list = []
    for spacer in protein_list:
        store_spacer = spacer[0]
        for index in range(1,len(spacer)):
            protein_hits_dict[spacer[index][1].replace(" ","_")] = [spacer[index][3], store_spacer]  #makes a dictionary of proteins with name, sequence in each element
            protein_hits_list.append([spacer[index][1],spacer[index][3]])  #name, protein sequence
    BLAST_file = "subject_list.fasta"  
    query_file = "query.fasta"
    
    #Now, split the fasta file entry by entry, so only aligning down the list
    families = []     #blast check for identities
    protein_counter = 0
    for protein_num in protein_hits_list:
        #write a short file for the query
        with open(query_file, "w") as compiled_file:  #need to write a file for the blastp input (can't pass for multiple)
            name = protein_hits_list[protein_counter][0].replace(" ","_")
            AAseq = protein_hits_list[protein_counter][1]
            compiled_file.write(">{0}\n{1}\n".format(name,AAseq))   
        
        #write out the rest of the proteins in a separate subject file
        with open(BLAST_file, "w") as compiled_file2:  #need to write a file for the blastp input (can't pass for multiple)
            for protein in protein_hits_list[protein_counter+1:]:
                name = protein[0].replace(" ","_")
                AAseq = protein[1]
                compiled_file2.write(">{0}\n{1}\n".format(name,AAseq))     
                
        blast_cmd = "blastp -query {0} -subject {1} -outfmt 6".format(query_file,BLAST_file)
        handle = subprocess.Popen(blast_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        output, error = handle.communicate()
        
        temp_list = [protein_num[0].replace(" ","_")]
        for line in output.split("\n"):
            if line.strip() != '':
                result = line.split('\t')
                E_value = float(result[-2])
                target = result[1]
                #query, target, ident, length, bit = result[0], result[1], result[2], result[3], float(result[-1]) 
                if E_value <= E_value_limit:  #default is 1e-3
                    temp_list.append(target)
        fam_num = 0
        no_family_match = True
        for family in families:
            ## look at protein 1 vs. 2-99, save list of matches E-value <= 1e-3 as 1, put blank if none
            ## look at protein 2 vs. 3-99, obtain  list of matches E-value <= 1e-3 
            #   check group 1 if contains 2, if yes, add to previous group (1), if in no groups make new group (2)
            for member in family:
                if protein_num[0] == member:
                    families[fam_num].append(temp_list)
                    no_family_match = False   
            fam_num += 1
        if no_family_match == True and len(temp_list) > 1:
            families.append(temp_list) 
        protein_counter += 1

    #Export in-island blast results to fasta format
    with open("full_protein_list.txt","w") as ifile:
        for line in protein_list:
            y = 1
            for x in line:
                if y > 1:
                    ifile.write(">"+str(x[1]).replace(" ","_")+"\n" + str(x[3]) + "\n")
                y += 1    
                
    #check for duplicates in the alignment matrix (by sequence) and remove them, remove family if only one entry remains after duplicates
    family_no = 0
    for family in families:
        num_members = len(family)
        member_no = 0
        for member in family:
            member_seq = protein_hits_dict[member][0]
            for x in range(num_members-1,member_no,-1):  #search backwards so the index numbers don't change as being removed.
                check_seq = protein_hits_dict[family[x]][0]
                if member_seq == check_seq:
                    del families[family_no][x]
                    num_members -= 1  #if removed, one less
            member_no += 1 #keeps track of position of the member being checked
        family_no += 1 
    ordered_families = sorted(families, key=len, reverse=True)
    families = []
    for family in ordered_families:
        if len(family) > 1:          #remove single element families
            families.append(family)          
    print("Finished assigning protein families.")
    return families,protein_hits_dict,query_file,BLAST_file

def families_print(families,protein_hits_dict,families_limit):

    #Output the lists of protein families with details
    families_file = "potential_families.txt"
    fam_num = 1
    with open(families_file, "w") as fileobj:    
        fileobj.write("GI #\tAccession #\tCRISPR #\tSpacer #\tSpacer Locus Pos.\tSpacer Genome Pos.\tSpacer Sequence\tPAM Region\tPHASTER Island\tProtein Name\tAA Sequence")
        for family in families:
            if fam_num <= families_limit:
                fileobj.write("\nCandidate Family {0}".format(fam_num))
                for member in family:
                    member_seq = protein_hits_dict[member][0]
                    details = protein_hits_dict[member][1]
                    member_details = [str(x) for x in details]
                    fileobj.write("\n{0}\t{1}\t{2}".format("\t".join(member_details),member,member_seq))
                fileobj.write("\n")  #add a space to separate families
            fam_num += 1

    #Also create a file that has just the first member of each family in fasta format to search for conserved domains, etc.
    with open("family_representatives.fasta", "w") as domain_file:
        family_no = 1
        for family in families:
            name = "Family {0} Representative".format(family_no)
            seq = protein_hits_dict[family[0]][0]
            domain_file.write(">{0}\n{1}\n".format(name, seq))
            family_no += 1

def families_search(families,families_limit,current_dir,search='',E_value_limit=1e-3):
    
    if families == []:
            print("No families found.")             ####This part should be recoded with CDD
            return
    if not os.path.exists("Family_BLASTs"):
        os.mkdir("Family_BLASTs")
    print("Waiting for BLAST of families against NCBI database....")
    with open("family_representatives.fasta", 'rU') as inputfile:
        fasta2 = SeqIO.parse(inputfile, 'fasta')
        family_no = 1
        for record in fasta2:
            if family_no <= families_limit:
                if search != '':
                    ignore_results = "NOT {0}".format(search)
                else:
                    ignore_results = ''
                blastp2 = NCBIWWW.qblast("blastp", 'nr', record, entrez_query=ignore_results, expect=E_value_limit, hitlist_size=10, format_type="Text")
                with open(current_dir+"Family_BLASTs/family_{0}_BLAST.txt".format(family_no), "w") as save_file: 
                    for lines in blastp2:
                        save_file.write(lines)    
                family_no += 1
    blastp2.close()
    print("Completed BLAST of family representatives to NCBI database")


def families_alignment(families,protein_hits_dict,current_dir):

    #Now use either BLAST or Clustal Omega to do the alignment and report results
    fam_num = 1
    results_dir = current_dir+"Candidate_family_alignments/"
    if not os.path.exists(results_dir):
        os.mkdir("Candidate_family_alignments")
    for family in families:
        if len(family) > 2:
            #Use Clustal Omega to compare 3+
            output_file = results_dir + "family_{0}.aln".format(fam_num)
            with open(current_dir+"clustal_temp.fasta","w") as holder:
                for member in family:
                    holder.write(">{0}\n{1}\n".format(member,protein_hits_dict[member][0]))        
            clustal_cmd = "clustalo -i {0} -o {1} --force".format("clustal_temp.fasta",output_file)
            clo = subprocess.Popen(clustal_cmd.split())
            clo.communicate()
        else:
            #Use BLAST to compare two sequences
            output_file = results_dir + "family_{0}.txt".format(fam_num)
            with open("Blast_q_temp.fasta","w") as holder:        
                holder.write(">{0}\n{1}\n".format(family[0],protein_hits_dict[family[0]][0]))
            with open("Blast_s_temp.fasta","w") as holder:  
                holder.write(">{0}\n{1}\n".format(family[1],protein_hits_dict[family[1]][0]))
            blast_cmd = "blastp -query {0} -subject {1} -outfmt 0".format("Blast_q_temp.fasta","Blast_s_temp.fasta",output_file)
            handle = subprocess.Popen(blast_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
            output, error = handle.communicate()
            with open(output_file,"w") as result:
                result.write(output.decode())
        fam_num += 1
    print("Family alignments complete.")

def anti_CRISPR_cluster_tool(protein_list,E_value_limit,families_limit,search='',skip_family_search=True,skip_family_create=True,skip_alignment=True):

    if protein_list == []:
        print("No candidate proteins found. Exiting...\n")
        return
    
    #BLAST all the protein results against themselves to see if any proteins group together
    families,protein_hits_dict,query_file,BLAST_file = family_cluster(protein_list,E_value_limit)
 
    #Output information about the families that were determine
    families_print(families,protein_hits_dict,families_limit)
 
    #Compare families to what's on NCBI (this sould be CDD, not blast)
    if not skip_family_search:
        families_search(families)
    
    #Align the families using ClustelO or BLAST
    if not skip_alignment:
        families_alignment(families,protein_hits_dict)
    
    #Clean up leftover files
    for xx in [query_file,BLAST_file,"clustal_temp.fasta","Blast_q_temp.fasta","Blast_s_temp.fasta"]:
        if os.path.isfile(xx):
            os.remove(xx)  #get rid of the temp files
