#This file contains the CRISPR Type definitions by protein
#The format below is preserved, but can be edited appropriately to allow changing of the definitions
#Current implementation is based on:
#        Makarova, et al. 2015 (Nat. Rev. Micro.)
#        Makarova, et al. 2011 (Nat. Rev. Micro.)   - Not all of these names are used
#        Makarova, et al. 2015 (Biology Direct)
#        Shmakov, et al. 2015 (Mol. Cell)


#Cas11 could include:  Csa5 (I-A), Cse2 (I-E), Csm2 (III-A) and Cmr5 (III-B) 
#Cas5 includes: cas5a, cas5d, cas5e, cas5h, cas5p, cas5, cmx5, CasD
#Cas7 includes: csa2, csd2, cse4, csh2, casC, csp1, cst2
#Cas8a1 includes: cmx1, cst1, csx8, csx13, NCXXC-CXXC
#Cas8a2 includes: Csa4, Csx9
#Cas9 could include: Csx12


CRISPR_types = {

                "Type I-A":[["Cas6"],["Csa5"],["Cas7"],["Cas5"],["Cas8a1","Cas8a2","Csa4","Csx9"],["Cas3'"],["Cas3''"],["Cas2"],["Cas4","Csa1"],["Cas1"],["Cas4"]],
                "Type I-B":[["Cas6"],["Cas8b1","Csh1"],["Cas7","Cst2"],["Cas5"],["Cas3"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type I-C":[["Cas3"],["Cas5"],["Cas8c","Csd1","Csp2"],["Cas7"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type I-U":[["Cas3"],["Cas8c"],["Cas7"],["Cas5","GSU0054"],["Cas6"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type I-D":[["Cas3'"],["Cas3''"],["Cas10d","Csc3"],["Cas7","Csc2"],["Cas5","Csc1"],["Cas6"],["Cas4"],["Cas1"],["Cas2"]],
                "Type I-E":[["Cas3"],["Cas8e","Cse1"],["Cse2","CasB"],["Cas7","Cse4","CasC"],["Cas5","CasD"],["Cas6","Cse3","CasE"],["Cas1"],["Cas2"]],
                "Type I-F":[["Cas1"],["Cas2","Cas3"],["Cas8f","Csy1"],["Cas5","Csy2"],["Cas7","Csy3"],["Cas6f","Csy4"]],
                "Type I-G":[["Cas6"],["Cas8a1","Cst1"],["Cas7","Cst2"],["Cas5","Cas5t"],["Cas3"],["Cas4"],["Cas1"],["Cas2"]],
                
                "Type II-A":[["Cas9","Csn1"],["Cas1"],["Cas2"],["Csn2"]],
                "Type II-B":[["Cas9","Csn1"],["Cas1"],["Cas2"],["Cas4","Csa1"],],
                "Type II-C":[["Cas9","Csn1"],["Cas1"],["Cas2"]],
                
                "Type III-A":[["Cas6"],["Cas10","Csm1"],["Csm2"],["Cas7","Csm3"],["Cas5","Csm4"],["Cas7","Csm5"],["Csm6"],["Cas1"],["Cas2"]],
                "Type III-B":[["Cas7","Cmr1"],["Cas10"],["Cas5","Cmr4"],["Cmr5"],["Cas6"],["Cas7","Cmr6"],["Cas1"],["Cas2"]],
                "Type III-C":[["Cas7","Cmr1"],["Cas7","Cmr6"],["Cas10"],["Cas7","Cmr4"],["Cmr5"],["Cas5","Cmr3"]],
                "Type III-D":[["Cas10"],["Cas7","Csm3"],["Cas5","Csx10"],["Csm2"],["Cas7","Csm3"],["Cas7","Csm3"],["all1473"],["Cas7","Csm3"]],
                
                "Type IV-A":[["dinG","Csf4"],["Csf1"],["Cas7","Csf2"],["Cas5","Csf3"]],
                
                "Type V-A":[["Cas12a","Cpf1"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type V-B":[["Cas12b","C2c1"],["Cas4","Csa1"],["Cas1"],["Cas2"]],
                "Type V-C":[["Cas1"],["Cas12c","C2c3"]],
                
                "Type VI-A":[["Cas13a1","C2c2"],["Cas1"],["Cas2"]],
                "Type VI-B":[["Cas13b","C2c2"],["Csx27"],["Csx28"],["Cas1"],["Cas2"]],
                "Type VI-C":[["Cas13a2","C2c2"],["Cas1"],["Cas2"]]
                }
                
                
Cas_proteins =  {
                
                "Csa1":["Type I-A","Type I-B","Type I-C","Type I-D","Type II-B","Type V-A","Type V-B"],
                "Csa4":["Type I-A"],
                "Csa5":["Type I-A"],
                "Csc1":["Type I-D"],
                "Csc2":["Type I-D"],
                "Csc3":["Type I-D"],
                "Csd1":["Type I-C"],
                "Cse1":["Type I-E"],
                "Cse2":["Type I-E"],
                "Cse3":["Type I-E"],
                "Cse4":["Type I-E"],
                "Csh1":["Type I-B"],
                "Csf1":["Type IV-A"],
                "Csf2":["Type IV-A"],
                "Csf3":["Type IV-A"],
                "Csf4":["Type IV-A"],
                "Csn2":["Type II-A"],
                "Csp2":["Type I-C"],
                "Csy1":["Type I-F"],
                "Csy2":["Type I-F"],
                "Csy3":["Type I-F"],
                "Csy4":["Type I-F"],
                "Csn1":["Type II-A","Type II-B","Type II-C"],
                "Csm1":["Type III-A"],
                "Csm2":["Type III-A","Type III-D"],
                "Csm3":["Type III-A","Type III-D"],
                "Csm4":["Type III-A"],
                "Csm5":["Type III-A"],
                "Csm6":["Type III-A"],
                "Cst1":["Type I-G"],
                "Cst2":["Type I-G","Type I-B"],
                "Csx9":["Type I-A"],
                "Csx10":["Type III-D"],
                "Cmr1":["Type III-B","Type III-C"],
                "Cmr3":["Type III-C"],
                "Cmr4":["Type III-B","Type III-C"],
                "Cmr5":["Type III-B","Type III-C"],
                "Cmr6":["Type III-B","Type III-C"],
                "GSU0054":["Type I-U"],
                "all1473":["Type III-D"],
                "dinG":["Type IV-A"],
                "Cpf1":["Type V-A"],
                "C2c1":["Type V-B"],
                "C2c3":["Type V-B"],
                "C2c2":["Type VI-A", "Type VI-B","Type VI-C"],
                "Cas13a1":["Type VI-A"],
                "Cas13b":["Type VI-B"],
                "Cas13a2":["Type VI-C"],
                "Cas1":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type II-A","Type II-B","Type II-C","Type III-A","Type III-B","Type V-A","Type V-B","Type V-C","Type VI-A"],
                "Cas2":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type II-A","Type II-B","Type II-C","Type III-A","Type III-B","Type V-A","Type V-B","Type VI-A"],
                "Cas3":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G"],                    
                "Cas3'":["Type I-A","Type I-D"],
                "Cas3''":["Type I-A","Type I-D"],
                "Cas4":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-G","Type II-B","Type V-A","Type V-B"],               
                "Cas5":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D","Type IV-A"],                              
                "Cas5t":["Type I-G"],
                "Cas6":["Type I-A","Type I-B","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B"],                                             
                "Cas6f":["Type I-F"],
                "Cas7":["Type I-A","Type I-B","Type I-C","Type I-U","Type I-D","Type I-E","Type I-F","Type I-G","Type III-A","Type III-B","Type III-C","Type III-D","Type IV-A"],
                "Cas8a1":["Type I-A","Type I-G"],   
                "Cas8a2":["Type I-A"],   
                "Cas8b1":["Type I-B"],                                                                                
                "Cas8c": ["Type I-C","Type I-U"],                                                                                                                                           
                "Cas8e": ["Type I-E"],                                                                                                                                                                                                       
                "Cas8f": ["Type I-F"],                                                                                                                                                                                                                                                                   
                "Cas9":["Type II-A","Type II-B","Type II-C"],                                                                                                                                                                                                                                                                                                                
                "Cas10":["Type III-A","Type III-B","Type III-C","Type III-D"],
                "Cas10d":["Type I-D"],
                "CasB":["Type I-E"],
                "CasD":["Type I-E"],
                "CasC":["Type I-E"],
                "CasE":["Type I-E"]
                                                                                                                                                                                                                                                                                                           
              }               
                
Cas_synonym_list = {

                "Csa1":"Cas4",
                "Csa4":"Cas8a",
                "Csx9":"Cas8a",
                "Csh1":"Cas8b1",
                "Csd1":"Cas8c",
                "Csp2":"Cas8c",
                "GSU0054":"Cas5",
                "Csc3":"Cas10d",
                "Csc2":"Cas7",
                "Csc1":"Cas5",
                "Cse1":"Cas8e",
                "Cse3":"Cas6",
                "Cse4":"Cas7",
                "CasB":"Cse2",
                "CasC":"Cas7",
                "CasD":"Cas5",
                "CasE":"Cas6",
                "Csy1":"Cas8f",
                "Csy2":"Cas5",
                "Csy3":"Cas7",
                "Csy4":"Cas6f",
                "Cst1":"Cas8a1",
                "Cst2":"Cas7",
                "Cst5t":"Cas5",
                "Csn1":"Cas9",
                #skipping the Csm/Cmr proteins, as it is helpful to see them as is
                "Csf4":"dinG",
                "Csf2":"Cas7",
                "Csf3":"Cas5",
                "C2c2":"Cas13"
                
                }                

#According to  Lange, et al. NAR, 2013: "CRISPRmap: an automated classification of repeat con- servation in prokaryotic adaptive immune systems" 
Repeat_families_to_types = {                
                                "F1": ["I-B", "III-A", "III-B"],
                                "F2": ["I-E"],
                                "F3": ["I-C"], 
                                "F4": ["I-C", "I-E", "II-B"],
                                "F5": ["I-F"], 
                                "F6": ["I-A"], 
                                "F7": ["I-A"], 
                                "F8": ["I-F"], 
                                "F9": ["III-B"], 
                                "F10": ["I-B", "III-B"], 
                                "F11": ["III-B"], 
                                "F12": ["II-B", "III-A"], 
                                "F13": ["I-A", "III-B"], 
                                "F14": ["I-A", "I-D", "III-A"],
                                "F15": ["I-A", "III-B"], 
                                "F16": ["III-A"], 
                                "F17": ["?"], 
                                "F18": ["I-E", "II-B"], 
                                "F19": ["?"], 
                                "F20": ["I-B"], 
                                "F21": ["I-E"], 
                                "F22": ["I-E"], 
                                "F23": ["I-D", "II-B"], 
                                "F24": ["III-A", "III-B"], 
                                "F25": ["I-A", "II-B", "III-A"],
                                "F26": ["?"], 
                                "F27": ["II-A", "II-B"], 
                                "F28": ["I-A"], 
                                "F29": ["III-A"], 
                                "F30": ["?"], 
                                "F31": ["III-A"], 
                                "F32": ["I-C"], 
                                "F33": ["I-C", "I-E", "II-B"],
                                "F34": ["II-B"], 
                                "F35": ["II-A"], 
                                "F36": ["?"], 
                                "F37": ["I-C", "III-B"], 
                                "F38": ["I-A", "III-B"], 
                                "F39": ["I-A", "I-B", "II-B"],
                                "F40": ["I-B"]        
                            }         
                             
#Here, 1 means the Cas genes are expected to be upstream, while -1 is downstream (to be consistent with up_down), 0 means it is common to see both orientations                                                                                
Expected_array_directions = {
                                "Type I-A":1,
                                "Type I-B":1,
                                "Type I-C":1,
                                "Type I-U":1,
                                "Type I-D":1,
                                "Type I-E":1,
                                "Type I-F":1,
                                "Type I-G":1,
                                
                                "Type II-A":1,
                                "Type II-B":-1,
                                "Type II-C":-1,
                                
                                "Type III-A":1,
                                "Type III-B":1,
                                "Type III-C":1,
                                "Type III-D":1,
                                
                                "Type IV-A":1,
                                
                                "Type V-A":1,
                                "Type V-B":1,
                                "Type V-C":1,
                                
                                "Type VI-A":1,
                                "Type VI-B":1,
                                "Type VI-C":1
                                
                            }
                            
                                                        
                        