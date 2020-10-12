import os
from Bio import SeqIO
from CleaveFastaPeptidesFunctions import *


"""
this function take a fasta file of protein sequences and on each seq operate with cleave
to produce the peptides generated from it by a given cleaving rule (protease, efficiency)
then it creates a fasta file of all unique peptides created 
with data regarding how many original sequences in source file contained them and othe calculations
also, it creates a file containing bad seq (probably out of frame and thus containing '*')

input arguments:
path - the source file path
fasta_input - fasta source file
rule : list of dictionaries of proteases rules
        each dict contain:
        rule_name - str - just a name to be used in further files or grafs generated
        rule - str - a string with a regular expression describing the C-terminal site of cleavage.    
        missed_cleavages : int, optional - The maximal number of allowed missed cleavages. Defaults to 0.
        min_length - int, minimum allowed length of peptides generated
query_seq : str, optional (default value is 'CDR3') name of sequence to extract from description for each record in fasta file
clone_identifier: str, optional (default value is 'CDR3', as for right now, clonotyping based on CDR3 element alone) str to extract from description for each record in fasta file

the function first generate the data bases for unique sequences and query peprides and then prints each to a fasta file
the data base is created as a dictinaty, and then written to a fasta file when the key is written as record.seq and value is written as record.description


for peptides file
key - peptide sequence
value - a list of dictionaries, each dictionary contains data regarding a source sequence from which peptide has being produced.
    keys for each dictionaris:
        source_id - id of source in source file
        clone - clone of source sequence
        chain - chain type of source sequence
        isotype - isotype of source sequence
        reads - list of subsets and reads in that format: [(subset,reads),()..] of the source sequence (extracted from source description)
        query_seq - important part of the source sequence based on which further analysis is done (extracted from source description)
        overlaps - a list of dictionaris of all possible overlaps that peptide has with qurey_seq
            kyes:
                is_overlap - 1\0 based on occurance of overlap (at least partial overlap)
                N_tail - sequence in N terminus that does not overlap
                C_tail - sequence in C terminus that does not overlap
                beginning - number - position of peptide beginning in source sequence relative to query sequence
                end - number - position of peptide end in source sequence relative to query sequence
                coverage - percentage of query seq that peptide overlaps
"""
path = 'C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/'
file = 'simulation.fasta'

rules = {'trypsin':r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
         'asp-n':r'\w(?=D)',
         'lys-c':r'K(?!$)',
         'arg-c':r'R(?!$)',
         'glutamyl endopeptidase': r'E(?!$)',
         'chymotrypsin high specificity':r'([FY](?=[^P]))|(W(?=[^MP]))'}
rule = [{'name':'trypsin','rule':rules['trypsin'],'miss_cleavages':0,'min_length':0}]


#p = generate_cleaved_peptide_files(rule[0]['name'],path,file,rule)
def generate_cleaved_peptide_files(rule_name, path, fasta_input, rules_list):
    
    n = 0 #counter for sources number
    total_peptides = 0 #counter for peptides number
    b = 0 #counter for bad sources
    peptides = {}
    input_path = path
    
    #define output folder
    if not os.path.isdir(path + "cleaved_peptides_from_" + fasta_input[:-6] + "/" + rule_name):
        os.makedirs(path + "cleaved_peptides_from_" + fasta_input[:-6] + "/" + rule_name) 
    output_path = path + "cleaved_peptides_from_" + fasta_input[:-6] + "/" + rule_name + "/"
     
    #create fasta file to contain bad sequences (out of frame, containing '*')
        #define output folder
    if not os.path.isdir(path + "cleaved_peptides_from_" + fasta_input[:-6] + "/bad_sequences"):
        os.makedirs(path + "cleaved_peptides_from_" + fasta_input[:-6] + "/bad_sequences") 
    bad_seq_path = path + "cleaved_peptides_from_" + fasta_input[:-6] + "/bad_sequences/"
    
    badseqlog = open(bad_seq_path + "bad_sequnces.fasta", "w")
    
    #define files names (enzymes used and misscleavages allowed for each)
    file_name_str = ''
    for rule in rules_list:
        file_name_str = file_name_str + rule['name'] + '_' + str(rule['miss_cleavages']) + '_' + 'mc_'
    out_pep_file = file_name_str + 'peptides' + '.fasta'

    print('\n')
    print('In silico cleavage for source file: ' + fasta_input)
    print (file_name_str[:-1]) 
    
    
    #iterate over all sequences in source file
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        
#        if n%1000 == 0:
#            print('\n')
#            print(str(n) + ' sequences cleaved')
            
        n+=1
        
        #throw sequences containing '*' to a different file
        if '*' in record.seq:
            badseqlog.write(">" + str(record.id) + " |" + str(record.description)  + "\n" + str(record.seq) + "\n") 
            b+=1
        
        else:
#            print(str(record.id))    
#==============================================================================
#           creates the Peptide_dict and queries dict
#==============================================================================
            pep_list = multi_cleave_peptides(str(record.seq),rules_list)
            total_peptides += len(pep_list)
            pep_set = set(pep_list)
                
            temp_dict = dict((k,record.id) for k in pep_set) #a temp dict of the sequences produced by cleaveing the curret sequence
            for pep in temp_dict.keys(): #iterate over peptides generetade by cleave to check if overlap query sequence
                #creates peptide_dict 
                if pep in peptides.keys():
                    peptides[pep].append({'source_id':temp_dict[pep]}) #if seq already exist from previous cleave, append the seq id to a list (the list is a value in peptides dictionary) 
                else:
                    peptides.update({pep:[{'source_id':temp_dict[pep]}]}) #if seq is new to peptide dictionary - add it
            
    badseqlog.close()
    
    unique_peptides = len(peptides)

    
    print(str(n) + ' record(s) in source file')
    print(str(b) + ' bad record(s) were found and written to a file')
    print(str(total_peptides) + ' total peptides created')
    print(str(unique_peptides) + ' unique peptides were created')
    print('writing peptides to initial database...')
    #creating data-bade fasta files:
    output_db_file(peptides, output_path + out_pep_file)
    print(out_pep_file + ' was written')

    del(peptides)

    return {'outpath':output_path,
            'peptides_file':out_pep_file,
            'total_sequences':n,
            'bad_sequences':b,
            'total_peptides':total_peptides,
            'unique_peptides':unique_peptides,
            }    