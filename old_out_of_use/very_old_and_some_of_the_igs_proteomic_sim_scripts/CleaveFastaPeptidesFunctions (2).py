import re
import itertools as it
import peptide
from collections import deque
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

regex = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
seq = 'AAATGCATCGGACTGTGGAAACCCACTATCAGGCATTCAGGCATCAGCATGCGCCCTATTGCACGC'
cleavage_sites_list = find_unedited_cleavage_sites(regex,seq)



# =============================================================================
# rules = {'trypsin':r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
#          'asp-n':r'\w(?=D)',
#          'lys-c':r'K(?!$)',
#          'arg-c':r'R(?!$)',
#          'glutamyl endopeptidase': r'E(?!$)',
#          'chymotrypsin high specificity':r'([FY](?=[^P]))|(W(?=[^MP]))'}
# =============================================================================


#a = find_unedited_cleavage_sites(regex,seq)
def find_unedited_cleavage_sites(regex,sequnce,cleave_before_pattern=False):
                    
    initial_list = [x.end() for x in regex.finditer(sequnce)]
    initial_set = set(initial_list)
    
    #find other overlapping sites
    if cleave_before_pattern:
        for site in initial_list:
            match1 = regex.match(sequnce,site+1)
            match2 = regex.match(sequnce,site+2)
            if match1:
                initial_set.add(match1.start())
            if match2:
                initial_set.add(match2.start())
                
    else:
        for site in initial_list:
            match1 = regex.match(sequnce,site-1)
            match2 = regex.match(sequnce,site-2)
            if match1:
                initial_set.add(match1.end())
            if match2:
                initial_set.add(match2.end()) 
        
    #return only sites representing in-frame codons
    cleavage_sites = [x for x in initial_set if not x%3 and x!=len(sequence)]
    cleavage_sites.sort()
    return cleavage_sites
    
#peps = cleave_rna_as_peptide(seq,a,0)
def cleave_rna_as_peptide(sequence, cleavage_sites_list, missed_cleavages=0, min_length=None):
    
    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
        
    for i in it.chain(cleavage_sites_list,[None]):
        cleavage_sites.append(i)
#        print(cleavage_sites)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append((seq,cleavage_sites[j],cleavage_sites[j]+len(seq)-1))
    return peptides
    


#peps = cleave_rna_as_peptide(seq,a,0)
def cleave_edited_rna_as_peptide(sequence, cleavage_sites_list, a_to_g_sites, c_to_t_sites, missed_cleavages=0, min_length=None):
    
    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
        
    a_to_g_sites_it = it.chain(a_to_g_sites)
    c_to_t_sites_it = it.chain(a_to_g_sites)
    
    for i in it.chain(cleavage_sites_list,[None]):
        cleavage_sites.append(i)
        print(cleavage_sites)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append((seq,cleavage_sites[j],cleavage_sites[j]+len(seq)-1))
    return peptides


"""

"""
    

"""
return a list of all permutations of a list [1,2] --> [[],[1],[2],[1,2]]
"""
#a_to_g_perm = creat_editing_sites_permutations_one_type(a_to_g_sites)
#c_to_t_perm = creat_editing_sites_permutations_one_type(c_to_t_sites)
def creat_editing_sites_permutations_one_type(lst):  
    editning_perm = []
    for l in range(len(lst)+1):
        [editning_perm.append(list(subset)) for subset in it.combinations(lst, l)]
    return editning_perm
            
            
"""
from two lists 
return a list of tuples, each is a permutation of elements from the two lists.
"""   
#all_perm = creat_editing_sites_permutations_two_type(a_to_g_perm,c_to_t_perm)
def creat_editing_sites_permutations_two_type(lst1,lst2):
    return[(l[0],l[1]) for l in it.product(lst1,lst2)]
    
  
"""
change a sub-sequence given:
the subsequence coordinates relative to the sequence
a_to_g editing coordinates relative to the sequence
c_to_t editing coordinates relative to the sequence
"""
#k = edit_rna_as_peptide('AGTACAGTTTGC', [10,12], [10,13,15], [14,21])
def edit_rna_as_peptide(sequence, sequence_coordinates, a_to_g_list, c_to_t_list):
    splited_ag = [sequence[i:j] for i,j in zip([0] + [x+1-sequence_coordinates[0] for x in a_to_g_list], [x-sequence_coordinates[0] for x in a_to_g_list] + [None])]
    ag_swaped_seq = 'G'.join(splited_ag)
    splited_ct = [ag_swaped_seq[i:j] for i,j in zip([0] + [x+1-sequence_coordinates[0] for x in c_to_t_list], [x-sequence_coordinates[0] for x in c_to_t_list] + [None])]
    new_seq = 'T'.join(splited_ct)
    return new_seq


    
    

"""
CleavageSites
return number of cleavage sites in input peptide (based on input cleavage rule)
"""
def cleavage_sites(pep,rules_list):
    
    cleavage_sites = []
    conc_rule = ''
    
    for rule in rules_list:
        conc_rule = conc_rule + '|' + rule['rule']
    conc_rule = conc_rule[1:]
#    print(conc_rule)
    
    for i in re.finditer(conc_rule, pep):
        cleavage_sites.append(i.start()+1)
    return cleavage_sites


"""
CountReads
count total reads in reads list - [(subset,reads),(),...]
"""
def CountReads(reads_list):
    return (sum(subset_reads[1] for subset_reads in reads_list))



"""
get_clone
get clone ID from description
"""
def get_clone(description, identifier = 'clone: '):
    try:
        return re.findall('(?<=' + identifier + ')[^\s]+',description)[0]
    except:
        return 'unknown'
    

"""
multi_cleave_peptides
operate on a with multiple proteases one after the other
input
peptide - protein sequence
rules_list - list of dictionary each holds the values:
    name - name of protease
    rule - regular expression representing cleavage rule
    miss_cleavages - how many miss cleavages allowed for the specific protease
    overlap_cleave - argument supplied for cleave function - see detailes in function description (usually = False)  
"""

def multi_cleave_peptides(peptide,rules_list):
    
    #theres one rule and therfore operate regularly with cleav
    if len(rules_list) == 1:
        return cleave_rna_as_peptide(peptide, rules_list[0]['rule'], rules_list[0]['miss_cleavages'], rules_list[0]['min_length'])
    
    else:
        pep_list = []
        conc_rules = ''
        min_length = 999999999
        miss_cleavages = 0
            
        #concatenate all rules
        for rule in rules_list:
            conc_rules = conc_rules + '|' + rule['rule']
        conc_rules = conc_rules[1:]
#        print(conc_rules)

        #effective miss_cleavages is the sum of all
        for rule in rules_list:
            miss_cleavages += rule['miss_cleavages']
#        print(str(miss_cleavages))

        #min length is the min among all
        for rule in rules_list:
            if (rule['min_length'] == None) or (rule['min_length'] == 0):
                min_length = None
                break
            if (rule['min_length'] < min_length) and (rule['min_length'] < 999999999):
                    min_length = rule['min_length']

        #create initial pool of peptide
        initial_pep_list = cleave_rna_as_peptide(peptide,conc_rules,miss_cleavages,min_length)
        
        #create a list containing only certain number of miss_cleavages allowd for each rule separately
        for pep in initial_pep_list:
            append_pep_to_list = 1
            for rule in rules_list:
                if append_pep_to_list:
                    if len(cleavage_sites(pep,[rule])) > rule['miss_cleavages']:
#                        print(cleavage_sites(pep,[rule]))
                        append_pep_to_list = 0
#                        print(rule['name'] + ' ' + pep)
                        
            if append_pep_to_list:    
                pep_list.append(pep)
                    
        return pep_list
                    


"""
get_reads
finds in descriptions all sets of numbers in format (X,Y) and return a list of all sets found
"""
def get_reads(description,quant_format = '\((\d+),(\d+)\)'):
    int_reads_list = []
    reads_list = re.findall(quant_format,description)
    for read in reads_list:
        int_read = (int(read[0]),int(read[1]))
        int_reads_list.append(int_read)
    return int_reads_list

 
    
"""
count_reads
for a peptide, after DB was created:
count the number of reads from each subset and create a dictionaty where keys are subsets and values are total read for subset
(peptide could be related to many sources from 1 subset)
"""
def count_reads(sources):
    subsets_reads_dict = {}
    for source in sources: #iterate over source dictionaries
        if source['reads'] == 'reads_not_found':
            pass
        else:
            for subset_data in source['reads']:
                if subset_data[0] in subsets_reads_dict:
                    subsets_reads_dict[subset_data[0]] += subset_data[1]
                else:
                    subsets_reads_dict.update({subset_data[0]:subset_data[1]})
    return subsets_reads_dict
    
    

"""
subsets_reads
from subset_reads_dict created by count_reads
create a string contaning all subsets and total reads of peptide for that subset (sum of all peptide sources read in that subset)
"""
def subsets_reads(subsets_reads_dict):
    Total_reads = 0
    subsets_reads_string = ''
    if len(list(subsets_reads_dict.keys())):
        for key in subsets_reads_dict:
            subsets_reads_string += '(' + str(key) + ',' + str(subsets_reads_dict[key]) + ') '
            Total_reads += subsets_reads_dict[key]
    else:
        subsets_reads_string = 'Unknown reads in subsets'
        
    return subsets_reads_string,Total_reads

  
"""
get_isotype
find isotype in record description 
"""
def get_isotype(description, identifier = 'chain: IG\S '):
    try:
        return re.findall('(?<=' + identifier + ')[A-Z]\d?',description)[0]
    except:
        return 'unknown'

"""
get_query_peptide
find query sequence in record description
if continuation is set to true, the function will return all the sequence feom the last occurrence in whole sequence
of the query seq found in the description, up until the end of the whole seq
"""
def get_query_peptide(description, seq, query_name = 'CDR3', firstdelimiter=': ', continuation=True):
    query = re.findall('(?<=' + query_name + firstdelimiter + ')([A-Z_*]+)(?=\s|$)',description)[0]
    
    if not continuation:
        return query
    else:
        if seq.rfind(query) == -1:
            return query
        else:
            return seq[seq.rfind(query):]

                      
"""
get_chain
find all chains specified in record description
"""
def get_chain(description, name = 'chain', firstdelimiter=': '):  
    try:
        return re.findall('(?<=' + name + firstdelimiter + ')IG\S',description)[0]
    except:
        return 'unknown'
    
    

"""
output_db_file
input - peptides_dict containing sequences as keys and description in some format as values
print a fasta file of all keys as sequences and values as description
"""
def output_db_file(peptides_dict, file_name):
    k=1
    with open(file_name, "w") as handle:
        for pep in peptides_dict:
            rec = SeqRecord(Seq(pep,generic_protein), id = "seq_" + str(k), description = str(peptides_dict[pep]))
            SeqIO.write(rec, handle, "fasta") 
            k += 1    
    handle.close()

    

    
    