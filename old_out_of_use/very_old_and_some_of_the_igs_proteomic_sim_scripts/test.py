import re
import itertools as it
import sys
import peptide
from collections import deque
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_rna


sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/' +'test.txt', 'w')

regex = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
sequence = 'AAATGCATCGGACTGTGGAAACCCACTATCAGGCATTCAGGCATCAGCATGCGCCCTATTGCACGC'
cleavage_sites_list = find_cleavage_sites(regex,sequence)
missed_cleavages = 1
min_length = None
a_to_g = [0,1,6,18,19,20,30,38,45,57,62]
c_to_t = [3,23,25,41,51,53]
look_ahead_for_editnig = 3
look_behing_for_editing = 3

#here it begins


peptides = []
edited_peptides = []    
ml = missed_cleavages+2 #size of normal window (2 boundaries + missed cleavage sites)
trange = range(ml) 
cleavage_sites = deque([0], maxlen=ml) #the normal, dynamic window for digestion
cl = 1
k = 0        
l = 0
all_a_to_g_sites_it = it.chain(a_to_g)
all_c_to_t_sites_it = it.chain(c_to_t)
a_to_g_sites = deque([])
c_to_t_sites = deque([])
a_to_g_current_site = None
c_to_t_current_site = None

n=0
#run over all cleavage site and append to deque with size = ml (2 boundaries + missed cleavage sites)
for i in it.chain(cleavage_sites_list,[None]):

    print('========================================================================================================')
    
    #current window in sequence
    cleavage_sites.append(i)
    print('\ncleavage sites window is: ' + str(cleavage_sites))
    
    #window bounderies
    if cleavage_sites[-1] == None:
        downstream_bound = len(sequence)-1
    else:
        downstream_bound = cleavage_sites[-1] + look_ahead_for_editnig - 1
    if cleavage_sites[0] - look_behing_for_editing < 0:
        up_stream_bound = 0
    else:
        up_stream_bound = cleavage_sites[0] - look_behing_for_editing
        
    window_range = (up_stream_bound,downstream_bound)
    print('window coor: ' + str(window_range))
    
    #editing sites in current window
    a_to_g_current_site, a_to_g_sites = update_editing_sites_in_window(window_range,a_to_g_sites,all_a_to_g_sites_it,a_to_g_current_site)
    c_to_t_current_site, c_to_t_sites = update_editing_sites_in_window(window_range,c_to_t_sites,all_c_to_t_sites_it,c_to_t_current_site)    
    print('a to g sites are ' + str(a_to_g_sites))
    print('c to t sites are ' + str(c_to_t_sites))
    
    #find all eiditing permutations, a permutation is ([some a_to_g sites permutation],[some c_to_t sites permutation])
    editing_permutarions = find_permutations_in_range(creat_editing_sites_permutations_one_type(a_to_g_sites),creat_editing_sites_permutations_one_type(c_to_t_sites))
    
    
    
    
    
    #for each pemutation produce new sub_sequence (with window range) and check for new cleavage sites relateive to whole sequence
    for perm in editing_permutarions:
        print('\n')
        new_seq = edit_rna_as_peptide(sequence[window_range[0]:window_range[1]+1],window_range, perm[0], perm[1])
        new_cleavage_sites = [x+window_range[0] for x in find_cleavage_sites(regex,new_seq)]
        print(new_seq)
        print('permutation:' + str(perm))
        print('new cleavage sites:' + str(new_cleavage_sites))
        
        #the relevant cleavage sites for creating peptides in original cleavage sites window.
        if i == None: #we are at the end of the sequence and last original cleavage site was already appended to cleavage_sites
            relevant_cleavage_sites = deque([x for x in new_cleavage_sites if x>=cleavage_sites[0]])
        else:
            relevant_cleavage_sites = deque([x for x in new_cleavage_sites if (x>=cleavage_sites[0] and x<=cleavage_sites[-1])])
        if cleavage_sites[0] == 0: #we are at the beginning of the sequence and first coordinate to cleave from is 0
            relevant_cleavage_sites.appendleft(0)
        print('relevant cleavage sites:' + str(relevant_cleavage_sites))
        
        n+=1
    
    
    
    
    if cl < ml:
        cl += 1
    for j in trange[:cl-1]:
        print('j='+str(j))
        seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
        print('peptide coordinates: ' + str(cleavage_sites[j]) + ',' +str(cleavage_sites[-1]))
        print(seq)
        if seq:
            if min_length is None or len(seq) >= min_length:
                peptides.append((seq,cleavage_sites[j],cleavage_sites[j]+len(seq)-1))


def find_cleavage_sites(regex,sequnce,cleave_before_pattern=False):
                    
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
            

"""
change a sub-sequence given:
the subsequence coordinates relative to the sequence
a_to_g editing coordinates relative to the sequence
c_to_t editing coordinates relative to the sequence
"""
#k = edit_rna_as_peptide(mini_seq, a, lst1, lst2)
def edit_rna_as_peptide(sequence, sequence_coordinates, a_to_g_list, c_to_t_list):
    splited_ag = [sequence[i:j] for i,j in zip([0] + [x+1-sequence_coordinates[0] for x in a_to_g_list], [x-sequence_coordinates[0] for x in a_to_g_list] + [None])]
    ag_swaped_seq = 'G'.join(splited_ag)
    splited_ct = [ag_swaped_seq[i:j] for i,j in zip([0] + [x+1-sequence_coordinates[0] for x in c_to_t_list], [x-sequence_coordinates[0] for x in c_to_t_list] + [None])]
    new_seq = 'T'.join(splited_ct)
    return new_seq

        
def find_permutations_in_range(a_to_g_sites,c_to_t_sites):
    return[(l[0],l[1]) for l in it.product(a_to_g_sites,c_to_t_sites)]
    
    
def creat_editing_sites_permutations_one_type(lst):  
    editning_perm = []
    for l in range(len(lst)+1):
        [editning_perm.append(list(subset)) for subset in it.combinations(lst, l)]
    return editning_perm



def update_editing_sites_in_window(window_range,editing_sites_list,editing_sites_it,site=None):
    
    #cleane left side sites if out of range
    if len(editing_sites_list):
        k = editing_sites_list[0]
        while k < window_range[0]:
            editing_sites_list.popleft()
            k = editing_sites_list[0]
    
    #append current site if in range
    if site != None and site > editing_sites_list[-1]:
        if site <= window_range[1]:
            editing_sites_list.append(site)
    
    #keep appending sites in range until reaching an out of range site
    try:
        for next_site in editing_sites_it:
            if next_site <= window_range[1]:
                editing_sites_list.append(next_site)
            else:
                break
        #return current site reached by the sites iterator and all a list of all editing sites in range
        return next_site, editing_sites_list
    
    except:
#        print('iterator exhusted')
        return None, editing_sites_list
    


        
#peps = cleave_rna_as_peptide(seq,a,1)
def cleave_edited_rna_as_peptide(sequence, start_coor,cleavage_sites_list, a_to_g_sites, c_to_t_sites, missed_cleavages=0, min_length=None):
    
    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
        
    for i in it.chain(cleavage_sites_list):
        cleavage_sites.append(i)
        print(cleavage_sites)
        if cl < ml:
            cl += 1
        
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]-start_coor:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append((seq,cleavage_sites[j],cleavage_sites[j]+len(seq)-1))
                        
    return peptides



#a_list, c_list = c_a_list(sequence)
def c_a_list(sequence):
    
    a_list = []
    c_list = []
    
    for i in range(len(sequence)):
        if sequence[i:i+1] == 'A':
            a_list.append(i)
        elif sequence[i:i+1] == 'C':
            c_list.append(i)
            
    return a_list, c_list

# =============================================================================
# cleavage sites window is: deque([21, 33, 54], maxlen=3)
# window coor: (18, 56)
# a to g sites are deque([18, 19, 20, 30, 38, 45])
# c to t sites are deque([23, 25, 41, 51, 53])
# 
# 
# AGGCCTATTATCAGGCATTCAGGTATCGGCATGCGCCCT
# permutation:([19, 20, 45], [23, 25, 41])
# new cleavage sites:[33, 54]
# relevant cleavage sites:deque([33, 54])
# 
# def cleave_edited_rna_as_peptide(sequence, start_coor,cleavage_sites_list, a_to_g_sites, c_to_t_sites, missed_cleavages=0, min_length=None):
#     
#     peptides = []
#     ml = missed_cleavages+2
#     trange = range(ml)
#     cleavage_sites = deque([0], maxlen=ml)
#     cl = 1
#         
#     for i in it.chain(cleavage_sites_list,[len(sequence)+start_coor]):
#         cleavage_sites.append(i)
#         print(cleavage_sites)
#         if cl < ml:
#             cl += 1
#         
#         for j in trange[:cl-1]:
#             seq = sequence[cleavage_sites[j]-start_coor:cleavage_sites[-1]]
#             if seq:
#                 if min_length is None or len(seq) >= min_length:
#                     peptides.append((seq,cleavage_sites[j],cleavage_sites[j]+len(seq)-1))
#                         
#     return peptides
# 
# =============================================================================



