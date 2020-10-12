import re
import sys
import itertools as it
from rna_peptide import rna_peptide
from collections import deque
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis

all_mm = ['AG','AC','AT','CA','CG','CT','GA','GC','GT','TA','TG','TC']
empty_dict = {'AG':[],
              'AC':[],
              'AT':[],
              'CA':[],
              'CG':[],
              'CT':[],
              'GA':[],
              'GC':[],
              'GT':[],
              'TA':[],
              'TG':[],
              'TC':[]}


def create_edited_rna_peptides(sequence,sites_dict,digestion_rule, missed_cleavages = 0,
                               look_ahead_for_editnig=3,look_behind_for_editing=6,min_length = None, max_mass = None, max_sites_per_pep = None):
    
    
    rna_peptides_list = []
    under_min_length = []
    over_max_sites_peps = []
    over_mass_peps =[]
    over_mass_pep_combs_dict = {}
    
    seq_length = len(sequence)
    cleavage_sites, binary_for_fixed_cs = create_all_cleavage_sites(sequence,digestion_rule,sites_dict,look_ahead_for_editnig,look_behind_for_editing)
    cleavage_sites_lengeth = len(cleavage_sites)
    original_sites, original_sites_boo = create_all_cleavage_sites(sequence,digestion_rule,empty_dict,look_ahead_for_editnig,look_behind_for_editing)


    #running over all possible windows from a site until mc+1 fixed sites ahead
    for i in range(cleavage_sites_lengeth-1):
        fixed_sites_cnt = 0
        j=1
        
        #running over all j (1+j -> index of cleavage site at the end of the subsequence window) untill betwin i-th cs and (i+j)-th cleavage site there are mc number of fixed site or i+j site is sequence end
        #in this method all possible windows for subsequences are examined
        while (fixed_sites_cnt < missed_cleavages+1) and (cleavage_sites_lengeth-i-1>=j):    
            enter_seq = True
            examin_each_comb_mass = False
            edge = None
            append_zero = False
            append_seq_end = False
            subsequence_coor = (cleavage_sites[i],cleavage_sites[i+j]-1)
            coor_str = '_'+str(subsequence_coor[0])+'_'+str(subsequence_coor[1])
            subsequence = sequence[subsequence_coor[0]:subsequence_coor[1]+1]
            
            if binary_for_fixed_cs[i+j]: #count site if site ahead is fixed.
                fixed_sites_cnt+=1
            
            if min_length != None: #rule out peptide if shorter than min_length
                if len(subsequence) < min_length:
                    enter_seq = False
                    under_min_length.append(subsequence_coor)
        
# =============================================================================
            #determining range for editing perutations:
            if binary_for_fixed_cs[i]: #site is fixed and therfore upstream editing dosent matter
                edit_range_start = cleavage_sites[i]
            else:    
                if cleavage_sites[i]-look_behind_for_editing < look_behind_for_editing: #sequence is close to beginning, looking for permutations from beginning
                    edit_range_start = 0
                else:
                    edit_range_start = cleavage_sites[i]-look_behind_for_editing
             
            if binary_for_fixed_cs[i+j]: #site is fixed and therfore downstream editing dosent matter
                edit_range_end = cleavage_sites[i+j]
            else:    
                if cleavage_sites[i+j]+look_ahead_for_editnig > seq_length: #sequence is close to end, looking for permutations till end
                    edit_range_end = seq_length
                else:
                    edit_range_end = cleavage_sites[i+j]+look_ahead_for_editnig
             
            permutation_range = (edit_range_start,edit_range_end-1) #this is the range in which examine editing combinations (could extand further than subsequence (peptide) coordinates)
# =============================================================================
            
            # getting all editing sites in permutations window
            sites_in_perm_range = {}
            [sites_in_perm_range.update({mm:find_sites_in_window(permutation_range,sites_dict[mm])}) for mm in all_mm]
            number_of_sites_in_perm_range = sum([len(sites_in_perm_range[mm]) for mm in all_mm])

            #if too many editing sites in peptide, dont examine peptide
            if max_sites_per_pep != None:
                if number_of_sites_in_perm_range > max_sites_per_pep:
                    enter_seq = False
                    over_max_sites_peps.append(subsequence_coor)
                    
                    
            #last condition to rule out peptide - lower bound mass > max_mass allowed
            #getting all editing sites in susequence window (only for getting maximal\minimal mass range)
            sites_in_subseq_range = {}
            [sites_in_subseq_range.update({mm:find_sites_in_window(subsequence_coor,sites_dict[mm])}) for mm in all_mm]
            minimal_mass, maximal_mass = get_mass_range(subsequence,subsequence_coor,sites_in_subseq_range)
            
            if max_mass != None:
                if minimal_mass > max_mass:
                    enter_seq = False
                    over_mass_peps.append(subsequence_coor)
                    #if max_mass is betwin peptide minimal mass and maxiamal mass each editing comb mass is calc before appended
                elif maximal_mass > max_mass: 
                    examin_each_comb_mass = True
# =============================================================================
            if enter_seq:    
                sebsequence_in_perm_range = sequence[edit_range_start:edit_range_end] #on this subsequence we examine different combinations (could extand further than subsequence (peptide) coordinates)
                
                permutations_dict_perm_range = {}
                [permutations_dict_perm_range.update({mm:creat_editing_sites_permutations_one_type(sites_in_perm_range[mm])}) for mm in all_mm]
                #all eiditing sites permutation in permutation_range
                permutations = find_permutations_in_range(permutations_dict_perm_range)
                
                #determining if peptide sequence is at the edge of full sequence
                if edit_range_start == 0:
                    append_zero = True
                    edge = 5
                if edit_range_end == seq_length:
                    append_seq_end = True
                    edge = 3
                
                for perm in permutations:
                    
                    #edited peptide including up\downstrame relevant codons according to cleavage rule
                    extended_edited_sebsequence = edit_rna_as_peptide(sebsequence_in_perm_range,permutation_range,perm)
 
                    #the (edited or not) cleavage sites withing the range of subsequence (including_bounderies)
                    perm_cs = [k+edit_range_start for k in find_cleavage_sites(digestion_rule,extended_edited_sebsequence,append_zero=append_zero, append_seq_end=append_seq_end) 
                                if (k+edit_range_start>=cleavage_sites[i] and k+edit_range_start<=cleavage_sites[i+j])]
                    
                    #the original cleavage sites within the peptide range
                    original_sites_in_pep_range = [k for k in original_sites if (k>=subsequence_coor[0] and k<=subsequence_coor[1])]
                    
                    #if cleavage site at start is fixed, append to cs list for perm
                    if edit_range_start == cleavage_sites[i] and edit_range_start not in perm_cs:
                        perm_cs.insert(0,edit_range_start)
                    #if cleavage site at end is fixed, append to cs list for perm
                    if edit_range_end == cleavage_sites[i+j] and edit_range_end not in perm_cs:
                        perm_cs.append(edit_range_end)
                     
                     
                    #if bouderies still exists as cleavage sites and no more than mc cleavage sites in perm_cs - a valid seq
                    if set([cleavage_sites[i], cleavage_sites[i+j]]).issubset(perm_cs) and len(perm_cs)<=missed_cleavages+2:
                        
                        #translating subsequence according to sites from perm that are within permutations range
                        perm_sites_in_subseq_range = ()
                        for k, mm in enumerate(all_mm):
                            perm_sites_in_subseq_range+=(find_sites_in_window(subsequence_coor,perm[k]),)
#                        [perm_sites_in_subseq_range+=(find_sites_in_window(subsequence_coor,perm[i])) for i, mm in enumerate(all_mm)]
                        edited_sub_seq = edit_rna_as_peptide(subsequence,subsequence_coor,perm_sites_in_subseq_range)
                        
                        #determine if cleavage sites defining peptide were created due to editing
                        N_terminus = 'no_change'
                        if cleavage_sites[i] not in original_sites:
                            N_terminus = 'new_cleavage_site'
                        C_terminus = 'no_change'
                        if cleavage_sites[i+j] not in original_sites:
                            C_terminus = 'new_cleavage_site'
                        cancelled_original_cs  = False
                        if not(all(cs in perm_cs for cs in original_sites_in_pep_range)):
                            cancelled_original_cs = True
                            

                        if not examin_each_comb_mass: #no need to examine each edited peptides mass - append peptide to peptides list
                            rna_pep_perm = rna_peptide(edited_sub_seq,subsequence_coor, extended_edited_sebsequence, permutation_range,edited_sites = perm, cleavage_sites = perm_cs, edge = edge, N_terminus = N_terminus, C_terminus=C_terminus, sites_in_pep_range = number_of_sites_in_perm_range, cancelled_cs = cancelled_original_cs)
                            rna_peptides_list.append(rna_pep_perm)
#                            rna_pep_perm.print_peptide_data() 
                        
                        else: #each peptide mass should be examined before appended to list
#                            print(''.join(str(Seq(edited_sub_seq, generic_rna)).split('*')))
                            comb_mass = ProteinAnalysis(''.join(str(Seq(edited_sub_seq, generic_rna).translate()).split('*'))).molecular_weight()
#                            print(perm)
#                            print(comb_mass)
                            if comb_mass <= max_mass:
#                                print('in')
                                rna_pep_perm = rna_peptide(edited_sub_seq,subsequence_coor, extended_edited_sebsequence, permutation_range, edited_sites = perm, cleavage_sites = perm_cs, edge = edge, N_terminus = N_terminus, C_terminus=C_terminus, sites_in_pep_range = number_of_sites_in_perm_range, cancelled_cs = cancelled_original_cs)
                                rna_peptides_list.append(rna_pep_perm)
#                                rna_pep_perm.print_peptide_data()
                            else:   #peptide mass is out of range - append to designated dictionary
                                if coor_str in over_mass_pep_combs_dict:
                                    over_mass_pep_combs_dict[coor_str].append(perm)
                                else:
                                    over_mass_pep_combs_dict.update({coor_str:[perm]})

# =============================================================================
            j+=1

    return rna_peptides_list, under_min_length, over_max_sites_peps, over_mass_peps, over_mass_pep_combs_dict


"""
given editing coordinations and a sequence and a lookahead lookbehind range
run over all codons, edit and if codos could be a cleavage site append to list.
create another binary list - for each codon fixed or mutable
"""
def create_all_cleavage_sites(sequence,digestion_rule,sites_dict,look_ahead_for_editnig,look_behind_for_editing):
    
    cleavage_sites = deque([0])
    binary_for_fixed_sites = deque([1])
    length = len(sequence)
    
    #checking that sequence is devided to whol codons
    if len(sequence)%3:
        return None, None
    else:   
        for i in range(0,length-2,3): #run over all codons not including last wich is a fixed site
            
            #sequence is close to beginning, looking for permutations from beginning
            if i <= look_behind_for_editing:
                coor = (0,i+look_ahead_for_editnig-1)
                codon_to_asses_in_window = i
            #sequence is close to end, looking for permutations till end
            elif length-i<= look_ahead_for_editnig:
                coor = (i-look_behind_for_editing,length-1)
                codon_to_asses_in_window = look_behind_for_editing 
            else:
                coor = (i-look_behind_for_editing,i+look_ahead_for_editnig-1)
                codon_to_asses_in_window = look_behind_for_editing
            
            seq = (sequence[coor[0]:coor[1]+1])
            
            #all editing permutations in range
            perm_dict = {}
            [perm_dict.update({mm:creat_editing_sites_permutations_one_type(find_sites_in_window(coor,sites_dict[mm]))}) for mm in all_mm]
            permutations = find_permutations_in_range(perm_dict)
            
            is_cs, is_fixed = determine_site(seq,coor,digestion_rule,codon_to_asses_in_window,permutations)
            
            if is_cs:
                cleavage_sites.append(i)
                binary_for_fixed_sites.append(is_fixed)
                
        cleavage_sites.append(length)
        binary_for_fixed_sites.append(1)
        
        return cleavage_sites,binary_for_fixed_sites
            

def determine_site(subsequence,coor,digestion_rule,codon_to_asses,permutations):
    """
    subsequence - expected a codon with some up\downstream relevant nucs according to cleavage rule
    determine if codon is a cleavage site (1/0) and if it is mutable by some editing (1/0)
    """  
    
    cleavage_site = set()
    
    for perm in permutations:
        seq = edit_rna_as_peptide(subsequence,coor,perm)
#        print(seq)
        cleavage_sites = find_cleavage_sites(digestion_rule,seq)
#        print(cleavage_sites)
        
        if codon_to_asses in cleavage_sites:
            cleavage_site.add(1)
        else:
            cleavage_site.add(0)
        
        if len(cleavage_site) == 2:
            break
        
    if len(cleavage_site)==2:
        return 1,0
    elif list(cleavage_site)[0] == 1:
        return 1,1
    else:
        return 0,1
       
            
def find_cleavage_sites(digestion_rule,sequence,append_zero=False, append_seq_end = False):
    """
    find all in-frame cleavage sites given a sequence and a digestion_rule
    also append 0 and sequence_length if needed
    """
    
    if append_zero:
        initial_list = [0] + [x.end() for x in digestion_rule.finditer(sequence)]
    else:
        initial_list = [x.end() for x in digestion_rule.finditer(sequence)]
    initial_set = set(initial_list)

    for site in initial_list:
        match1 = digestion_rule.match(sequence,site-1)
        match2 = digestion_rule.match(sequence,site-2)
        if match1:
            initial_set.add(match1.end())
        if match2:
            initial_set.add(match2.end()) 
 
    #return only sites representing in-frame codons
    cleavage_sites = [x for x in initial_set if not x%3 and x!=len(sequence)]
    
    if append_seq_end:
        cleavage_sites.append(len(sequence))
    
    cleavage_sites.sort()
    return cleavage_sites
            


#k = edit_rna_as_peptide(tsubseq, [273,314], [305,306], [])
def edit_rna_as_peptide(sequence, sequence_coordinates, permutation):
    """
    change a sub-sequence given:
    permutation is a nested list of all editing sites, 
    ****sites lists order is identical to the order in all_mm list****
    for each mm, get list of sites from permutation and swap all sites with second letter in missmatch
    """
    new_seq = sequence
    for i,mm in enumerate(all_mm):
        new_nuc = mm[1]
        sites = sorted(permutation[i])
        splited_by_old = [new_seq[i:j] for i,j in zip([0] + sorted([x+1-sequence_coordinates[0] for x in sites]), sorted([x-sequence_coordinates[0] for x in sites]) + [None])]
        new_seq = new_nuc.join(splited_by_old)

    return new_seq

  

def find_permutations_in_range(permutations_per_mm):    
    """
    return all combinations of elements from 12 lists of mm types
    each permutation is a tuple of 12 lists with the miss matches present that permutation by the order of all_mm
    """
    return[(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11]) for l in it.product(*[permutations_per_mm[mm] for mm in all_mm])]
    
   
"""
return all permutation of a list
"""
def creat_editing_sites_permutations_one_type(lst):  
    editning_perm = []
    for l in range(len(lst)+1):
        [editning_perm.append(list(subset)) for subset in it.combinations(lst, l)]
    return editning_perm



"""
return a list of all elements from a given list that are within a given range
"""
#a = find_sites_in_window((18,56),[1,2,4,5,16,17,29,30,31,54])
def find_sites_in_window(editing_window,editing_sites):
    editing_sites_in_window = []
    editing_sites_it = it.chain(editing_sites)
    for i in editing_sites_it:
        if i<=editing_window[1]:
            if i >= editing_window[0]:
                editing_sites_in_window.append(i)
        else:
            break
    return editing_sites_in_window


def get_mass_range(sequence,coor,sites_in_subseq_range):
    
    aa_seq_of_max_mass = []
    aa_seq_of_min_mass = []

    for i in range(0,len(sequence),3):
        
        max_aa_mass = 0
        min_aa_mass = 99999
        
        codon = sequence[i:i+3]
        editing_window = [i+coor[0], i+2+coor[0]]
        codon_perm_dict = {}
        [codon_perm_dict.update({mm:creat_editing_sites_permutations_one_type(find_sites_in_window(editing_window,sites_in_subseq_range[mm]))}) for mm in all_mm]
        permutations = find_permutations_in_range(codon_perm_dict)

#        print('i=' + str(i))
        for perm in permutations:
            edited_sub_seq = edit_rna_as_peptide(codon,editing_window,perm)
            
            aa = str(Seq(edited_sub_seq, generic_rna).translate())
            aa_mass = ProteinAnalysis(aa).molecular_weight()

            if aa_mass>max_aa_mass:
                max_aa_mass = aa_mass
                heavy_aa = aa
            if aa_mass<min_aa_mass:
                min_aa_mass = aa_mass
                light_aa = aa
                
        aa_seq_of_max_mass.append(heavy_aa)
        aa_seq_of_min_mass.append(light_aa)
        
    maximal_mass = ProteinAnalysis(''.join(aa_seq_of_max_mass)).molecular_weight()
    minimal_mass = ProteinAnalysis(''.join(aa_seq_of_min_mass)).molecular_weight()
    
    return minimal_mass, maximal_mass


#sequence = 'ATGGCTGCTTTGAGACAGCCCCAGGTCGCGGAGCTGCTGGCCGAGGCCCGGCGA'
#sites_dict = {'AG':[0,12,14,31,53],
#              'AC':[],
#              'AT':[],
#              'CA':[],
#              'CG':[],
#              'CT':[4,7,19,20,40],
#              'GA':[],
#              'GC':[],
#              'GT':[],
#              'TA':[],
#              'TG':[],
#              'TC':[]}


#sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/' +'test3.txt', 'w')
#sys.stdout = open('C:/Users/shosh/OneDrive/Desktop/' +'new.txt', 'w')

#input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
##input_fasta = 'sim.fasta'
#input_fasta = 'in_frame_rna_from_orfs_squ.fasta'
#digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
#sequence = 'ATGGCTGCTTTGAGACAGCCCCAGGTCGCGGAGCTGCTGGCCGAGGCCCGGCGAGCCTTCCGGGAGGAGTTCGGGGCCGAGCCCGAGCTGGCCGTGTCAGCGCCGGGCCGCGTCAACCTCATCGGGGAACACACGGACTACAACCAGGGCCTGGTGCTGCCTATGGCTCTGGAGCTCATGACGGTGCTGGTGGGCAGCCCCCGCAAGGATGGGCTGGTGTCTCTCCTCACCACCTCTGAGGGTGCCGATGAGCCCCAGCGGCTGCAGTTTCCACTGCCCACAGCCCAGCGCTCGCTGGAGCCTGGGACTCCTCGGTGGGCCAACTATGTCAAGGGAGTGATTCAGTACTACCCAGCTGCCCCCCTCCCTGGCTTCAGTGCAGTGGTGGTCAGCTCAGTGCCCCTGGGGGGTGGCCTGTCCAGCTCAGCATCCTTGGAAGTGGCCACGTACACCTTCCTCCAGCAGCTCTGTCCAGACTCGGGCACAATAGCTGCCCGCGCCCAGGTGTGTCAGCAGGCCGAGCACAGCTTCGCAGGGATGCCCTGTGGCATCATGGACCAGTTCATCTCACTTATGGGACAGAAAGGCCACGCGCTGCTCATTGACTGCAGGTCCTTGGAGACCAGCCTGGTGCCACTCTCGGACCCCAAGCTGGCCGTGCTCATCACCAACTCTAATGTCCGCCACTCCCTGGCCTCCAGCGAGTACCCTGTGCGGCGGCGCCAATGTGAAGAAGTGGCCCGGGCGCTGGGCAAGGAAAGCCTCCGGGAGGTACAACTGGAAGAGCTAGAGGCTGCCAGGGACCTGGTGAGCAAAGAGGGCTTCCGGCGGGCCCGGCACGTGGTGGGGGAGATTCGGCGCACGGCCCAGGCAGCGGCCGCCCTGAGACGTGGCGACTACAGAGCCTTTGGCCGCCTCATGGTGGAGAGCCACCGCTCACTCAGAGACGACTATGAGGTGAGCTGCCCAGAGCTGGACCAGCTGGTGGAGGCTGCGCTTGCTGTGCCTGGGGTTTATGGCAGCCGCATGACGGGCGGTGGCTTCGGTGGCTGCACGGTGACACTGCTGGAGGCCTCCGCTGCTCCCCACGCCATGCGGCACATCCAGGAGCACTACGGCGGGACTGCCACCTTCTACCTCTCTCAAGCAGCCGATGGAGCCAAGGTGCTGTGCTTG'
#sequence = 'AACATGAAAGGTGACATGTTTTGTAGAAAGAAAAAACATAAAGAATGGCTAAAAATGGCTAAAGCTGGCCAATACACCTTAAAAGAGAAAATTCAAGAAGAGTTTTTGGAATGCAAAATCTGTTTTGAACCGTATGTAAAACCGAAGGCATTACCTTGCCTTCACTCATTCTGTGCCGAGTGTCTAAAGGACTACGTGCGAAAGAACCCCAACAAAAATGCAGTGCGTTTTTGCTGTCCAATCTGTCGCAAAGAAATCCCGATGCCTGCCGGTGGCATCGACGATTTCCAAGATAATTTTTGGTTATTGAGCTTGTCAAATTCCTTGGAAGAAGGAGAGGAAGACTGTCGTGTATCGTGCAATGGGAAGGCCGTGCGTAGCAATGCCTGGCCTACTCCAAAATGGAAACAACCCAAGAATGTTTCAACCTCATCTAAACCTTTGTATCCCCCATCATTCCAGGAACTTAAACTCAGGGGGTTTGAATGGTATTTTGGTAAAGTGAGTCGTAATGCTTCAGAAGAATGGCTTCTCCATCCTGGGCTACAGAAGGGAACATTCCTTATTCGACAAGGGGAGGCTCTTCCTGACACGTATACCTTGTCAGTTCGTGACTGTGACGAGCTGAGGGGTTACCTGGTGAAACATTACAAGATTCTAACCAAAAAAGCAAGTGATGGGGAGAAGGAAGTTTATTACATCACCCCAAAGCGAACTTTCCGTTCGCTTGAGGAACTGGTGAATCATTATTCAATTTCAGACGGCCTTTGCTGTAAACTAACACAAATATGTAACAAACCAAGGTCTTTACTCTGGGCAATGGAAAGAGGAAAACCGGATGATTTTATGACCACTAAGGACACTTTGCAATTGGTCAAGAAAATTGGCAGTGGCCAATTTGCTGAAGTTTATTATGCCAAATGGAACAACCAAGTAGAGGCAGCCGTAAAGATGCAAAAAAAGGATTGCGTGACAACCTCGGCTTTCTTAGATGAAGCTCAGATCTTAAAGACAATTCAACATGTAAATATAATCAAGTTGTTGGCTGTGTGCAGTGACGAACCTGTTTATTTGGTGACAGAATATATGCCTAATGGCCGACTCTCACAATACCTTCGAGAAGGAAAGGGGAAACAGCTCGGGGTTAACAGTCTCCTCTGGCTTGCGGCTCAGATCGCAGACGGAATGGCTTACATGGAAAAGGAGAATTTTGTTCATAGAAATTTGGGTGCCCGAAACATTCTGGTTGCTGACCAAAACAAAGTGAAAATTGCTGGTTTTGGGATGACAAAAGTGGCCGATGATCCTGATTACAATTTCAGAAAAGGTTTGAAAATGGCTGTAAAATGGATGGCCCCTGAAGTGTTGTTGTACAATAAATACAGCACAAAGGCTGACGTGTGGTCGTTTGGAATTGTCCTAATGGAAATATTTTCATATGGCAAAGAACCCTACGATGGTATGGGAAGCAAGGAAGCATTTGAAAATGTCCAGTCAGGTTACCGAATGCCTTGTCCTCACTGTTGCCCTGCCGAGGTATACAATGTGGCACTGACCTGTTGGAATATCAACACCCAGCGTCGGCCATCTTTTGACTTTCTTAACAGCTTCCTTCATGACTGGCACTATACGTCC'
#sequence = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
#sequence = str(rec[0].seq)
#missed_cleavages = 3
#min_length = 6
#max_mass = 4600
#a_to_g_editings = a_list
#c_to_t_editings = c_list
#look_ahead_for_editnig = 3
#look_behind_for_editing = 6
#max_sites_per_pep = 20
