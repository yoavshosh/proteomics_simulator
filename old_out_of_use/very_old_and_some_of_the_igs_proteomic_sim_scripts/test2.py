import re
import itertools as it
import sys
from rna_peptide import rna_peptide
from collections import deque
from bisect import bisect_left
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_rna


sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/' +'test.txt', 'w')

digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
sequence = 'AAACGCTAGAAACCATGGAAACCCACTATCAGGCATTCAGGCATCAGCATGCGCCCTATTGCACGC'
#sequence = 'ACTACCACCACTACTACTACTACACAAATTTTCTATGGCCATGTAATCATCATGTAGTAAGTAGTGAGAAGAGATGGTGCGGCAAACGAGGTATCTCCTGGTTGGGAATATACCCGAAAAAATCACCAAGGAGAAAATCACAGAACACTTTAAGAACTATGGAAAGATCCAGAATGTAAAGTTTCATCCAAAAAAAGAAAATGAACTTGGCATTACAGCAACTGTGGCATTTATGGACATTAAATGTGCCTCTAAGGCATATAGTTCAGAGAACAAAATAGAAAATGTAGTGCTTCGGACAGAATACAGTGAGGGCACTGCAACAGGAAGTGTTGTTACAAGACCTGTGGACCCTGGATTGCCTACAACTGGAGCCCACAGGGGACAGTACCTACATGCTCGAGGGCCATCTACTAGCTTCACAACACGGACTAAAGGTAGTCTTGATGAAGATGATACATACAGGGAATATTACTATGGTGGTAGCGGCAGCACTGGCAGTAGGGAAAACAGTCAATTTGAAAACAGATATCCCACATATGATGAACAAGTACAGCATCAGTCAAGGTCTTATCAAAGAAATAGGACCTCATTTTCCAAGCAGTATGAGCGAGGACCACAACGACAGAAGTCTAGCAGTGCGAGTGTCACGGGTGTTGTCACTAGCGGTGGTAGTGGCAGTAGTGTCACAGGGGTTGGGAATAGCAATACAGTTGGCAGTACAGTGACACCCAGCAGTGGTGGTAGTGTGGTTGGGAGTGGAGCTGTTGGAACCAGTGGCAGTGGGGGTGGGGGTGCTGGTAGTGGTGCTGGCAGTGGAGGAGGCAACAATAGTGGAGGTGGTGGCACCAACAATGGAGGGGACACTGGCAACACTGGGGTCAGTGTTGTGGTCTCTAATACCTTGCAAAAGCATTCCTTTGATGAGCGTTACCACGAACAGTATAGCAGCGACAGAGAAAGTTTTGAAAGGCCACGGATATCCACATTGTCTCCAACATGTTCCATCTCTCAGGCCAGCAGCCGACACAGTAGTGGTCGGCGGAAAGGTTCTGTCTCTGTGTCCCGTTCGTCCAGTCGCAGTCGGTCAAGGTCCAGTTCCAGCAGTGAGTATAGTCGGTCTTCATCAGGAAGCAGTGGATCAAGATCACATTCCAGAAGTCGTGCAAGCAGCCGGGGGAGTTTGCGGCGAATACCTGTTGGTTCAGCATCGAGTCAAAAACACTCTAAATCAGAAAGTCAAAGTATGATGTCCAGTAAATCCTCAACAGCAGGATCGCAATCGGAAAGAAAACCACTTGGGATTTGTGTAACGGGTTTACCATTAAGATCCAGTGATACTAGTCTTCGAGATGGTCTTTTTCATGAATACAAAAAACATGGAAAAGTTACTTCTGTACAGGTGCTAGGACAAGGCGAAGAAAGATATTCAGTGGTCAGTTTCCGGAAACCTGAAGATGCGGCAAAGGCATTAGAAGCTTCACAGGACAAAATGTTTTTTGGAAAGAAAATCAAAGTTATGGCTCATGAAGGAGTTGAAGTGGAAGACAATGAGTTCCGTCCCCCAGAAGCTGAATGTGATGAATACCACCCTAAGGCAACACGAACGTTATTTGTGGGAAACTTAGAAAAAGAAATAACGACACAAGAGCTTAGAGATCGTTTTAAGTCTTTTGGAGAAATTATAGATATTGATATCAAACGTCAGGGTGCAGTGTCAGCATATGCATTTGTACAGTATGCAGATATTACAAGTGTTGTAAAATCACTCCGAAAAATGGAGGGCGAACATATAGGAGCCAATAAAATAAAGCTTGGTTTTGGTAAAAGTATGCCAACCAATTGTGTGTGGCTCGATAACCTGTCTGAATCCGTAACAGAGAAGTTTCTTTGCCAACAATTTGGCCGTTATGGCGAGGTCAGTCATGCAGTGATTGATCGGTATAAAGGAAGGGCCTTGATATATTTCACAAGTATGGATACAGCTCAATATGCAGTCACTGAGATGAGAAATAGAATTTTGAAGAAGAAACGAGTGCAGACTGGTGACTTTAGACTCGGTGGACATCTGAATATTGATTTTGCTAGTCGAGATTGTCAAACAGCTTTCTTTGAAAAAATGGAGCGGACGGGTCAACTGCAACCTGGAGAACGTCCTGATGAACGAGGAGGGCGCCTCTACCGGCACAATAGGGGAGCCAATTTCGAGCAATATTATGATTCTCCTGCAACCCCTTCCGCTTACAAGGAAGATAGCGCAAGATACGAGACTGGCACTGTAACACCTGTTGCTACAACTGCACGAGTTGGACTGACCCCACCACAGTTCAGTACGACAAAACGTAATCGAGCAACAAATTTTCAAGCCTCTTCGACCAGACAAGGATCCTCATACAGAGGAAGACATTTTACAGAAAGTACTTACTCAGAAGAATTCCCGGCAAACCAGAGACAGAGACATATTGATGAATACAGTCAAGGAAGTGCTCCTCCTTATGCTGAGGAAGATTCCTATGAACAGGAGCTGAGAGAATATGGCTACAGTCAACGGGAACGCCGTGAGAGGGGTAACGTAACTCCACCGAGACGATACCACAGTCCTTATCGGGAAGATAGGCATTCAAATGCCGACAAAGAATCGACGCGGGGTAGCCGTTATGAAGAAGCCTTCACCAAAGACAGATATGATTATTACAGTAGTGACAGTCCTGTCGTCAGTGACAATGGTCATGATGAGAACTGGGTGGATCATCATACATTGGATCTTTCAGATCGTCATGCCACACCAACAAAAATGAAATATTCAAGTCGGAACGAGACACGCTTAAAGCAATCGTCTCGGGATCACAGCCCTGCCATTGAGTTTTCCTACAGGCACTCACCAGTGAGACGAAGAAGCAAGAGCAGGAGCAGAAGTAAAAGTGGGAGTCGGACAGAAAGTCCTAGTCGAAGTAGGAGCCGGACGCCTCTAAGGGAAAGGCGGAAAGAATCTGTGTTCCGGCCACGTTCTCCCCACACACCACAAACACCTCCTACACCAACAGAAGAAGAACCTCCCATTTATTCTAAAATCCGATCAATCAGCCCCACACCCCGGTTTGAGGAATTTGAGCGAACAGGAGAATTGAAAAAAATATGCAGGACAGATGATAAAATTGAAGGCATTAGGTCTTCGGATCATGAAAAGGTGTACAGTATAAAATGCTTGGGTACTCCATCTATGTCTTTCATGAAAGGGGAAACGACACAAGCCCGTAAAAAGAGTAGGGAGTTTGATGATAGCAAACATGATGTAATAGCAGGACGTCTGCATAAGTCCCAGGTGCCTGCAACGAACCGGCTGTTAGAATCTGCGCTCCAAGACAGACATCAAGTCTTGTTAAAGAAGGCAGATGAAAACCGTCTCAATAGACGCAAACTATTAGACCAGCCACATGGTAGCTCGAGTCATTCCAGTGGCCACGATGATTCAGATCCTGGTTCTTTTGCTGAGACTGATCTTGGCCATTTGCATCGAGAAAAGCGACTACTCCTAGAAAAATTAAAACAACTGGAAGACGCTGGAAGCCCCAGTGATAATGATTCTTTACCGACACATGATGATCTGCGAGAAGATGGACACGGGCGATTTGTACCGAGGAAAGGTCGACCTGAGGATCCTCATTGTAGCAGTTTACGGCATAGTCTGAGCGACATACACGCAAAGGAAAAACGAAAACTAGAAACTTCTCAGAGTTTTCGAAAACAGATGGAAGCCAGGCGCCTGTTAGAACAGACTCAGGGAGAGTCAAATTTGTTATCAAAACCAAAAAAGAGTGTTGAAAAAGTTGTGATGCTTGATGGTATGGATAGCGAAGGAGAGGAGGAAATAAGTATTACCTCTCCGAAAGCAATCCCCACTGCATTCAAGCCAAAGAGAAAGAAGAAACTCGAAGAAGACCCAGTGTCTGGTGTTCGTACTGGCAGGTTTTATCGCATGAAGAAAAGTGGAATGATGAGCAGTGATGACGAAAAAGATCATTCAAAATCTCTCGAAGCAGAAGAACCTGTTATGCGACGCACCAGTCGCGAGAGCAAAGACGAACTATCCGCTGTTATAAGCCGGTTTAAGGAGCAGTCAAGACATTCTGAATTAAGCCGGGAGGTTTCCGAGTCTGGGTCTCACAAGGAATCCAGAGAGTATAAAATGGTCCATATGACCTCGAAGGATGGAATGGAGACGGGTACAAAAATTATCAAAAAAGAACTATCAGTCGACGATGTCGTTCAGGATCCTCGACATGAAATAAAGAAAGAGGAACCTCTTTCTTTACCGTTACCCCGTTTTGCTCTTAAGCACAACACGTCACCTATTGAGTCTCCTCTTCCTTGTCCTTCACCACCTCTATTGCCTGTTGGGCCCAAGGCTCGTTCTCCTTGTTTCTCTCCACAATGCAGTAGCAAATCCCAGTCTCCGGCAATGTCACCTGCTGGTTCGTTATCGGATGTCTCTCATGATACAGTTAATTCTAATCTCCCCGAGTCTGGACTGAACAGGTTACACGAGCAAAAAGATGAAGATCAGCTTAGTAAAGATTTGGCTGCCAGTTCAAGCGGAATCAAAACCCACGTACCAAACCCTGATGGTGACGACAATTTATCAGACAGTTTGAGTGAACTTAGCTGTTCCTCTCTGGATGAAAAAATTCGACAGTTGGATAAAAAACTGAGTATGACTCCTGTTCCACGGCCAGTTGATACAATCACCCCTGGGCTTTATACGAAGTTTAAAATCAAAAAGAAGGAAACACCTCTGTCTGGTGGAAGCATGCTGACATCTAGTCTGCGGAGTGAGCCATCTGATATTGTCAAATCATTGCTTTCTCGGTCGAGTATATTTGATCAAGATTCAAAGCGGTTAGAGCAGATAAATGAAAAATACAAACCAAAAGAGATCAATATAAATATTGAAGGCTCCCCACCCAAAATGAATATCCGTACTAAAGCAGCGGCCAAAGAAATGCCCCCTACAATGCCTCCAAATTTCAGCCCTCTCATGCCCTTTGGTCCTTTTAGTTCCACTTGTAACATTACTTCCCCTTTTCAGCCTTTTACCTCTAACTTGTTCCCGCAGCACCCTGTGATACGCTTGTCAAATACCCCTTCTACTCCACAGTTTTGGGGAGGTTCAACTCCAACCATGAGCTCGCCAACAAAACCACTCGTGGATGTGGCAGCGCCAGTTTCAGTGTTAAAGAAGACTCCGCTCACTCCTGCCACGACAGAGCTGTCGGTCAAAACTGAAGAGTTACCCCATGCTTTGCTTTCTTCGCCACCACCATCTGTGGCAGAAATTAGTCTGCCCCCTGTACGGCAGCCAGTGGATACTCCACTTGCAAGTTCACACACTCCTATGGTGAGTCCATGTTCGGAACAAACTGAACCTGCAGTGAAGAAGGAACCTAGCGTCCCTTTAGTAACTGGTAAAAATGAAGTTGCCGTTTGGATGTCCCCAAGGAAAGATCTCAAAGTAGAAGAGATGAGTGAAGTTGATCCCAAAGAAAAAAGAGAATTTATCCCTGGTCCAGTATCTGTGCCGTCCACATCATCTGGTTCCTGTCTGGGAAAGAGGAAATCCTCAGATGATTTGGATTCCCCCAGGAATAAGATGATGAAAACTGAACCCAAAACTTGTCCTCTTGCAATGAAGCCTCAGGAAGTTCCAAAGAACGACAGCTTTGAGAAGACCGAAGCAAAAGACGAGCCCGTGCAGCAGCCTCGAAAGAAAGCCAAGCACGGCGAGTCTACAAAACCCAAAGAGGAAAAAGAGTGTGAGAAGAAATCTCCTGGGGTACCTTCAGTCAAGCCCAGCAAATCGCCTGAGAAGAAACTGGTAAAGACACCAAGTCGTGAAAAAAACGAAAAAAAGGAGCCAACAGAGGGCGAGAAGAAAGATGTGGAAAAGAAACCAGCTGAACTACCAACTCCGAAGAAAGATCAAAAAGCTTCAAAAGACAGTGATAAAAGCGCTGAGAAATGGAAATCAGGAGGCAAGGAACCAGTCTGTGATTTAGATGAGAAGGCTAAGTTGAAGAATAGCCACTCTGAACGGAAGCATGATAAAAGTGAAAAGAAAGAAAAGAACCGTACAAAAAGCTCAGAGGGGAAGCCTGACAAGGAGAGAAGCTCCTCAAAAGAAGGAGAGGCAAAGGAGTCTCACAAAAAGGACAAATTAAAATCTAAAACTGACGGTCGGTCAAAAACCAAATCGGAAGATAAAGGAGAAGGGGATGCTTCCAAAAATGGGCGGGGTGGCCATGAAAAAATAGAAGCTGAGAATCTATTAAAAGCCTCCAAGTCTCACCGGAAATCTGAAGAAGACAAAGGAGGAGATGACAAAACGGACAATGTAAAATGCTCTGAAGTTAACAAGATTGTTAAATCCCATCATCGTATTGAGCCACCAGAGGGCGAAGCGTCAACAAAAAACAATAAGACACATCACAAAACAGACCCAGAAAAAAGTGAACTTGAAACCTGCAAGGCGGGTAAGTCTCACCATAAGGCTGACAATGATAAAATAGATGGGGAAAGTGTGAAGTCTTCTGGAAAAACTCACCAGAAACCAGCAGAGGGCCTGGAAAAAGTTGAGGCTGATGCCGCCAAGAATAAATGTAACAAAGTGGAAACTGAAAAGACAGAACCAGAAAATACTAAGGCCAGTAAATCTCACCATAAATATGAATGTGAAAAGGGTGAGGCCGAAGGAACAAAACCTAGCAGCAAGTCTCACCACAAACACGACCAATTAAAGTCTGAGGCAGACAGCTTGAAATCGAGTAAGTTGCATCACAAGTCTGACAAAAGTAACCGCAGTTCTTCAGATCGGACGCGTAAAGATAGCGTTAAAAGCAGCTCTGATCGTAAGAAGAACGAAAAAGAGAAAATGGACAAATCTGAAAAACAGGAGAAGTCTGATCGGTTAGAAAAAATGGATAAACCTGACAAAGATACCAAGACCGAAAACAGGGAAAAGTTAGAAAAGACTGAAAAGTCCGAGAAAGAAAAATTTGATAAAGAAAAATCCGATAAAACGGAAAAGGAAAAATCTGACAAGGAAAAGTCCGAGAAAGTGGAAAAGGCGGAAAAAGAAAAGTCCGAAAAGGAAAAATCAGACAAGGAGAAAGATAAAGCTGAAAAGGAGAAGCAAGATAAGAGCAAACCGGATAAAGAAAAACCGGAAAAAGAAAAACCAGATAAGGACAAACAAGAAAAACCCGACAAAGACAAGCAGGAAAAACCTGAGAAAGAGAAACCTGACAAATCTGAAAAGGAAAAACCAGACAAAGAAAAACCCGATCGATCAGAAAAGGATCGGTCAGAAAAGGACAAGTCCGACAGGTCTGAAAAATCTGACAGATCTGAGAAAACTGAGAGCACCGAAAAGATGGACAAGTCCGAGAAAGAAAAGTCCGATAAGGAAAAGTGTGGAGATGAAAAGAAGAGCGAGTTCTGCCTCGATGCCTTTGGTCCGTATGTGTCCATGTACGACAAGGTGAAGAGACGATCCTGCAGTAACAAGGACAAAGACATGGAAGACATGCGCAAGAAACTGAGCCAACTTAAAAACACTCGTAAGAAAAGAGGCAGCAAACCAAGTCGAAGTGATGAGACAGAGAGCAGCGTACAGAACACTGATGATGAGAGTAGTACAAGCAATTCGTTTGTGGCTTCCAAGAAGGAGACGCTAACAGAAGAAAAAGATGAAACTACAAACAAAAAGAAAAAACGAAATGTAATCGAGAGCTCCTCTTCAGAAGAAGACCCTCCAACATTTATACCCACACCCTCCAAGAGAGAACTGAAGAAATCTTCTCCAAAGACCCGACCTATGCACAAAAAGTCCAAGGCCGTCTTAGAAGTCTCTTCTGATAGTGATACTGAAGAGAAAGATTATAAAAAGATCTCATCTAAACTGGCCAAGCAGAAGTTGGGCTCATCATCAACGTTGCTGCATGACCAGAGAAGCAACTCCACCGATGCTAAAAAAGACAGTCAGGATGAAATGATGACAGAAAAGACAGACAAACCGAATCGGAAAAAGTCTAGTAAAGAAAAGAAGAAGAAATGCGAAAAAGAGAAGAAGAAATCCAAGAAGTCTGCTGACAGGTCTTCTAAAGATGAGAAGAAGGTTTTACCAAAATCTTCACCTGTGTTCAATGATATCAGTGAAATGCCAGGCATGGACTCACAAGGTTCCCAGGTTGTTAACATCCTGTCGGACAACATATACATGTCCATAGCCCATGAAGATTCTCCCAAGGAGGATTATCCCCAGATGCAAAAATCAAATAAGCCATCAGTGGCAACAGCATCAACAGTGGAAGACGAAGAAGATGAGGAAGAAAATGCTTCTGATGTTGCTATTGGCACCACTCCTGGCACGAATACCGCAACTGCAGATGAAGACGAGTCCAAACCTTTTGCCCAGCGGAAATCAAGCAAAGATTCCTTCCCTCACAAAATAGAAGAAAGCCATCAGGTTGAGAGCTTGTCAGAGAAACTCTTCCTACCTATAAGTCTGAGTGAGGCTACAGCCAAGGACTCCGAGGAGAACCACCACGCAGATACAGAGGATGAGCCCCCTGACATGTCATCCAAATCGAACAAGAGCAAGAAGAAGAAAAAGAAAGAGAAGAAGAAACATGATACTAAGTCTAAGAAAAGTTCCAAGAAAAACATTGTGAAGATCTCAATGCTCGATGACAATGTTTTCCTAAATGATCAGCCTGAGGTTACGCAAGAGGAAGAGGATGAAGAAGAGGAAGAAGAAGGTACAAAGTTCTTCCAGATTTCTAGCCTTGAGTCTTTTTTTCTTTCTGGCTCAAAAGAAGAACCAAATGAGCTCATAAAAACCGAGCCCCTTGTCGCAGAGGAAGATATGAAAGAAGAAACCGATGTGGGTAGTCTAAAAAGGACGTCATCTTGCGAAGATATCAAGGTAGCAACTAACGATGACTATTTGGAAGAGCCTGTAAAAGAAGAGGCTTCCACAGATTACATGGAGCCTGAAGACCCGCTTCCTACCAAGGTTGAAGAGCCACCAAAGGAAGAACTGGAGGCTGAAGAAGATGAAGTTGAAGAGGAAATCGAGGAGATCAATAAACCAAATACCGACTCCTCCCCGCTGGAAGAAGACAACTCCAACAATTTAGAAAAATGTGAAAGTTTTCAAGATGACGTCTACGATTTTGACAAGTCGGAAAATATCACTGAGAAACTGACTGATTGTGAACAGTTGAAGAAATTTGAAAAAATCTACGACTCTGAGAAGGTGAAAAAGGTTTCTGATGAAGCTGTCGATTCTGACATTCCAAAAAGTGTTGGAGAAGAGTCCCCTGAGTTCCTCTTATCACCAGGAGAGGAGGTCACTGAGGAAGACCCTGATCAACCAGAACCCACACCTGAAGAGCAGTCAGATAACGATGTTCAAAAAGAAACAACAGACGATATTCAAGCTCAAAGCCCTGAGCACTCAGATACAAAAGAGTTGGACAGTCCTTATCAGTTAAAAGAACAAGAAGATACTGTATCGGATGTGTTTACTGATTTGAAGAAACTCGATAAAAGCTGCTCTGAAGTTGAAGACAACAGTGAAGATAAGATTGTTGCAAAAAAGCCTGAAAATGAGGATTTACTTGAAAAAGAAGAACCAGAAAAGACTGAAGAAGAATCGAAATCTCCTGAAACCTCAGAAGTTGAAAAAGATCCGACAGACGAATCTTGTGAGACAACTGTAGATGGTGAGGATCTGAAACGAGATGAAGGCAAAGATCAAAATTCCAGCAGTTTTGAGAAACTAGCCGAAACCAACGAAAAGGTAATTGGTATGTTCAACTTTGGCCAATTACAAGTTCATGTTGTTGACAAGAATACAGGTAACCTTCACAAGGTGCATGACGCCGAGTTATCCGATAAAGATAGTCGGGCCGCTCTTGGTGTGGATCAGGTTCTTGACAAAGCTAACAGCAACACTGATGCCAACAGGAATGTAGACCACATCTATGACTTTGACGATACCACTGAAGAAATGTACGGTAATACTTATCTCGATAAAGGAGATAGCAAATCCAGACAGTTCAACTACAAAAAGGTCTTTATGATGCAACCCCGTAAATCGGACAGAATTTATAATTCTGATGTCACCAAGAAAATCTGCGAGAAGATTGAGCAGAAGATTTACAGCCTAGAAACGCCAACCAAATTTGATACCATCTTCAAGAAAATTAACAAAAAGGCGTACAATATGGACGAGGCTAATTCTACTGAAGATAAGATTTATGAGGATGAGATTTGCAAGACGCGCAATCGAATTGCAGATCTCCACCACCCGCAGCAGCATACAGATAAAACAGATGGTGACGTCTATGATTTCGATGAACAAGACAAGCTTGAAGTTAATGAATTGTTAGATTCCAAGAAAACGGATCATATAATCCATGATCTTCACGATTACAAGAAAACCGAAGAGCATTTACCCAATTTTGATGGAGGCAAGAAAAGTACCGAAAAACTGGCAAGGTTCCAAGATTGGTTTACAAATGACGATAAGATATATTACTTTGACCCTGTTAAAGATTTTGGTGCCAAGCCTAACAAAAGCCGTTCCAAATGTGGTAAACACAAATCTGGTCAGGGGAAACGCCGCAATAAAGGTAAAATTGTTACCTTTACTTGGCCTCGTGAACTAGACTTAAAATCATCGGTTGAATCTGCAACGGACAATGATACTGTCACTGCCGATACGTTTGAGATGAAAGAGAGCGATTATGCTGAACTCGATGACCAAGAAGATAAACTTGTCATCGATACAAGTGACATGTGCGTCTATGATGAAAATTCAACCGAACCGATTGAAAAAGAGGCTGATAAAGCAAAGCTGTGTAGCTCGGAAACAAAAACCGAAAACGATTTAACAGCCAATGATGAAAAAAACGAAGATATCCAAGACGACGTGAAAGAAAAGGTGTCAAGTGAAGAGGAAGAAAAAGAGCAAGAAGAAACGGAGGCTCTGGATACCACGTCTTTTAAAAATAATCAAGACACTTTTGGAATGTCGATCAAACAGAAAGTGGTGGTGAATTATGGACTAACAAATCACAGTGGCACAACCCCCGACAAGATAGAATCGGTTATAGGCTTTGGAATGCTATCTGGCAACAAAGACCTGCTGTCGGCATCTAAATTGTCTGATGACTGTGACAGCAGCGATGAAGGTGCGCTCAGAATCGATGAAGATGAAGTAGACGAAAATACTGATAGGGATGAAATGACACTTGAAGACAAAACCCCAGAAAAAGAGAAAGAAATTGAAGAAGAAGAAAAAGTGATCAAGAAGCCTGAACCCTATGAAGAAGAGGAGGAGGAAGAAGC'
#sequence = 'AAACGCTAGAAACCATGGAAACCCACTATC'
missed_cleavages = 1
min_length = None
a_to_g_editings = [0,1,7,18,19,20,30,34,38,45,57,62]
#c_to_t_editings = [3,23,25,41,51,53]
#a_to_g_editings = [0,1,7,18,19,20,30]
#c_to_t_editings = [3,23,25]
#a_to_g_editings = a_list
#c_to_t_editings = c_list

look_ahead_for_editnig = 3
look_behind_for_editing = 6


#run over all peptides in initial_peps
def find_all_edited_peptides_from_seq(sequence,a_to_g_editings,c_to_t_editings,digestion_rule,missed_cleavages=0,min_length=None,look_ahead_for_editnig=3,look_behind_for_editing=6):
    
    ml = missed_cleavages+2 #size of normal window (2 boundaries + missed cleavage sites)
    trange = range(ml) 
    cleavage_sites = deque([0], maxlen=ml) #the normal, dynamic window for digestion
#    all_a_to_g_sites_it = it.chain(a_to_g_editings)
#    all_c_to_t_sites_it = it.chain(c_to_t_editings)
#    a_to_g_sites = deque([])
#    c_to_t_sites = deque([])
#    a_to_g_current_site = None
#    c_to_t_current_site = None
    original_cs_list = find_cleavage_sites(digestion_rule,sequence)
    initial_peps = cleave_rna_as_peptide(sequence, original_cs_list, missed_cleavages=missed_cleavages, min_length=min_length)
    peptides_list = []

    for pep in initial_peps:
        
        find_all_perm_for_pep(pep,sequence,digestion_rule,missed_cleavages,original_cs_list,
                              a_to_g_editings,c_to_t_editings,peptides_list,
                              look_ahead_for_editnig=look_ahead_for_editnig, look_behind_for_editing=look_behind_for_editing)
                              
    return peptides_list


def find_all_perm_for_pep(original_pep ,sequence, digestion_rule, missed_cleavages, original_cs_list,
                          a_to_g_editings, c_to_t_editings, peptides_list,
                          look_ahead_for_editnig = 3, look_behind_for_editing = 6,
                          cleave_before_pattern=False, temp_edited_peps = []):    
    
    #determine the extanded window in which try permutations
    editing_window = determine_rna_pep_window_range(original_pep,look_behind_for_editing,look_ahead_for_editnig)
    print('\n========================================================================================================')
    print('peptide: ' +  original_pep.seq + ' | peptide coor: ' + str(original_pep.coo) + ' | cleavage sites: ' + str(original_pep.cleavage_sites) + ' | editing window: ' + str(editing_window))
    print('\n')
    
    append_zero = False
    append_seq_end = False
    first_original_pep_cleavage_site = original_pep.cleavage_sites[0]
    last_original_pep_cleavage_site = original_pep.cleavage_sites[-1]
    
    #determining if seq beginning and seq end are in cleavage sites
    if 0 in original_pep.cleavage_sites:
        append_zero = True
    if len(sequence) in original_pep.cleavage_sites:
        append_seq_end = True
    
    #find all eiditing permutations, a permutation is ([some a_to_g sites permutation],[some c_to_t sites permutation])
    editing_permutarions = find_permutations_in_range(creat_editing_sites_permutations_one_type(find_editing_sites_in_window(editing_window,a_to_g_editings)),
                                                      creat_editing_sites_permutations_one_type(find_editing_sites_in_window(editing_window,c_to_t_editings)))

    #for each pemutation produce new sub_sequence (with window range) and check for new cleavage sites relateive to whole sequence
    for perm in editing_permutarions:
        
        #the edited sequence (extended to look_behind and look_ahead codons)
        new_pep = rna_peptide(edit_rna_as_peptide(sequence[editing_window[0]:editing_window[1]+1],editing_window, perm[0], perm[1]),editing_window, a2g = perm[0], c2t = perm[1])
        new_pep.find_cleavage_sites(digestion_rule,cleave_before_pattern=cleave_before_pattern,append_zero = append_zero, append_seq_end = append_seq_end)
        print(new_pep.seq + ' | coordinates: ' + str(new_pep.coo) + ' | permutation: ' + str(new_pep.a2g) + ' ; ' + str(new_pep.c2t) + ' | new cleavage sites:' + str(new_pep.cleavage_sites))
        
        edited_cleavage_sites_in_original_pep = [c for c in new_pep.cleavage_sites if (c>=original_pep.cleavage_sites[0] and c<=original_pep.cleavage_sites[-1])]
        edited_cs_in_original_pep_length = len(edited_cleavage_sites_in_original_pep)
        
        
        #boundaries are the same
        if first_original_pep_cleavage_site in new_pep.cleavage_sites and last_original_pep_cleavage_site in new_pep.cleavage_sites:         
            
#            edited_cleavage_sites_in_original_pep = new_pep.cleavage_sites[new_pep.cleavage_sites.index(original_pep.cleavage_sites[0],0):new_pep.cleavage_sites.index(original_pep.cleavage_sites[-1],1)+1]
            
            #number of cleavage_sites in new pep are the same and therfore it is a valid peptide
            if len(original_pep.cleavage_sites) == len(edited_cleavage_sites_in_original_pep):
                peptide_permutation = rna_peptide(new_pep.seq[original_pep.coo[0]-new_pep.coo[0]:len(new_pep.seq)-(new_pep.coo[1]-original_pep.coo[1])],original_pep.coo,a2g = perm[0],c2t = perm[1],cleavage_sites = edited_cleavage_sites_in_original_pep)
                peptide_permutation_key = (peptide_permutation.seq,peptide_permutation.coo,peptide_permutation.a2g,peptide_permutation.c2t,peptide_permutation.cleavage_sites)
            
                #peptide wasnt already appended for some other permutation which varies in an out of range (peptide coordinate) site
                if peptide_permutation_key not in temp_edited_peps: 
                    peptides_list.append(peptide_permutation)
                    temp_edited_peps.append(peptide_permutation_key)
                    print(str(peptide_permutation_key) + ' appended - regular')
                    peptide_permutation.print_peptide_data()
                    print('\n')
                else:
                    print('already exist\n')
              
            #number of cleavage sites in new pep is greater than in original pep --> new cleavage sites due to editing
            if len(original_pep.cleavage_sites) < edited_cs_in_original_pep_length:
                #creating all new peptides in the new pep range
                new_peptides_from_new_cleavage_sites(original_pep,new_pep, perm, peptides_list, edited_cleavage_sites_in_original_pep, first_original_pep_cleavage_site, last_original_pep_cleavage_site,
                                                     look_ahead_for_editnig = look_ahead_for_editnig, look_behind_for_editing = look_behind_for_editing, temp_edited_peps = temp_edited_peps)
                
            #from new cleavage sites to other downsteam cleavage sites, create all peptides.
            for i in range(edited_cs_in_original_pep_length):
                if edited_cleavage_sites_in_original_pep[i] not in original_pep.cleavage_sites and edited_cs_in_original_pep_length-1-i<=missed_cleavages:
                    print('!!!')
                    #get index of cs in original cs list from which to create a new 
                    index = bisect_left(original_cs_list,edited_cleavage_sites_in_original_pep[i])
                    print(str(index))
                    cleavage_sites_in_new_window = list(set(edited_cleavage_sites_in_original_pep[i:] + original_cs_list[index:index+missed_cleavages+1]))
                    cleavage_sites_in_new_window.sort()
                    print(cleavage_sites_in_new_window)
                    print(sequence[cleavage_sites_in_new_window[0]-1:cleavage_sites_in_new_window[-1]-1])
#                    new_peptides_new_ranges = cleave_rna_as_peptide(sequence[cleavage_sites_in_new_window[0]-1:cleavage_sites_in_new_window[-1]-1], original_cs_list, missed_cleavages=missed_cleavages, min_length=min_length)
                    
                
                
                
def new_peptides_from_new_cleavage_sites(original_pep, new_pep, perm, peptides_list, edited_cleavage_sites_in_original_pep,
                                         first_original_pep_cleavage_site, last_original_pep_cleavage_site,
                                         look_ahead_for_editnig = 3, look_behind_for_editing = 6, temp_edited_peps = []):
    
    cleavage_sites_for_new_peps_in_range = [l-new_pep.coo[0] for l in edited_cleavage_sites_in_original_pep if (l>first_original_pep_cleavage_site and l<last_original_pep_cleavage_site)]
    new_peps_in_range = cleave_rna_as_peptide(new_pep.seq,cleavage_sites_for_new_peps_in_range, start_coo = new_pep.coo[0], missed_cleavages=missed_cleavages, min_length=min_length)
    print('cleavage_sites sent: '+ str(cleavage_sites_for_new_peps_in_range))
            
    print('new sub peptide due to new cleavage site')
    for x in new_peps_in_range:
        
        #the extanded window of x pep (with lookahead\lookbehing)
        extended_x_range = [max(x.coo[0]-look_behind_for_editing,new_pep.coo[0]),min(x.coo[1]+look_ahead_for_editnig,new_pep.coo[1]-1)]
        a2g_sites = find_editing_sites_in_window(extended_x_range,perm[0])
        c2t_sites = find_editing_sites_in_window(extended_x_range,perm[1])
                    
        #the new final rna peptide coo
        final_rna_coo = [max(original_pep.coo[0],x.coo[0]),min(x.coo[1],original_pep.coo[1])]
#        final_rna_seq_range = [extended_x_range[0]-new_pep.coo[0],extended_x_range[1]-new_pep.coo[0]+1]
                
        #the new final rna peptide cleavage sites
        first_cs = []
        last_cs = []        
        if x.cleavage_sites[0] == 0:
            first_cs = [first_original_pep_cleavage_site]
        if x.cleavage_sites[-1] == len(new_pep.seq):
            last_cs = [last_original_pep_cleavage_site]
        final_rna_peptide_cleavage_sites = first_cs + [k for k in x.cleavage_sites if k not in [0,len(new_pep.seq)]] + last_cs
                    
        #the new real rna peptide to add
        peptide_permutation = rna_peptide(new_pep.seq[final_rna_coo[0]:final_rna_coo[1]+1],final_rna_coo ,a2g = a2g_sites, c2t = c2t_sites, cleavage_sites = final_rna_peptide_cleavage_sites)
        peptide_permutation_key = (peptide_permutation.seq,peptide_permutation.coo,peptide_permutation.a2g,peptide_permutation.c2t,peptide_permutation.cleavage_sites)
        print(x.seq + ' | coordinates: ' + str(x.coo) + ' | permutation: ' + str(a2g_sites) + ' ; ' + str(c2t_sites) + ' | new cleavage sites:' + str(x.cleavage_sites))
        peptide_permutation.print_peptide_data()
            
        #peptide wasnt already appended for some other permutation which varies in an out of range (peptide coordinate) site
        if peptide_permutation_key not in temp_edited_peps: 
            peptides_list.append(peptide_permutation)
            temp_edited_peps.append(peptide_permutation_key)
            print(str(peptide_permutation_key) + ' appended - new cleavage site - in original range')
            peptide_permutation.print_peptide_data()
            print('\n')
        else:
            print('already exist\n')
    print('new sub peps due to new cleavage sites ended\n')
                

def find_cleavage_sites(digestion_rule,sequnce,cleave_before_pattern=False,append_zero=False, append_seq_end = False):
    
    if append_zero:
        initial_list = [0] + [x.end() for x in digestion_rule.finditer(sequnce)]
    else:
        initial_list = [x.end() for x in digestion_rule.finditer(sequnce)]
    initial_set = set(initial_list)
    
    #find other overlapping sites
    if cleave_before_pattern:
        for site in initial_list:
            match1 = digestion_rule.match(sequnce,site+1)
            match2 = digestion_rule.match(sequnce,site+2)
            if match1:
                initial_set.add(match1.start())
            if match2:
                initial_set.add(match2.start())
                
    else:
        for site in initial_list:
            match1 = digestion_rule.match(sequnce,site-1)
            match2 = digestion_rule.match(sequnce,site-2)
            if match1:
                initial_set.add(match1.end())
            if match2:
                initial_set.add(match2.end()) 

        
    #return only sites representing in-frame codons
    cleavage_sites = [x for x in initial_set if not x%3 and x!=len(sequence)]
    
    if append_seq_end:
        cleavage_sites.append(len(sequnce))
    
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



#a = find_editing_sites_in_window((18,56),[1,2,4,5,16,17,29,30,31,54])
def find_editing_sites_in_window(editing_window,editing_sites):
    editing_sites_in_window = []
    editing_sites_it = it.chain(editing_sites)
    for i in editing_sites_it:
        if i<=editing_window[1]:
            if i >= editing_window[0]:
                editing_sites_in_window.append(i)
        else:
            break
    return editing_sites_in_window
    
        

def update_editing_sites_in_window(editing_window,editing_sites_list,editing_sites_it,site=None):
    
    #cleane left side sites if out of range
    if len(editing_sites_list):
        k = editing_sites_list[0]
        while k < editing_window[0]:
            editing_sites_list.popleft()
            k = editing_sites_list[0]
    
    #append current site if in range
    if site != None and site > editing_sites_list[-1]:
        if site <= editing_window[1]:
            editing_sites_list.append(site)
    
    #keep appending sites in range until reaching an out of range site
    try:
        for next_site in editing_sites_it:
            if next_site <= editing_window[1]:
                editing_sites_list.append(next_site)
            else:
                break
        #return current site reached by the sites iterator and all a list of all editing sites in range
        return next_site, editing_sites_list
    
    except:
#        print('iterator exhusted')
        return None, editing_sites_list
    


        
#peps = cleave_rna_as_peptide(sequence,initial_cleavage_sites_list,1)
def cleave_rna_as_peptide(sequence,cleavage_sites_list, start_coo=0, missed_cleavages=0, min_length=None):
    
    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
        
    for i in it.chain(cleavage_sites_list,[len(sequence)]):
        cleavage_sites.append(i)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append(rna_peptide(seq,[start_coo+cleavage_sites[j],start_coo+cleavage_sites[j]+len(seq)-1],cleavage_sites = [k+start_coo for k in list(cleavage_sites)[j:]]))
                        
    return peptides

def determine_rna_pep_window_range(pep,look_behind_for_editing,look_ahead_for_editnig):
    
    #peptide bounderies
    if pep.cleavage_sites[-1] == len(sequence)-1:
        downstream_bound = len(sequence)-1
    else:
        downstream_bound = pep.cleavage_sites[-1] + look_ahead_for_editnig - 1
    if pep.cleavage_sites[0] - look_behind_for_editing < 0:
        up_stream_bound = 0
    else:
        up_stream_bound = pep.cleavage_sites[0] - look_behind_for_editing
        
    return (up_stream_bound,downstream_bound)



"""
return 2 lists of 'A' indeces and 'C' indeces' in sequence
"""
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


peptides_list = find_all_edited_peptides_from_seq(sequence,a_to_g_editings,c_to_t_editings,digestion_rule,missed_cleavages=missed_cleavages,min_length=min_length,look_ahead_for_editnig=look_ahead_for_editnig,look_behind_for_editing=look_behind_for_editing)

print('\n\n peptides so far:')
for pep in peptides_list:
    pep.print_peptide_data()

