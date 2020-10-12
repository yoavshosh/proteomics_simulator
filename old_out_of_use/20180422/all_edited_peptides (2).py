import re
import itertools as it
import sys
from rna_peptide import rna_peptide
from collections import deque



sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/' +'test.txt', 'w')

digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
#sequence = 'AAACGCTAGAAACCATGGAAACCCACTATCAGGCATTCAGGCATCAGCATGCGCCCTATTGCACGC'
sequence = 'ACTACCACCACTACTACTACTACACAAATTTTCTATGGCCATGTAATCATCATGTAGTAAGTAGTGAGAAGAGATGGTGCGGCAAACGAGGTATCTCCTGGTTGGGAATATACCCGAAAAAATCACCAAGGAGAAAATCACAGAACACTTTAAGAACTATGGAAAGATCCAGAATGTAAAGTTTCATCCAAAAAAAGAAAATGAACTTGGCATTACAGCAACTGTGGCATTTATGGACATTAAATGTGCCTCTAAGGCATATAGTTCAGAGAACAAAATAGAAAATGTAGTGCTTCGGACAGAATACAGTGAGGGCACTGCAACAGGAAGTGTTGTTACAAGACCTGTGGACCCTGGATTGCCTACAACTGGAGCCCACAGGGGACAGTACCTACATGCTCGAGGGCCATCTACTAGCTTCACAACACGGACTAAAGGTAGTCTTGATGAAGATGATACATACAGGGAATATTACTATGGTGGTAGCGGCAGCACTGGCAGTAGGGAAAACAGTCAATTTGAAAACAGATATCCCACATATGATGAACAAGTACAGCATCAGTCAAGGTCTTATCAAAGAAATAGGACCTCATTTTCCAAGCAGTATGAGCGAGGACCACAACGACAGAAGTCTAGCAGTGCGAGTGTCACGGGTGTTGTCACTAGCGGTGGTAGTGGCAGTAGTGTCACAGGGGTTGGGAATAGCAATACAGTTGGCAGTACAGTGACACCCAGCAGTGGTGGTAGTGTGGTTGGGAGTGGAGCTGTTGGAACCAGTGGCAGTGGGGGTGGGGGTGCTGGTAGTGGTGCTGGCAGTGGAGGAGGCAACAATAGTGGAGGTGGTGGCACCAACAATGGAGGGGACACTGGCAACACTGGGGTCAGTGTTGTGGTCTCTAATACCTTGCAAAAGCATTCCTTTGATGAGCGTTACCACGAACAGTATAGCAGCGACAGAGAAAGTTTTGAAAGGCCACGGATATCCACATTGTCTCCAACATGTTCCATCTCTCAGGCCAGCAGCCGACACAGTAGTGGTCGGCGGAAAGGTTCTGTCTCTGTGTCCCGTTCGTCCAGTCGCAGTCGGTCAAGGTCCAGTTCCAGCAGTGAGTATAGTCGGTCTTCATCAGGAAGCAGTGGATCAAGATCACATTCCAGAAGTCGTGCAAGCAGCCGGGGGAGTTTGCGGCGAATACCTGTTGGTTCAGCATCGAGTCAAAAACACTCTAAATCAGAAAGTCAAAGTATGATGTCCAGTAAATCCTCAACAGCAGGATCGCAATCGGAAAGAAAACCACTTGGGATTTGTGTAACGGGTTTACCATTAAGATCCAGTGATACTAGTCTTCGAGATGGTCTTTTTCATGAATACAAAAAACATGGAAAAGTTACTTCTGTACAGGTGCTAGGACAAGGCGAAGAAAGATATTCAGTGGTCAGTTTCCGGAAACCTGAAGATGCGGCAAAGGCATTAGAAGCTTCACAGGACAAAATGTTTTTTGGAAAGAAAATCAAAGTTATGGCTCATGAAGGAGTTGAAGTGGAAGACAATGAGTTCCGTCCCCCAGAAGCTGAATGTGATGAATACCACCCTAAGGCAACACGAACGTTATTTGTGGGAAACTTAGAAAAAGAAATAACGACACAAGAGCTTAGAGATCGTTTTAAGTCTTTTGGAGAAATTATAGATATTGATATCAAACGTCAGGGTGCAGTGTCAGCATATGCATTTGTACAGTATGCAGATATTACAAGTGTTGTAAAATCACTCCGAAAAATGGAGGGCGAACATATAGGAGCCAATAAAATAAAGCTTGGTTTTGGTAAAAGTATGCCAACCAATTGTGTGTGGCTCGATAACCTGTCTGAATCCGTAACAGAGAAGTTTCTTTGCCAACAATTTGGCCGTTATGGCGAGGTCAGTCATGCAGTGATTGATCGGTATAAAGGAAGGGCCTTGATATATTTCACAAGTATGGATACAGCTCAATATGCAGTCACTGAGATGAGAAATAGAATTTTGAAGAAGAAACGAGTGCAGACTGGTGACTTTAGACTCGGTGGACATCTGAATATTGATTTTGCTAGTCGAGATTGTCAAACAGCTTTCTTTGAAAAAATGGAGCGGACGGGTCAACTGCAACCTGGAGAACGTCCTGATGAACGAGGAGGGCGCCTCTACCGGCACAATAGGGGAGCCAATTTCGAGCAATATTATGATTCTCCTGCAACCCCTTCCGCTTACAAGGAAGATAGCGCAAGATACGAGACTGGCACTGTAACACCTGTTGCTACAACTGCACGAGTTGGACTGACCCCACCACAGTTCAGTACGACAAAACGTAATCGAGCAACAAATTTTCAAGCCTCTTCGACCAGACAAGGATCCTCATACAGAGGAAGACATTTTACAGAAAGTACTTACTCAGAAGAATTCCCGGCAAACCAGAGACAGAGACATATTGATGAATACAGTCAAGGAAGTGCTCCTCCTTATGCTGAGGAAGATTCCTATGAACAGGAGCTGAGAGAATATGGCTACAGTCAACGGGAACGCCGTGAGAGGGGTAACGTAACTCCACCGAGACGATACCACAGTCCTTATCGGGAAGATAGGCATTCAAATGCCGACAAAGAATCGACGCGGGGTAGCCGTTATGAAGAAGCCTTCACCAAAGACAGATATGATTATTACAGTAGTGACAGTCCTGTCGTCAGTGACAATGGTCATGATGAGAACTGGGTGGATCATCATACATTGGATCTTTCAGATCGTCATGCCACACCAACAAAAATGAAATATTCAAGTCGGAACGAGACACGCTTAAAGCAATCGTCTCGGGATCACAGCCCTGCCATTGAGTTTTCCTACAGGCACTCACCAGTGAGACGAAGAAGCAAGAGCAGGAGCAGAAGTAAAAGTGGGAGTCGGACAGAAAGTCCTAGTCGAAGTAGGAGCCGGACGCCTCTAAGGGAAAGGCGGAAAGAATCTGTGTTCCGGCCACGTTCTCCCCACACACCACAAACACCTCCTACACCAACAGAAGAAGAACCTCCCATTTATTCTAAAATCCGATCAATCAGCCCCACACCCCGGTTTGAGGAATTTGAGCGAACAGGAGAATTGAAAAAAATATGCAGGACAGATGATAAAATTGAAGGCATTAGGTCTTCGGATCATGAAAAGGTGTACAGTATAAAATGCTTGGGTACTCCATCTATGTCTTTCATGAAAGGGGAAACGACACAAGCCCGTAAAAAGAGTAGGGAGTTTGATGATAGCAAACATGATGTAATAGCAGGACGTCTGCATAAGTCCCAGGTGCCTGCAACGAACCGGCTGTTAGAATCTGCGCTCCAAGACAGACATCAAGTCTTGTTAAAGAAGGCAGATGAAAACCGTCTCAATAGACGCAAACTATTAGACCAGCCACATGGTAGCTCGAGTCATTCCAGTGGCCACGATGATTCAGATCCTGGTTCTTTTGCTGAGACTGATCTTGGCCATTTGCATCGAGAAAAGCGACTACTCCTAGAAAAATTAAAACAACTGGAAGACGCTGGAAGCCCCAGTGATAATGATTCTTTACCGACACATGATGATCTGCGAGAAGATGGACACGGGCGATTTGTACCGAGGAAAGGTCGACCTGAGGATCCTCATTGTAGCAGTTTACGGCATAGTCTGAGCGACATACACGCAAAGGAAAAACGAAAACTAGAAACTTCTCAGAGTTTTCGAAAACAGATGGAAGCCAGGCGCCTGTTAGAACAGACTCAGGGAGAGTCAAATTTGTTATCAAAACCAAAAAAGAGTGTTGAAAAAGTTGTGATGCTTGATGGTATGGATAGCGAAGGAGAGGAGGAAATAAGTATTACCTCTCCGAAAGCAATCCCCACTGCATTCAAGCCAAAGAGAAAGAAGAAACTCGAAGAAGACCCAGTGTCTGGTGTTCGTACTGGCAGGTTTTATCGCATGAAGAAAAGTGGAATGATGAGCAGTGATGACGAAAAAGATCATTCAAAATCTCTCGAAGCAGAAGAACCTGTTATGCGACGCACCAGTCGCGAGAGCAAAGACGAACTATCCGCTGTTATAAGCCGGTTTAAGGAGCAGTCAAGACATTCTGAATTAAGCCGGGAGGTTTCCGAGTCTGGGTCTCACAAGGAATCCAGAGAGTATAAAATGGTCCATATGACCTCGAAGGATGGAATGGAGACGGGTACAAAAATTATCAAAAAAGAACTATCAGTCGACGATGTCGTTCAGGATCCTCGACATGAAATAAAGAAAGAGGAACCTCTTTCTTTACCGTTACCCCGTTTTGCTCTTAAGCACAACACGTCACCTATTGAGTCTCCTCTTCCTTGTCCTTCACCACCTCTATTGCCTGTTGGGCCCAAGGCTCGTTCTCCTTGTTTCTCTCCACAATGCAGTAGCAAATCCCAGTCTCCGGCAATGTCACCTGCTGGTTCGTTATCGGATGTCTCTCATGATACAGTTAATTCTAATCTCCCCGAGTCTGGACTGAACAGGTTACACGAGCAAAAAGATGAAGATCAGCTTAGTAAAGATTTGGCTGCCAGTTCAAGCGGAATCAAAACCCACGTACCAAACCCTGATGGTGACGACAATTTATCAGACAGTTTGAGTGAACTTAGCTGTTCCTCTCTGGATGAAAAAATTCGACAGTTGGATAAAAAACTGAGTATGACTCCTGTTCCACGGCCAGTTGATACAATCACCCCTGGGCTTTATACGAAGTTTAAAATCAAAAAGAAGGAAACACCTCTGTCTGGTGGAAGCATGCTGACATCTAGTCTGCGGAGTGAGCCATCTGATATTGTCAAATCATTGCTTTCTCGGTCGAGTATATTTGATCAAGATTCAAAGCGGTTAGAGCAGATAAATGAAAAATACAAACCAAAAGAGATCAATATAAATATTGAAGGCTCCCCACCCAAAATGAATATCCGTACTAAAGCAGCGGCCAAAGAAATGCCCCCTACAATGCCTCCAAATTTCAGCCCTCTCATGCCCTTTGGTCCTTTTAGTTCCACTTGTAACATTACTTCCCCTTTTCAGCCTTTTACCTCTAACTTGTTCCCGCAGCACCCTGTGATACGCTTGTCAAATACCCCTTCTACTCCACAGTTTTGGGGAGGTTCAACTCCAACCATGAGCTCGCCAACAAAACCACTCGTGGATGTGGCAGCGCCAGTTTCAGTGTTAAAGAAGACTCCGCTCACTCCTGCCACGACAGAGCTGTCGGTCAAAACTGAAGAGTTACCCCATGCTTTGCTTTCTTCGCCACCACCATCTGTGGCAGAAATTAGTCTGCCCCCTGTACGGCAGCCAGTGGATACTCCACTTGCAAGTTCACACACTCCTATGGTGAGTCCATGTTCGGAACAAACTGAACCTGCAGTGAAGAAGGAACCTAGCGTCCCTTTAGTAACTGGTAAAAATGAAGTTGCCGTTTGGATGTCCCCAAGGAAAGATCTCAAAGTAGAAGAGATGAGTGAAGTTGATCCCAAAGAAAAAAGAGAATTTATCCCTGGTCCAGTATCTGTGCCGTCCACATCATCTGGTTCCTGTCTGGGAAAGAGGAAATCCTCAGATGATTTGGATTCCCCCAGGAATAAGATGATGAAAACTGAACCCAAAACTTGTCCTCTTGCAATGAAGCCTCAGGAAGTTCCAAAGAACGACAGCTTTGAGAAGACCGAAGCAAAAGACGAGCCCGTGCAGCAGCCTCGAAAGAAAGCCAAGCACGGCGAGTCTACAAAACCCAAAGAGGAAAAAGAGTGTGAGAAGAAATCTCCTGGGGTACCTTCAGTCAAGCCCAGCAAATCGCCTGAGAAGAAACTGGTAAAGACACCAAGTCGTGAAAAAAACGAAAAAAAGGAGCCAACAGAGGGCGAGAAGAAAGATGTGGAAAAGAAACCAGCTGAACTACCAACTCCGAAGAAAGATCAAAAAGCTTCAAAAGACAGTGATAAAAGCGCTGAGAAATGGAAATCAGGAGGCAAGGAACCAGTCTGTGATTTAGATGAGAAGGCTAAGTTGAAGAATAGCCACTCTGAACGGAAGCATGATAAAAGTGAAAAGAAAGAAAAGAACCGTACAAAAAGCTCAGAGGGGAAGCCTGACAAGGAGAGAAGCTCCTCAAAAGAAGGAGAGGCAAAGGAGTCTCACAAAAAGGACAAATTAAAATCTAAAACTGACGGTCGGTCAAAAACCAAATCGGAAGATAAAGGAGAAGGGGATGCTTCCAAAAATGGGCGGGGTGGCCATGAAAAAATAGAAGCTGAGAATCTATTAAAAGCCTCCAAGTCTCACCGGAAATCTGAAGAAGACAAAGGAGGAGATGACAAAACGGACAATGTAAAATGCTCTGAAGTTAACAAGATTGTTAAATCCCATCATCGTATTGAGCCACCAGAGGGCGAAGCGTCAACAAAAAACAATAAGACACATCACAAAACAGACCCAGAAAAAAGTGAACTTGAAACCTGCAAGGCGGGTAAGTCTCACCATAAGGCTGACAATGATAAAATAGATGGGGAAAGTGTGAAGTCTTCTGGAAAAACTCACCAGAAACCAGCAGAGGGCCTGGAAAAAGTTGAGGCTGATGCCGCCAAGAATAAATGTAACAAAGTGGAAACTGAAAAGACAGAACCAGAAAATACTAAGGCCAGTAAATCTCACCATAAATATGAATGTGAAAAGGGTGAGGCCGAAGGAACAAAACCTAGCAGCAAGTCTCACCACAAACACGACCAATTAAAGTCTGAGGCAGACAGCTTGAAATCGAGTAAGTTGCATCACAAGTCTGACAAAAGTAACCGCAGTTCTTCAGATCGGACGCGTAAAGATAGCGTTAAAAGCAGCTCTGATCGTAAGAAGAACGAAAAAGAGAAAATGGACAAATCTGAAAAACAGGAGAAGTCTGATCGGTTAGAAAAAATGGATAAACCTGACAAAGATACCAAGACCGAAAACAGGGAAAAGTTAGAAAAGACTGAAAAGTCCGAGAAAGAAAAATTTGATAAAGAAAAATCCGATAAAACGGAAAAGGAAAAATCTGACAAGGAAAAGTCCGAGAAAGTGGAAAAGGCGGAAAAAGAAAAGTCCGAAAAGGAAAAATCAGACAAGGAGAAAGATAAAGCTGAAAAGGAGAAGCAAGATAAGAGCAAACCGGATAAAGAAAAACCGGAAAAAGAAAAACCAGATAAGGACAAACAAGAAAAACCCGACAAAGACAAGCAGGAAAAACCTGAGAAAGAGAAACCTGACAAATCTGAAAAGGAAAAACCAGACAAAGAAAAACCCGATCGATCAGAAAAGGATCGGTCAGAAAAGGACAAGTCCGACAGGTCTGAAAAATCTGACAGATCTGAGAAAACTGAGAGCACCGAAAAGATGGACAAGTCCGAGAAAGAAAAGTCCGATAAGGAAAAGTGTGGAGATGAAAAGAAGAGCGAGTTCTGCCTCGATGCCTTTGGTCCGTATGTGTCCATGTACGACAAGGTGAAGAGACGATCCTGCAGTAACAAGGACAAAGACATGGAAGACATGCGCAAGAAACTGAGCCAACTTAAAAACACTCGTAAGAAAAGAGGCAGCAAACCAAGTCGAAGTGATGAGACAGAGAGCAGCGTACAGAACACTGATGATGAGAGTAGTACAAGCAATTCGTTTGTGGCTTCCAAGAAGGAGACGCTAACAGAAGAAAAAGATGAAACTACAAACAAAAAGAAAAAACGAAATGTAATCGAGAGCTCCTCTTCAGAAGAAGACCCTCCAACATTTATACCCACACCCTCCAAGAGAGAACTGAAGAAATCTTCTCCAAAGACCCGACCTATGCACAAAAAGTCCAAGGCCGTCTTAGAAGTCTCTTCTGATAGTGATACTGAAGAGAAAGATTATAAAAAGATCTCATCTAAACTGGCCAAGCAGAAGTTGGGCTCATCATCAACGTTGCTGCATGACCAGAGAAGCAACTCCACCGATGCTAAAAAAGACAGTCAGGATGAAATGATGACAGAAAAGACAGACAAACCGAATCGGAAAAAGTCTAGTAAAGAAAAGAAGAAGAAATGCGAAAAAGAGAAGAAGAAATCCAAGAAGTCTGCTGACAGGTCTTCTAAAGATGAGAAGAAGGTTTTACCAAAATCTTCACCTGTGTTCAATGATATCAGTGAAATGCCAGGCATGGACTCACAAGGTTCCCAGGTTGTTAACATCCTGTCGGACAACATATACATGTCCATAGCCCATGAAGATTCTCCCAAGGAGGATTATCCCCAGATGCAAAAATCAAATAAGCCATCAGTGGCAACAGCATCAACAGTGGAAGACGAAGAAGATGAGGAAGAAAATGCTTCTGATGTTGCTATTGGCACCACTCCTGGCACGAATACCGCAACTGCAGATGAAGACGAGTCCAAACCTTTTGCCCAGCGGAAATCAAGCAAAGATTCCTTCCCTCACAAAATAGAAGAAAGCCATCAGGTTGAGAGCTTGTCAGAGAAACTCTTCCTACCTATAAGTCTGAGTGAGGCTACAGCCAAGGACTCCGAGGAGAACCACCACGCAGATACAGAGGATGAGCCCCCTGACATGTCATCCAAATCGAACAAGAGCAAGAAGAAGAAAAAGAAAGAGAAGAAGAAACATGATACTAAGTCTAAGAAAAGTTCCAAGAAAAACATTGTGAAGATCTCAATGCTCGATGACAATGTTTTCCTAAATGATCAGCCTGAGGTTACGCAAGAGGAAGAGGATGAAGAAGAGGAAGAAGAAGGTACAAAGTTCTTCCAGATTTCTAGCCTTGAGTCTTTTTTTCTTTCTGGCTCAAAAGAAGAACCAAATGAGCTCATAAAAACCGAGCCCCTTGTCGCAGAGGAAGATATGAAAGAAGAAACCGATGTGGGTAGTCTAAAAAGGACGTCATCTTGCGAAGATATCAAGGTAGCAACTAACGATGACTATTTGGAAGAGCCTGTAAAAGAAGAGGCTTCCACAGATTACATGGAGCCTGAAGACCCGCTTCCTACCAAGGTTGAAGAGCCACCAAAGGAAGAACTGGAGGCTGAAGAAGATGAAGTTGAAGAGGAAATCGAGGAGATCAATAAACCAAATACCGACTCCTCCCCGCTGGAAGAAGACAACTCCAACAATTTAGAAAAATGTGAAAGTTTTCAAGATGACGTCTACGATTTTGACAAGTCGGAAAATATCACTGAGAAACTGACTGATTGTGAACAGTTGAAGAAATTTGAAAAAATCTACGACTCTGAGAAGGTGAAAAAGGTTTCTGATGAAGCTGTCGATTCTGACATTCCAAAAAGTGTTGGAGAAGAGTCCCCTGAGTTCCTCTTATCACCAGGAGAGGAGGTCACTGAGGAAGACCCTGATCAACCAGAACCCACACCTGAAGAGCAGTCAGATAACGATGTTCAAAAAGAAACAACAGACGATATTCAAGCTCAAAGCCCTGAGCACTCAGATACAAAAGAGTTGGACAGTCCTTATCAGTTAAAAGAACAAGAAGATACTGTATCGGATGTGTTTACTGATTTGAAGAAACTCGATAAAAGCTGCTCTGAAGTTGAAGACAACAGTGAAGATAAGATTGTTGCAAAAAAGCCTGAAAATGAGGATTTACTTGAAAAAGAAGAACCAGAAAAGACTGAAGAAGAATCGAAATCTCCTGAAACCTCAGAAGTTGAAAAAGATCCGACAGACGAATCTTGTGAGACAACTGTAGATGGTGAGGATCTGAAACGAGATGAAGGCAAAGATCAAAATTCCAGCAGTTTTGAGAAACTAGCCGAAACCAACGAAAAGGTAATTGGTATGTTCAACTTTGGCCAATTACAAGTTCATGTTGTTGACAAGAATACAGGTAACCTTCACAAGGTGCATGACGCCGAGTTATCCGATAAAGATAGTCGGGCCGCTCTTGGTGTGGATCAGGTTCTTGACAAAGCTAACAGCAACACTGATGCCAACAGGAATGTAGACCACATCTATGACTTTGACGATACCACTGAAGAAATGTACGGTAATACTTATCTCGATAAAGGAGATAGCAAATCCAGACAGTTCAACTACAAAAAGGTCTTTATGATGCAACCCCGTAAATCGGACAGAATTTATAATTCTGATGTCACCAAGAAAATCTGCGAGAAGATTGAGCAGAAGATTTACAGCCTAGAAACGCCAACCAAATTTGATACCATCTTCAAGAAAATTAACAAAAAGGCGTACAATATGGACGAGGCTAATTCTACTGAAGATAAGATTTATGAGGATGAGATTTGCAAGACGCGCAATCGAATTGCAGATCTCCACCACCCGCAGCAGCATACAGATAAAACAGATGGTGACGTCTATGATTTCGATGAACAAGACAAGCTTGAAGTTAATGAATTGTTAGATTCCAAGAAAACGGATCATATAATCCATGATCTTCACGATTACAAGAAAACCGAAGAGCATTTACCCAATTTTGATGGAGGCAAGAAAAGTACCGAAAAACTGGCAAGGTTCCAAGATTGGTTTACAAATGACGATAAGATATATTACTTTGACCCTGTTAAAGATTTTGGTGCCAAGCCTAACAAAAGCCGTTCCAAATGTGGTAAACACAAATCTGGTCAGGGGAAACGCCGCAATAAAGGTAAAATTGTTACCTTTACTTGGCCTCGTGAACTAGACTTAAAATCATCGGTTGAATCTGCAACGGACAATGATACTGTCACTGCCGATACGTTTGAGATGAAAGAGAGCGATTATGCTGAACTCGATGACCAAGAAGATAAACTTGTCATCGATACAAGTGACATGTGCGTCTATGATGAAAATTCAACCGAACCGATTGAAAAAGAGGCTGATAAAGCAAAGCTGTGTAGCTCGGAAACAAAAACCGAAAACGATTTAACAGCCAATGATGAAAAAAACGAAGATATCCAAGACGACGTGAAAGAAAAGGTGTCAAGTGAAGAGGAAGAAAAAGAGCAAGAAGAAACGGAGGCTCTGGATACCACGTCTTTTAAAAATAATCAAGACACTTTTGGAATGTCGATCAAACAGAAAGTGGTGGTGAATTATGGACTAACAAATCACAGTGGCACAACCCCCGACAAGATAGAATCGGTTATAGGCTTTGGAATGCTATCTGGCAACAAAGACCTGCTGTCGGCATCTAAATTGTCTGATGACTGTGACAGCAGCGATGAAGGTGCGCTCAGAATCGATGAAGATGAAGTAGACGAAAATACTGATAGGGATGAAATGACACTTGAAGACAAAACCCCAGAAAAAGAGAAAGAAATTGAAGAAGAAGAAAAAGTGATCAAGAAGCCTGAACCCTATGAAGAAGAGGAGGAGGAAGAAGC'
#sequence = 'AAACGCTAGAAACCATAGAAACCCACTATC'
missed_cleavages = 1
min_length = None
#a2g_sites = [1,7,18,19,20,30,34,38,45,57,62]
#c2t_sites = [23,25,41,51,53]
a2g_sites = [0, 1, 2, 7, 10, 11, 14, 16, 18, 20, 24, 27]
c2t_sites = [3, 5, 12, 13, 21, 22, 23, 25, 29]
#a_to_g_editings = a_list
#c_to_t_editings = c_list
look_ahead_for_editnig = 3
look_behind_for_editing = 6


def create_edited_rna_peptides(sequence,a2g_sites,c2t_sites,digestion_rule, missed_cleavages = 0,
                               look_ahead_for_editnig=3,look_behind_for_editing=6):

    rna_peptides_list = []

    cleavage_sites, binary_for_fixed_cs = create_all_cleavage_sites(sequence,digestion_rule,a2g_sites,c2t_sites,look_ahead_for_editnig,look_behind_for_editing)
    
#    print(list(cleavage_sites))
#    print(list(binary_for_fixed_cs))
    
    cleavage_sites_lengeth = len(cleavage_sites)
    seq_length = len(sequence)
    
    #running over all posiblle windows from a site until mc+1 fixed sites ahead
    for i in range(cleavage_sites_lengeth-1):
        
        fixed_sites_cnt = 0
        j=1
        
        #running over all j (1+j -> index of cleavage site at the end of the subsequence window) untill betwin i-th cs and (i+j)-th cleavage site there are mc number of fixed site or i+j site is sequence end
        #in this method all possible windows for subsequences are examined
        while (fixed_sites_cnt < missed_cleavages+1) and (cleavage_sites_lengeth-i-1>=j):
            
            append_zero=False, 
            append_seq_end = False
#            print('new_range')
#            print(str(i) + ',' + str(j) + '  _' + str(fixed_sites_cnt))
            if binary_for_fixed_cs[i+j]: #count site if site ahead is fixed.
                fixed_sites_cnt+=1
            
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
            
            permutation_range = (edit_range_start,edit_range_end-1)
            sebsequence_in_perm_range = sequence[edit_range_start:edit_range_end]
            
            subsequence_coor = [cleavage_sites[i],cleavage_sites[i+j]-1]
            subsequence = sequence[subsequence_coor[0]:subsequence_coor[1]+1]
            
            #all eiditing sites permutation in permutation_range
            permutations = permutations = find_permutations_in_range(creat_editing_sites_permutations_one_type(find_sites_in_window(permutation_range,a2g_sites)),
                                                                     creat_editing_sites_permutations_one_type(find_sites_in_window(permutation_range,c2t_sites)))
            
            #determining if sequence beginning or end is in cleavage_sites
            if edit_range_start == 0:
                append_zero = True
            if edit_range_end == seq_length:
                append_seq_end = True
                
            for perm in permutations:
                
                edited_sebsequence_to_edit = edit_rna_as_peptide(sebsequence_in_perm_range,permutation_range,perm[0],perm[1])
#                print(edited_sebsequence_to_edit)
                
                #the edited cleavage sites withing the range of subsequence (including_bounderies)
                perm_cs = [k+edit_range_start for k in find_cleavage_sites(digestion_rule,edited_sebsequence_to_edit,append_zero=append_zero, append_seq_end = append_seq_end) if (k+edit_range_start>=cleavage_sites[i] and k+edit_range_start<=cleavage_sites[i+j])]
#                print(perm_cs)
                
                #if bouderies still exists as cleavage sites and no more than mc cleavage sites in perm_cs - a valid seq
                if set([cleavage_sites[i], cleavage_sites[i+j]]).issubset(perm_cs) and len(perm_cs)<=missed_cleavages+2:
                    edited_sub_seq = edit_rna_as_peptide(subsequence,subsequence_coor,
                                                         find_sites_in_window(subsequence_coor,perm[0]),
                                                         find_sites_in_window(subsequence_coor,perm[1]))
                    
                    rna_pep_perm = rna_peptide(edited_sub_seq,subsequence_coor,a2g = perm[0], c2t = perm[1], cleavage_sites = perm_cs)
                    rna_peptides_list.append(rna_pep_perm)
#                    rna_pep_perm.print_peptide_data()            
            j+=1
            
    return rna_peptides_list
            



"""
given editing coordinations and a sequence and a lookahead lookbehind range
run over all codons, edit and if codos could be a cleavage site append to list.
create another binary list - for each codon fixed or mutable
"""
#a,b = create_all_cleavage_sites(sequence,digestion_rule,a2g,c2t,look_ahead_for_editnig,look_behind_for_editing)
def create_all_cleavage_sites(sequence,digestion_rule,a2g,c2t,look_ahead_for_editnig,look_behind_for_editing):
    
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
            permutations = find_permutations_in_range(creat_editing_sites_permutations_one_type(find_sites_in_window(coor,a2g)),
                                                      creat_editing_sites_permutations_one_type(find_sites_in_window(coor,c2t)))
            
#            print(seq)
#            print(permutations)
            
            is_cs, is_fixed = determine_site(seq,coor,digestion_rule,codon_to_asses_in_window,permutations)
            
            if is_cs:
                cleavage_sites.append(i)
                binary_for_fixed_sites.append(is_fixed)
                
        cleavage_sites.append(length)
        binary_for_fixed_sites.append(1)
#        print(list(cleavage_sites))
#        print(list(binary_for_fixed_sites))
        
        return cleavage_sites,binary_for_fixed_sites
            
                
"""
for a subsequence determine codon_to_asses is a cleavage site and if so if mutable
"""        
#a,b = determine_site('AAACGC',[0,5],digestion_rule,3,permutations)
def determine_site(subsequence,coor,digestion_rule,codon_to_asses,permutations):
    
    cleavage_site = set()
    
    for perm in permutations:
        seq = edit_rna_as_peptide(subsequence,coor,perm[0],perm[1])
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
       
            
"""
find all in-frame cleavage sites given a sequence and a digestion_rule
also append 0 and sequence_length if needed
"""
#k = find_cleavage_sites(digestion_rule,sequnce,append_zero=False, append_seq_end = False)
def find_cleavage_sites(digestion_rule,sequnce,append_zero=False, append_seq_end = False):
    
    if append_zero:
        initial_list = [0] + [x.end() for x in digestion_rule.finditer(sequnce)]
    else:
        initial_list = [x.end() for x in digestion_rule.finditer(sequnce)]
    initial_set = set(initial_list)

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

  

"""
return all combinations of elements from to lists
"""      
#k = find_permutations_in_range(creat_editing_sites_permutations_one_type(a2g_sites),creat_editing_sites_permutations_one_type(c2t_sites))
def find_permutations_in_range(a_to_g_permutations,c_to_t_permutations):
    return[(l[0],l[1]) for l in it.product(a_to_g_permutations,c_to_t_permutations)]
    
   
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



#pep_list = create_edited_rna_peptides(sequence,a2g_sites,c2t_sites,digestion_rule, missed_cleavages = missed_cleavages, look_ahead_for_editnig=look_ahead_for_editnig, look_behind_for_editing=look_behind_for_editing)
        
        
        

        

    
