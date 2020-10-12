import pandas as pd
import numpy as np
import pandas as pd
import re
import functools
from Bio import SeqIO 
from edited_peptides_from_seq import create_edited_rna_peptides, all_mm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


#sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/' +'test3.txt', 'w')
#digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
#missed_cleavages = 2
#look_ahead_for_editnig = 3
#look_behind_for_editing = 6
#min_length = None
#max_length = None
#max_sites_per_pep = None
#max_mass = 4600
#input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
#input_fasta = 'sim.fasta'
#input_fasta = 'in_frame_rna_from_orfs_squ.fasta'
#sequence = 'ATGGCTGCTTTGAGACAGCCCCAGGTCGCGGAGCTGCTGGCCGAGGCCCGGCGAGCCTTCCGGGAGGAGTTCGGGGCCGAGCCCGAGCTGGCCGTGTCAGCGCCGGGCCGCGTCAACCTCATCGGGGAACACACGGACTACAACCAGGGCCTGGTGCTGCCTATGGCTCTGGAGCTCATGACGGTGCTGGTGGGCAGCCCCCGCAAGGATGGGCTGGTGTCTCTCCTCACCACCTCTGAGGGTGCCGATGAGCCCCAGCGGCTGCAGTTTCCACTGCCCACAGCCCAGCGCTCGCTGGAGCCTGGGACTCCTCGGTGGGCCAACTATGTCAAGGGAGTGATTCAGTACTACCCAGCTGCCCCCCTCCCTGGCTTCAGTGCAGTGGTGGTCAGCTCAGTGCCCCTGGGGGGTGGCCTGTCCAGCTCAGCATCCTTGGAAGTGGCCACGTACACCTTCCTCCAGCAGCTCTGTCCAGACTCGGGCACAATAGCTGCCCGCGCCCAGGTGTGTCAGCAGGCCGAGCACAGCTTCGCAGGGATGCCCTGTGGCATCATGGACCAGTTCATCTCACTTATGGGACAGAAAGGCCACGCGCTGCTCATTGACTGCAGGTCCTTGGAGACCAGCCTGGTGCCACTCTCGGACCCCAAGCTGGCCGTGCTCATCACCAACTCTAATGTCCGCCACTCCCTGGCCTCCAGCGAGTACCCTGTGCGGCGGCGCCAATGTGAAGAAGTGGCCCGGGCGCTGGGCAAGGAAAGCCTCCGGGAGGTACAACTGGAAGAGCTAGAGGCTGCCAGGGACCTGGTGAGCAAAGAGGGCTTCCGGCGGGCCCGGCACGTGGTGGGGGAGATTCGGCGCACGGCCCAGGCAGCGGCCGCCCTGAGACGTGGCGACTACAGAGCCTTTGGCCGCCTCATGGTGGAGAGCCACCGCTCACTCAGAGACGACTATGAGGTGAGCTGCCCAGAGCTGGACCAGCTGGTGGAGGCTGCGCTTGCTGTGCCTGGGGTTTATGGCAGCCGCATGACGGGCGGTGGCTTCGGTGGCTGCACGGTGACACTGCTGGAGGCCTCCGCTGCTCCCCACGCCATGCGGCACATCCAGGAGCACTACGGCGGGACTGCCACCTTCTACCTCTCTCAAGCAGCCGATGGAGCCAAGGTGCTGTGCTTG'
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


def edited_peptides_from_rna_seq(input_path,input_fasta, mm_headers,
                                 digestion_rule, missed_cleavages, look_ahead_for_editnig = 3 ,look_behind_for_editing = 6,
                                 min_length = None, max_mass = None, max_sites_per_pep = None):

    
    peps_dfs_list = []
    under_min_length_dict = {}
    over_max_sites_dict = {}
    over_mass_dict = {}
    over_mass_combs_dict = {}
    n=0
    tot_editing_sites = 0
    
    #iterate over all sequences in input_fasta file containing in-frame mRNA's
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        sites_dict = {}
        unique_peptides = 0
        over_mass_combs = 0
        n+=1
        [sites_dict.update({mm:sorted(eval(find_by_regex_in_header(record.description,mm_headers[mm])))}) for mm in all_mm]

        #get editing sites lists from record header
        tot_editing_sites += sum([len(sites_dict[mm]) for mm in all_mm])

        #create all peptides from all editing combinations
        rna_peptides_list, under_min_length, over_max_sites_peps, over_mass_peps, over_mass_combs_per_seq_dict = create_edited_rna_peptides(str(record.seq),sites_dict,digestion_rule, missed_cleavages = missed_cleavages,
                                                                                                                                            look_ahead_for_editnig=look_ahead_for_editnig, look_behind_for_editing=look_behind_for_editing,
                                                                                                                                            min_length=min_length, max_mass = max_mass, max_sites_per_pep = max_sites_per_pep)
        if len(under_min_length):
            under_min_length_dict.update({record.id:under_min_length})
        if len(over_max_sites_peps):
            over_max_sites_dict.update({record.id:over_max_sites_peps})
        if len(over_mass_peps):
            over_mass_dict.update({record.id:over_mass_peps})
        if len(over_mass_combs_per_seq_dict):
            over_mass_combs_dict.update({record.id:over_mass_combs_per_seq_dict})
            over_mass_combs = sum([len(over_mass_combs_per_seq_dict[coor]) for coor in over_mass_combs_per_seq_dict])    
        
        unattended_peps = len(under_min_length_dict) + len(over_max_sites_peps) + len(over_mass_peps)
        number_of_rna_peptides = len(rna_peptides_list)
            
        if number_of_rna_peptides: #for each sequence (if peptides are generated from it) create a dataframe with multiindex (pep,seq,coor)
            #list of tuples to create initial dataframe from
            pep_data_list = [(str(Seq(pep.seq, generic_dna).translate()), str(Seq(pep.seq, generic_dna).translate()).replace('I','X').replace('L','X'), str(Seq(pep.extended_seq, generic_dna).translate()), str(record.id), '_'+str(pep.coo[0])+'_'+str(pep.coo[1]), '_'+str(pep.perm_coo[0])+'_'+str(pep.perm_coo[1]), str(pep.edited_sites), pep.N_terminus, pep.C_terminus, pep.sites_in_pep_range, pep.cancelled_cs_in_pep) for pep in rna_peptides_list]
            del(rna_peptides_list)
            #initial dataframe
            peps_df_temp = pd.DataFrame(pep_data_list,columns = ['biological_peptide','peptide','biological_extended_peptide','seq_id','in_frame_coordinates_base0','permutation_coor_base0','editing_comb', 'N_terminus', 'C_terminus','sites_in_permutation_range','cancelled_cs_in_pep'])
            #groupby multiindex and create a list of all editing combinations for grouped rows
            peps_df = peps_df_temp.groupby(['peptide','seq_id','in_frame_coordinates_base0'])['editing_comb'].aggregate(lambda x: [eval(c) for c in list(x)])
            #set name of new column in grouped df
            peps_df = pd.DataFrame({'editing_combinations_relative_to_coding_seq_base0': peps_df.values}, index = peps_df.index)
            #create another column of editing combinations intersection - sites that are obviosly edited
            peps_df['editing_combinations_intersection_base0'] = peps_df['editing_combinations_relative_to_coding_seq_base0'].apply(infer_editing_comb)
#            print(peps_df['editing_combinations_relative_to_coding_seq_base0'])
#            print(peps_df['editing_combinations_intersection_base0'])
            #retrive N\C terminus change from pepd_df_temp
            peps_df_temp = peps_df_temp.loc[:,('biological_peptide','peptide','biological_extended_peptide','seq_id','in_frame_coordinates_base0','permutation_coor_base0','N_terminus','C_terminus','sites_in_permutation_range','cancelled_cs_in_pep')]
            peps_df_temp.drop_duplicates(['peptide','seq_id','in_frame_coordinates_base0'],inplace = True)
            peps_df = peps_df.reset_index()
            peps_df = pd.merge(peps_df,peps_df_temp, on = ['peptide','seq_id','in_frame_coordinates_base0'])
#            print(peps_df['editing_combinations_intersection_base0'])
            peps_df['editing_combinations_intersection_base1'] = peps_df.apply(lambda row: tuple([[site+1 for site in row['editing_combinations_intersection_base0'][i]] for i in range(len(row['editing_combinations_intersection_base0']))]) , axis = 1)
            #append df to dataframes list
            peps_dfs_list.append(peps_df)
            unique_peptides = len(set(list(peps_df['peptide'])))    
#           store.append('all_peps', peps_df)
                   
        print('recoed ' + str(n) + ', ' + record.id + ' | rna peps: ' + str(number_of_rna_peptides) + ' | unique peps: ' + str(unique_peptides) + ' | unattended peps: ' + str(unattended_peps) + ' | over mass combs (if pep attended): ' + str(over_mass_combs))
        
    #concat all dataframes and sort by index 
    print('\nDigestion and editing of Sequnces finished')        
    print('Concatenating tabels from all sequences')
    final_peps_df = pd.concat(peps_dfs_list)
    del(peps_dfs_list)
    print('Sorting final tabel')
    final_peps_df.sort_index(inplace = True)
    return final_peps_df, under_min_length_dict, over_max_sites_dict, over_mass_dict, over_mass_combs_dict, tot_editing_sites

    
"""
find regex (regex - pre compiled regex) in header
"""
def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return []
    


def infer_editing_comb(editing_combis):
    """
    find intersection of all editing combinations represented by peptide
    """
    return tuple([sorted(list(set.intersection(*map(set,[x[i] for x in editing_combis])))) for i in range(len(all_mm))])

def get_editing_sites(comps_inf_peps_df):
    """
    from dataframe containing list of inferred editing sites per peptide
    create a nested dictionary of record_id:{a2g_sites_per_comp:sites, c2t_detectable_sites:sites}
    and retrive the count of all sites represented by peptides in the input dataframe
    """
    comps_inf_peps_df = comps_inf_peps_df.reset_index().set_index(['seq_id','peptide','in_frame_coordinates_base0'])
    detectable_sites_dict = {}
    detectable_sites_num = 0
    
    for comp, comp_df in comps_inf_peps_df.groupby(level=0):
        temp_all_detectable_sites_per_seq = [[],[],[],[],[],[],[],[],[],[],[],[]]
        number_of_mm_types = len(temp_all_detectable_sites_per_seq)
        
        for index, row in comp_df.iterrows():
            for i in range(number_of_mm_types):
                temp_all_detectable_sites_per_seq[i] += row['editing_combinations_intersection_base1'][i]
            
        all_detectable_sites_per_seq = ()
        for l in temp_all_detectable_sites_per_seq:
            s = list(set(l))
            all_detectable_sites_per_seq += (s,)
            detectable_sites_num+=len(s)
            
        detectable_sites_dict.update({comp:all_detectable_sites_per_seq})
    
#    print(str(detectable_sites_num) + ' editing sites are detectable by at least one informative peptide')
    return detectable_sites_dict, detectable_sites_num 



def process_peptides_df(peps_df):
    """
    get peps_df with arbitrary index columns (no important data in index)
    and components dataframe (contain data of orfs etc)
    and return final_peps_df. for each pep determine if informative, for wich sites and if unedited version of peptide is also informative
    """
    
    peps_df = peps_df.set_index('peptide') #index is just peptide column
    peptide_occurrences = peps_df.groupby(level=0).size().reset_index(name = 'peptide_count').set_index('peptide') #count of occurances per peptide
    peps_df = peps_df.join(peptide_occurrences)
    peps_df['informative_peptide'] = peps_df.apply(lambda x: True if x.peptide_count == 1 else False, axis = 1) #peptide is informative if occure in only one location in one sequence
    peps_df['edited'] = peps_df.apply(lambda x: True if x.editing_combinations_intersection_base1 != ([],[],[],[],[],[],[],[],[],[],[],[]) else False, axis = 1) #determine if peptide is edited based on editing combinations data
    
    #detemining if unedited version is detectable as well
    peps_df = peps_df.reset_index()
    peps_df.set_index(['seq_id','in_frame_coordinates_base0'], inplace = True)
    unedited_peptides =  peps_df[~peps_df['edited']]['informative_peptide'].reset_index(name = 'informative_original_version').set_index(['seq_id','in_frame_coordinates_base0'])  #original peptides (unedited) df
    peps_df = peps_df.join(unedited_peptides)
    peps_df.fillna('no_original_version',inplace = True)
 
    peps_df['editing_sites_inferred'] = peps_df.apply(lambda row: sum([len(x) for x in row['editing_combinations_intersection_base0']]), axis = 1)
    
    peps_df = peps_df.reset_index()
    peps_df = peps_df.sort_values('seq_id')

    return peps_df
