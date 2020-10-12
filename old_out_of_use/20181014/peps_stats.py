import numpy as np
import pandas as pd
import re
import functools
from Bio import SeqIO 

#strand_regex = re.compile('(?<=strand:\s)[^\s]+')
#prot_start_nuc_regex = re.compile('(?<=prot_start_nuc:\s)[^\s]+')
#prot_end_nuc_regex = re.compile('(?<=prot_start_nuc:\s)[^\s]+')
#prot_start_regex = re.compile('(?<=prot_start:\s)[^\s]+')
#prot_end_regex = re.compile('(?<=prot_end:\s)[^\s]+')
#original_orf_start_regex = re.compile('(?<=original_orf_start:\s)[^\s]+')
#original_orf_end_regex = re.compile('(?<=original_orf_end:\s)[^\s]+')
#input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
#input_fasta = 'in_frame_rna_rec_only_from_orfs_squ.fasta'


def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'


"""
create dataframe from input fasta file whre indeces are components id's and rows are description data
"""
#comps_df = get_components_df(input_path, input_fasta, strand_regex, prot_start_nuc_regex, prot_end_nuc_regex, prot_start_regex, prot_end_regex)
def get_components_df(input_path, input_fasta, strand_regex, prot_start_nuc_regex, prot_end_nuc_regex, prot_start_regex, prot_end_regex,original_orf_start_regex,original_orf_end_regex):

    comps_df = pd.DataFrame(index = ['component'], columns = ['strand','prot_start_nuc_base1','prot_end_nuc_base1','prot_start','prot_end','original_orf_start','original_orf_end'])
    
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        strand = find_by_regex_in_header(record.description, strand_regex)
        prot_start_nuc = int(find_by_regex_in_header(record.description, prot_start_nuc_regex))
        prot_end_nuc = int(find_by_regex_in_header(record.description, prot_end_nuc_regex))
        prot_start = find_by_regex_in_header(record.description, prot_start_regex)
        prot_end = find_by_regex_in_header(record.description, prot_end_regex)
        original_orf_start = int(find_by_regex_in_header(record.description, original_orf_start_regex))
        original_orf_end = int(find_by_regex_in_header(record.description, original_orf_end_regex))
        
        comps_df.loc[record.id] = pd.Series({'strand':strand,'prot_start_nuc_base1':prot_start_nuc,'prot_end_nuc_base1':prot_end_nuc,'prot_start':prot_start,'prot_end':prot_end,'original_orf_start':original_orf_start,'original_orf_end':original_orf_end})
        
    return comps_df


"""
used by process_peptides_df (and detectable_editing_sitse_df) to calculate real coordinates of infered editing sites for each row
based on data joined from get_components_df results into the informative peptides datafrme (with index as component)
editing_type - index of list of site of editing type whithin wditing combinations represented by peptide:
0 - a2g 
1 - c2t
"""
def get_inferred_editings_as_list(infered_editing_comb, strand, prot_start_nuc, prot_end_nuc, editing_type = 0):
    
    infered_editing_sites_list = [site for site in infered_editing_comb[editing_type]]
    
    if strand == '+':
        infered_editing_sites_list = [int(site) + int(prot_start_nuc) for site in infered_editing_sites_list]
    elif strand == '-':
        infered_editing_sites_list = [int(prot_end_nuc) - int(site) for site in infered_editing_sites_list]
    
    return infered_editing_sites_list
       

"""
from comps dataframe containing informative peptides
get a dictionary where keys are components id and values are lists of detectable sites
"""
#det_sites_dict = count_edited_sites(comps_inf_peps_df)
def get_editing_sites(comps_inf_peps_df):
    
    comps_inf_peps_df = comps_inf_peps_df.reset_index().set_index(['seq_id','peptide','in_frame_coordinates_base0'])
    detectable_sites_dict = {}
    detectable_sites_num = 0
    
    for comp, comp_df in comps_inf_peps_df.groupby(level=0):
        for index, row in comp_df.iterrows():
            a2g_sites_per_comp = list(set([site for site in row['real_inferred_a2g_sites_coor_base1']]))
            c2t_sites_per_comp = list(set([site for site in row['real_inferred_c2t_sites_coor_base1']]))
            
        detectable_sites_dict.update({comp:{'a2g_detectable_sites':a2g_sites_per_comp,'c2t_detectable_sites':c2t_sites_per_comp}})
        detectable_sites_num += len(a2g_sites_per_comp) + len(c2t_sites_per_comp)
    
#    print(str(detectable_sites_num) + ' editing sites are detectable by at least one informative peptide')
    return detectable_sites_dict, detectable_sites_num 


"""
get peps_df with arbitrary index columns (no important data in index)
and components dataframe (contain data of orfs etc)
and return final_peps_df. for each pep determine if informative, for wich sites and if unedited version of peptide is also informative
"""
def process_peptides_df(peps_df, comps_df):
    
    peps_df = peps_df.set_index('peptide') #index is just peptide column
    peptide_occurrences = peps_df.groupby(level=0).size().reset_index(name = 'peptide_count').set_index('peptide') #count of occurances per peptide
    peps_df = peps_df.join(peptide_occurrences)
    peps_df['informative_peptide'] = peps_df.apply(lambda x: True if x.peptide_count == 1 else False, axis = 1) #peptide is informative if occure in only one location in one sequence
    peps_df['edited'] = peps_df.apply(lambda x: True if x.inferred_editing_combination_relative_to_sense_orf_base0 != ([],[]) else False, axis = 1) #determine if peptide is edited based on editing combinations data
    
    #detemining if unedited version is detectable as well
    peps_df = peps_df.reset_index()
    peps_df.set_index(['seq_id','in_frame_coordinates_base0'], inplace = True)
    unedited_peptides =  peps_df[~peps_df['edited']]['informative_peptide'].reset_index(name = 'informative_original_version').set_index(['seq_id','in_frame_coordinates_base0'])  #original peptides (unedited) df
    peps_df = peps_df.join(unedited_peptides)
    peps_df.fillna('no_original_version',inplace = True)

    #joining colmns from components_dataframe before real sites calculation
    peps_df = peps_df.reset_index().set_index('seq_id')
    peps_df = peps_df.join(pd.DataFrame(comps_df['strand']))
    peps_df = peps_df.join(pd.DataFrame(comps_df['prot_start_nuc_base1']))
    peps_df = peps_df.join(pd.DataFrame(comps_df['prot_end_nuc_base1']))
    peps_df = peps_df.join(pd.DataFrame(comps_df['prot_start']))
    peps_df = peps_df.join(pd.DataFrame(comps_df['prot_end']))
    
    #calculate real coordinates of editing sites inffered by each peptide (relative to original rna seq)
    peps_df['real_inferred_a2g_sites_coor_base1'] = peps_df.apply(lambda x: get_inferred_editings_as_list(x.inferred_editing_combination_relative_to_sense_orf_base0, x.strand, x.prot_start_nuc_base1, x.prot_end_nuc_base1, editing_type = 0), axis = 1)
    peps_df['real_inferred_c2t_sites_coor_base1'] = peps_df.apply(lambda x: get_inferred_editings_as_list(x.inferred_editing_combination_relative_to_sense_orf_base0, x.strand, x.prot_start_nuc_base1, x.prot_end_nuc_base1, editing_type = 1), axis = 1)
    #calculate number of inffered sites per peptide 
    peps_df['editing_sites_inferred'] = peps_df.apply(lambda x: len(x.real_inferred_a2g_sites_coor_base1) + len(x.real_inferred_c2t_sites_coor_base1), axis = 1)
    
    peps_df = peps_df.rename_axis('seq_id').reset_index()
    peps_df = peps_df.sort_values('seq_id')

    return peps_df

#
#peps_df = peps_df.set_index('peptide') #index is just peptide column
#peptide_occurrences = peps_df.groupby(level=0).size().reset_index(name = 'peptide_count').set_index('peptide') #count of occurances per peptide
#peps_df = peps_df.join(peptide_occurrences)
#peps_df['informative_peptide'] = peps_df.apply(lambda x: True if x.peptide_count == 1 else False, axis = 1) #peptide is informative if occure in only one location in one sequence
#peps_df['edited'] = peps_df.apply(lambda x: True if x.infered_editing_comb != ([],[]) else False, axis = 1) #determine if peptide is edited based on editing combinations data
#
##detemining if unedited version is detectable as well
#peps_df = peps_df.reset_index()
#peps_df.set_index(['seq_id','in_frame_coordinates_base0'], inplace = True)
#unedited_peptides =  peps_df[~peps_df['edited']]['informative_peptide'].reset_index(name = 'informative_original_version').set_index(['seq_id','in_frame_coordinates_base0'])  #original peptides (unedited) df
#peps_df = peps_df.join(unedited_peptides)
#peps_df.fillna('no_original_version',inplace = True)
#
##joining colmns from components_dataframe before real sites calculation
#peps_df = peps_df.reset_index().set_index('seq_id')
#peps_df = peps_df.join(pd.DataFrame(comps_df['strand']))
#peps_df = peps_df.join(pd.DataFrame(comps_df['prot_start_nuc']))
#peps_df = peps_df.join(pd.DataFrame(comps_df['prot_end_nuc']))
#peps_df = peps_df.join(pd.DataFrame(comps_df['prot_start']))
#peps_df = peps_df.join(pd.DataFrame(comps_df['prot_end']))



# =============================================================================
# def get_informative_peptides(peps):
#     return peps[pd.DataFrame(index = peps.index.get_level_values(0), data = peps.groupby(level = 0).size()).values == 1]
#    
# def get_uninformative_peptides(peps):
#     return peps[pd.DataFrame(index = peps.index.get_level_values(0), data = peps.groupby(level = 0).size()).values > 1]
# 
# 
# def create_detectable_editing_sites_df(peps_df, comps_additional_data):
#     
#     #get informative peptide dataframe
#     inf_peps = get_informative_peptides(peps_df)
#     #set index to seq_id
#     comps_inf_peps_df = inf_peps.reset_index().set_index('seq_id')
#     comps_inf_peps_df.sort_index(inplace = True)
#     
#     #getting only rows representing peptides containing editing sites
#     comps_inf_peps_df = comps_inf_peps_df[comps_inf_peps_df['infered_editing_comb'] != ([],[])]
#     
#     #joining colmns from components_dataframe before real sites calculation
#     comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['strand']))
#     comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['prot_start_nuc']))
#     comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['prot_end_nuc']))
#     comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['prot_start']))
#     comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['prot_end']))
#     
#     #calculate real coordinates (relative to original rna seq)
#     comps_inf_peps_df['real_infered_sites_coor_base1'] = comps_inf_peps_df.apply(lambda x: get_infered_editings_as_list(x.infered_editing_comb, x.strand, x.prot_start_nuc, x.prot_end_nuc), axis = 1)
# 
#     comps_inf_peps_df = comps_inf_peps_df.reset_index().set_index(['index','peptide','in_frame_coordinates_base0'])
#     
#     return comps_inf_peps_df
# =============================================================================

