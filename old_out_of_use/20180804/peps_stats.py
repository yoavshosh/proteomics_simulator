import numpy as np
import pandas as pd
import re
import functools
from Bio import SeqIO 

strand_regex = re.compile('(?<=strand:\s)[^\s]+')
start_nuc_regex = re.compile('(?<=start_nuc:\s)[^\s]+')
end_nuc_regex = re.compile('(?<=end_nuc:\s)[^\s]+')
prot_beginning_regex = re.compile('(?<=prot_beginning:\s)[^\s]+')
prot_end_regex = re.compile('(?<=prot_end:\s)[^\s]+')
input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
input_fasta = 'in_frame_rna_rec_only_from_orfs_squ.fasta'


def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'

def get_informative_peptides(peps):
    return peps[pd.DataFrame(index = peps.index.get_level_values(0), data = peps.groupby(level = 0).size()).values == 1]
    
def get_uninformative_peptides(peps):
    return peps[pd.DataFrame(index = peps.index.get_level_values(0), data = peps.groupby(level = 0).size()).values > 1]



"""
create dataframe from input fasta file whre indeces are components id's and rows are description data
"""
#comps_df = get_components_df(input_path, input_fasta, strand_regex, start_nuc_regex, end_nuc_regex, prot_beginning_regex, prot_end_regex)
def get_components_df(input_path, input_fasta, strand_regex, start_nuc_regex, end_nuc_regex, prot_beginning_regex, prot_end_regex):

    comps_df = pd.DataFrame(index = ['component'], columns = ['strand','start_nuc','end_nuc','prot_beginning','prot_end'])
    
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        strand = find_by_regex_in_header(record.description, strand_regex)
        start_nuc = int(find_by_regex_in_header(record.description, start_nuc_regex))
        end_nuc = int(find_by_regex_in_header(record.description, end_nuc_regex))
        prot_beginning = find_by_regex_in_header(record.description, prot_beginning_regex)
        prot_end = find_by_regex_in_header(record.description, prot_end_regex)
        
        comps_df.loc[record.id] = pd.Series({'strand':strand,'start_nuc':start_nuc,'end_nuc':end_nuc,'prot_beginning':prot_beginning,'prot_end':prot_end})
        
    return comps_df


"""
used by detectable_editing_sitse_df to calculate real coordinates of infered editing sites for each row
based on data joined from get_components_df results into the informative peptides datafrme (with index as component)
"""
def get_infered_editings_as_list(infered_editing_comb, strand, start_nuc, end_nuc):
    
    infered_editing_sites_list = []
    
    for editing_type in infered_editing_comb:
        [infered_editing_sites_list.append(site) for site in editing_type]
    
    if strand == '+':
        infered_editing_sites_list = [int(site + start_nuc) for site in infered_editing_sites_list]
    elif strand == '-':
        infered_editing_sites_list = [int(end_nuc - site) for site in infered_editing_sites_list]
    
    return infered_editing_sites_list
 
    
"""

"""
def create_detectable_editing_sites_df(peps_df, comps_additional_data):
    
    #get informative peptide dataframe
    inf_peps = get_informative_peptides(peps_df)
    #set index to seq_id
    comps_inf_peps_df = inf_peps.reset_index().set_index('seq_id')
    comps_inf_peps_df.sort_index(inplace = True)
    
    #getting only rows representing peptides containing editing sites
    comps_inf_peps_df = comps_inf_peps_df[comps_inf_peps_df['infered_editing_comb'] != ([],[])]
    
    #joining colmns from components_dataframe before real sites calculation
    comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['strand']))
    comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['start_nuc']))
    comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['end_nuc']))
    comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['prot_beginning']))
    comps_inf_peps_df = comps_inf_peps_df.join(pd.DataFrame(comps_additional_data['prot_end']))
    
    #calculate real coordinates (relative to original rna seq)
    comps_inf_peps_df['real_infered_sites_coor_base1'] = comps_inf_peps_df.apply(lambda x: get_infered_editings_as_list(x.infered_editing_comb, x.strand, x.start_nuc, x.end_nuc), axis = 1)

    comps_inf_peps_df = comps_inf_peps_df.reset_index().set_index(['index','peptide','in_frame_coordinates_base0'])
    
    return comps_inf_peps_df
       

"""
from comps dataframe containing informative peptides
get a dictionary where keys are components id and values are lists of detectable sites
"""
#det_sites_dict = count_edited_sites(comps_inf_peps_df)
def get_editing_sites(comps_inf_peps_df):
    
    detectable_sites_dict = {}
    detectable_sites_num = 0
    
    for comp, comp_df in comps_inf_peps_df.groupby(level=0):
        sites_per_comp = []
        for index, row in comp_df.iterrows():
            [sites_per_comp.append(site) for site in row['real_infered_sites_coor_base1']]
        
        comp_detectabale_sites = list(set(sites_per_comp))
        detectable_sites_dict.update({comp:comp_detectabale_sites})
        detectable_sites_num += len(comp_detectabale_sites)
    
#    print(str(detectable_sites_num) + ' editing sites are detectable by at least one informative peptide')
    
    return detectable_sites_dict, detectable_sites_num 

        
                
def sites_detectable_in_both_versions(final_peps_df):
            
            
            


# =============================================================================
# #only_edited_peps_df = get_all_edited_peps(final_peps_df)
# def get_all_edited_peps(final_peps_df):
#     
#     only_edited_peps_dfs_list = []
#     
#     for pep, pep_df in final_peps_df.groupby(level=0):
#         if True in [x!=[([], [])] for x in pep_df['editing_combis'].values]:
#             only_edited_peps_dfs_list.append(pep_df)
# 
#     only_edited_peps_df = pd.concat(only_edited_peps_dfs_list) 
#     return only_edited_peps_df
# 
# 
# 
# def get_detectable_combinations(only_edited_peps_df):
#     
#     tot_editing_combis = 0
#     detectable_combis = 0
#     detectable_sites = 0
#     
#     
#     for pep, pep_df in only_edited_peps_df.groupby(level=0):
#         for index, row in pep_df.iterrows():
#             if row['editing_combis'] != [([],[])]:
#                 tot_editing_combis += len(row['editing_combis'])
#                 if len(pep_df) == 1:
#                     if row['infered_editing_comb'] != ([],[]):
#                         detectable_combis += 1
#                         
#     return tot_editing_combis, detectable_combis
# 
# =============================================================================
