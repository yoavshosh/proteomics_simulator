import re
import sys
import os
import pandas as pd
from peps_stats import *
#import logging
from edited_peptides_from_fasta_seqs import *
#from edited_peptides_from_seq import create_edited_rna_peptides
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_rna


if __name__ == '__main__':

    #parameters for running the program

#    input_fasta = 'in_frame_rna_rec_only_from_shahar_squ_orfs.fasta'
    input_fasta = 'in_frame_rna_rec_only_from_orfs_squ_sim.fa'
    input_path = 'C:/Users/shosh/Google_Drive/RNA_Editing/files/proteomics_simulation/'
    missed_cleavages = 0
    min_aa_number = 7
    max_mass = 4600 
    max_sites_per_pep = 20
    print_peps_fasta = False
#    input_path = sys.argv[1]
#    input_fasta = sys.argv[2]
#    missed_cleavages = int(sys.argv[3])
#    min_aa_number = int(sys.argv[4])
#    max_mass = int(sys.argv[5])
#    max_sites_per_pep = int(sys.argv[6])
#    print_peps_fasta = eval(sys.argv[7])
    
    strand_regex = re.compile('(?<=strand:\s)[^\s]+')
    prot_start_nuc_regex = re.compile('(?<=prot_start_nuc:\s)[^\s]+')
    prot_end_nuc_regex = re.compile('(?<=prot_end_nuc:\s)[^\s]+')
    prot_start_regex = re.compile('(?<=prot_start:\s)[^\s]+')
    prot_end_regex = re.compile('(?<=prot_end:\s)[^\s]+')
    original_orf_start_regex = re.compile('(?<=original_orf_start:\s)[^\s]+')
    original_orf_end_regex = re.compile('(?<=original_orf_end:\s)[^\s]+')
    a2g_header_regex = re.compile(r'(?<=a2g_base0:\s).*?]')
    c2t_header_regex = re.compile(r'(?<=c2t_base0:\s).*?]')
    digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
    
    look_ahead_for_editnig = 3
    look_behind_for_editing = 6
    min_length = min_aa_number*3

    
    param_str = '_'+str(missed_cleavages)+'mc_'+str(min_aa_number)+'minl_'+str(max_mass)+'maxm_'+str(max_sites_per_pep)+'maxes'
    #creating ouput folder
    if not os.path.exists(input_path + 'results_from_' + input_fasta.split('.')[0] + param_str+  '/'):
        os.makedirs(input_path + 'results_from_' + input_fasta.split('.')[0] + param_str+ '/')
    output_path = input_path + 'results_from_' + input_fasta.split('.')[0] + param_str + '/'
    #string for various output files
    input_name = '_'.join(input_fasta.split('.'))  #just a string for output files names
    
#    sys.stdout = open(input_path + input_fasta[:-6] + '_cleave_log.txt', 'w')
    
    print('General data for sequences:')
    peps_df, under_min_length_dict, over_max_sites_dict, over_mass_dict, over_mass_combs_dict, tot_editing_sites = edited_peptides_from_rna_seq(input_path, input_fasta, 
                                                                                                                                                a2g_header_regex = a2g_header_regex, c2t_header_regex = c2t_header_regex,
                                                                                                                                                digestion_rule = digestion_rule, missed_cleavages = missed_cleavages,
                                                                                                                                                look_ahead_for_editnig = look_ahead_for_editnig, look_behind_for_editing = look_behind_for_editing,
                                                                                                                                                min_length = min_length, max_mass = max_mass, max_sites_per_pep = max_sites_per_pep)

    #print dataframe to a fasta file
    print('Calculating additional peptides data')
    comps_df = get_components_df(input_path, input_fasta, strand_regex, prot_start_nuc_regex, prot_end_nuc_regex, prot_start_regex, prot_end_regex,original_orf_start_regex,original_orf_end_regex)
    final_peps_df = process_peptides_df(peps_df, comps_df)
    print('Writing final tabel to files (pickle + xls)')
    final_peps_df.to_pickle(output_path + 'peps_from_' + input_name + '.pickle')
    final_peps_df.to_excel(output_path + 'peps_from_' + input_name + '.xlsx', index = False) 

# =============================================================================
    #print all anattended peptides to different files
    print('Writing coordinates of peptides shorter than ' + str(min_aa_number) + ' AA')
    under_min_length_peps = open(output_path + 'under_min_length_peps_from_' + input_name + '.txt', "w")
    for k, v in under_min_length_dict.items():
        under_min_length_peps.write(k + ':' + str(v) + '\n')
    under_min_length_peps.close()
     
    print('Writing coordinates of peptides with more than ' + str(max_sites_per_pep) + ' editing sites (peptides not in final tabel)')
    over_max_sites_peps = open(output_path + 'over_max_sites_peps_from_' + input_name + '.txt', "w")
    for k, v in over_max_sites_dict.items():
        over_max_sites_peps.write(k + ':' + str(v) + '\n')
    over_max_sites_peps.close()
     
    print('Writing coordinates of peptides heavier than ' + str(max_mass) + 'Da')
    over_weight_peps = open(output_path + 'overweight_peps_from_' + input_name + '.txt', "w")
    for k, v in over_mass_dict.items():
        over_weight_peps.write(k + ':' + str(v) + '\n')
    over_weight_peps.close()
     
    print('Writing coordinates of peptides heavier than ' + str(max_mass) + 'Da, for peptides which some of their editing combinations are within mass range')
    over_weight_peps_combs = open(output_path + 'overweight_combinations_of_peps_from_' + input_name + '.txt', "w")
    for k, v in over_mass_combs_dict.items():
        over_weight_peps_combs.write(k + ':' + str(v) + '\n')
    over_weight_peps_combs.close()
# =============================================================================
    
    #print all detectable sites to file
    print('Filtering for detectable sites')
    detectable_editing_sitse_df = final_peps_df[np.logical_and(final_peps_df['edited'],final_peps_df['editing_sites_inferred']>0)]
    det_sits_dict, detectable_sites_num = get_editing_sites(detectable_editing_sitse_df)
    print('Writing detectable sites lists to file')
    det_sites_list_file = open(output_path + 'detectable_sites_from_' + input_name + '.txt', "w")
    for k, v in det_sits_dict.items():
        det_sites_list_file.write(k + ' | a2g detectable sites:' + str(v['a2g_detectable_sites']) + ' | c2t detectable sites:' + str(v['c2t_detectable_sites']) + '\n')
    det_sites_list_file.close()  
    
    print('Writing edited proteins to fasta file')
    create_fully_edited_proteins_fasta(input_path, input_fasta, output_path, final_peps_df)
    
    if print_peps_fasta:
        print('Writing final peptides tabel to fasta file')
        print_peptides_to_fasta(final_peps_df,output_path,input_fasta)
    
    print('Unique peptides generated: ' + str(len(final_peps_df)))
    print('Informative unique peptides: ' + str(len(final_peps_df[final_peps_df['informative_peptide']])))
    print(str(detectable_sites_num) + ' (out of ' + str(tot_editing_sites) + ') editing sites are detectable by at least one informative peptide (as a part of at least one editing combination)')
    print('Finished  :-)')