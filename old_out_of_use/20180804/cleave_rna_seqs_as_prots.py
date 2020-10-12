import re
import sys
import os
import pandas as pd
from peps_stats import *
#import logging
from edited_peptides_from_fasta_seqs import print_peptides_to_fasta, find_by_regex_in_header, infer_editing_comb, edited_peptides_from_rna_seq
#from edited_peptides_from_seq import create_edited_rna_peptides
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_rna


if __name__ == '__main__':

    #parameters for running the program
    input_path = sys.argv[1]
    input_fasta = sys.argv[2]
    strand_regex = re.compile('(?<=strand:\s)[^\s]+')
    start_nuc_regex = re.compile('(?<=start_nuc:\s)[^\s]+')
    end_nuc_regex = re.compile('(?<=end_nuc:\s)[^\s]+')
    prot_beginning_regex = re.compile('(?<=prot_beginning:\s)[^\s]+')
    prot_end_regex = re.compile('(?<=prot_end:\s)[^\s]+')
    a2g_header_regex = re.compile(r'(?<=a2g:\s).*?]')
    c2t_header_regex = re.compile(r'(?<=c2t:\s).*?]')
    digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
    missed_cleavages = 0
    look_ahead_for_editnig = 3
    look_behind_for_editing = 6
    min_length = 8*3
    max_mass = 4600 
    max_sites_per_pep = 20
    
    #creating ouput folder
    if not os.path.exists(input_path + 'results_from_' + input_fasta.split('.')[0] + '/'):
        os.makedirs(input_path + 'results_from_' + input_fasta.split('.')[0] + '/')
    output_path = input_path + 'results_from_' + input_fasta.split('.')[0] + '/'
    #string for various output files
    input_name = '_'.join(input_fasta.split('.'))
    
#    sys.stdout = open(input_path + input_fasta[:-6] + '_cleave_log.txt', 'w')
    
    final_peps_df, under_min_length_dict, over_max_sites_dict, over_mass_dict, over_mass_combs_dict, tot_editing_sites = edited_peptides_from_rna_seq(input_path, input_fasta, 
                                                                                                                                                      a2g_header_regex = a2g_header_regex, c2t_header_regex = c2t_header_regex,
                                                                                                                                                      digestion_rule = digestion_rule, missed_cleavages = missed_cleavages,
                                                                                                                                                      look_ahead_for_editnig = look_ahead_for_editnig, look_behind_for_editing = look_behind_for_editing,
                                                                                                                                                      min_length = min_length, max_mass = max_mass, max_sites_per_pep = max_sites_per_pep)

    #print dataframe to a fasta file
    print('writing final tabel to file')
    final_peps_df.to_pickle(output_path + 'peps_from_' + input_name + '.pickle')
    final_peps_df.to_excel(output_path + 'peps_from_' + input_name + '.xlsx') 

# =============================================================================
    #print all anattended peptides to different files
    print('writing short peptides coordinates')
    under_min_length_peps = open(output_path + 'under_min_length_peps_from_' + input_name + '.txt', "w")
    for k, v in under_min_length_dict.items():
        under_min_length_peps.write(k + ':' + str(v) + '\n')
    under_min_length_peps.close()
     
    print('writing coordinates of peptides with more than ' + str(max_sites_per_pep) + ' editing sites')
    over_max_sites_peps = open(output_path + 'over_max_sites_peps_from_' + input_name + '.txt', "w")
    for k, v in over_max_sites_dict.items():
        over_max_sites_peps.write(k + ':' + str(v) + '\n')
    over_max_sites_peps.close()
     
    print('writing coordinates of peptides heavier than ' + str(max_mass) + 'Da')
    over_weight_peps = open(output_path + 'overweight_peps_from_' + input_name + '.txt', "w")
    for k, v in over_mass_dict.items():
        over_weight_peps.write(k + ':' + str(v) + '\n')
    over_weight_peps.close()
     
    print('writing coordinates of peptides editing combinations which exceeds ' + str(max_mass) + 'Da, for peptides which some of their editing combinations are within mass range')
    over_weight_peps_combs = open(output_path + 'overweight_combinations_of_peps_from_' + input_name + '.txt', "w")
    for k, v in over_mass_combs_dict.items():
        over_weight_peps_combs.write(k + ':' + str(v) + '\n')
    over_weight_peps_combs.close()
# =============================================================================
    
    print('reading components raw data')
    comps_raw_df = get_components_df(input_path, input_fasta, strand_regex, start_nuc_regex, end_nuc_regex, prot_beginning_regex, prot_end_regex)
    print('creating detectable editing sites dataframe')
    detectable_editing_sitse_df = create_detectable_editing_sites_df(final_peps_df,comps_raw_df)
    print('writing detectable editing sites dataframe to file')
    detectable_editing_sitse_df.to_pickle(output_path + 'detectable_sites_from_' + input_name + '.pickle')
    detectable_editing_sitse_df.to_excel(output_path + 'detectable_sites_from_' + input_name + '.xlsx')
    
    print('getting all detectable sites')
    det_sits_dict, detectable_sites_num = get_editing_sites(detectable_editing_sitse_df)
    print(str(detectable_sites_num) + ' (out of ' + str(tot_editing_sites) + ') editing sites are detectable by at least one informative peptide')
    print('printing detectable sites list to file')
    det_sites_list_file = open(output_path + 'detectable_sites_from_' + input_name + '.txt', "w")
    for k, v in det_sits_dict.items():
        det_sites_list_file.write(k + ':' + str(v) + '\n')
    det_sites_list_file.close()    
    
    print('writing final tabel to fasta file')
    n, inf = print_peptides_to_fasta(final_peps_df,output_path,input_fasta)
    print('unique peptides generated and written to file: ' + str(n))
    print('informative unique peptides generated and written to file: ' + str(inf))
    