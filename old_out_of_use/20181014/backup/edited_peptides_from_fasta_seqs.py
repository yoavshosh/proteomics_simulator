import re
import sys
import timeit
import pandas as pd
from edited_peptides_from_seq import create_edited_rna_peptides, create_all_cleavage_sites
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna, IUPAC, generic_dna
from peps_stats import get_inferred_editings_as_list

#sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/' +'test3.txt', 'w')

#a2g_header_regex = re.compile(r'(?<=a2g:\s).*?]')
#c2t_header_regex = re.compile(r'(?<=c2t:\s).*?]')
#digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
#missed_cleavages = 1
#look_ahead_for_editnig = 3
#look_behind_for_editing = 6
#min_length = None
#max_length = None
#max_sites_per_pep = None
#input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
#input_fasta = 'sim.fasta'
#input_fasta = 'in_frame_rna_from_orfs_squ.fasta'

def edited_peptides_from_rna_seq(input_path,input_fasta, a2g_header_regex ,c2t_header_regex,
                                 digestion_rule,missed_cleavages,look_ahead_for_editnig = 3 ,look_behind_for_editing = 6,
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
        unique_peptides = 0
        over_mass_combs = 0
        n+=1
        
        #get editing sites lists from record header
        a2g_sites = eval(find_by_regex_in_header(record.description,a2g_header_regex))
        c2t_sites = eval(find_by_regex_in_header(record.description,c2t_header_regex))
        tot_editing_sites += len(a2g_sites) + len(c2t_sites)

        #create all peptides from all editing combinations
        rna_peptides_list, under_min_length, over_max_sites_peps, over_mass_peps, over_mass_combs_per_seq_dict = create_edited_rna_peptides(str(record.seq),a2g_sites,c2t_sites,digestion_rule, missed_cleavages = missed_cleavages,
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
            pep_data_list = [(str(Seq(pep.seq, generic_rna).translate()).replace('I','X').replace('L','X'), str(record.id), '_'+str(pep.coo[0])+'_'+str(pep.coo[1]), (pep.a2g,pep.c2t), pep.N_terminus, pep.C_terminus) for pep in rna_peptides_list]
            del(rna_peptides_list)
            #initial dataframe
            peps_df_temp = pd.DataFrame(pep_data_list,columns = ['peptide','seq_id','in_frame_coordinates_base0','editing_comb', 'N_terminus', 'C_terminus'])
            #groupby multiindex and create a list of all editing combinations for grouped rows
            peps_df = peps_df_temp.groupby(['peptide','seq_id','in_frame_coordinates_base0']).apply(lambda x: x['editing_comb'].tolist())
            #set name of new column in grouped df
            peps_df = pd.DataFrame({'editing_combinations_relative_to_sense_orf_base0': peps_df.values}, index = peps_df.index)
            #create another column of editing combinations intersection - sites that are obviosly edited
            peps_df['inferred_editing_combination_relative_to_sense_orf_base0'] = peps_df['editing_combinations_relative_to_sense_orf_base0'].apply(infer_editing_comb)
            #retrive N\C terminus change from pepd_df_temp
            peps_df_temp = peps_df_temp.loc[:,('peptide','seq_id','in_frame_coordinates_base0','N_terminus','C_terminus')]
            peps_df_temp.drop_duplicates(['peptide','seq_id','in_frame_coordinates_base0'],inplace = True)
            peps_df = peps_df.reset_index()
            peps_df = pd.merge(peps_df,peps_df_temp, on = ['peptide','seq_id','in_frame_coordinates_base0'])
#            peps_df.set_index(['peptide','seq_id','in_frame_coordinates_base0'],inplace = True)
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
excepts a final_peps_df in which peptides and print to fasta file
"""
def print_peptides_to_fasta(final_peps_df,o_path,input_fasta, write_infotmative_peptides_fasta = False):
    
    final_peps_df = final_peps_df.reset_index()
    final_peps_df.set_index(['peptide','seq_id','in_frame_coordinates_base0'], inplace = True)
    
    n=0
    inf = 0
    all_writer =  FastaWriter(open(o_path + 'all_peps_from_' + input_fasta , 'w'), wrap=None) #writer for all peptides
    all_writer.write_header()
    
    if write_infotmative_peptides_fasta:
        inf_writer =  FastaWriter(open(o_path + 'inf_peps_from_' + input_fasta , 'w'), wrap=None) #writer for informative peptides only
        inf_writer.write_header()
    
    #break dataframe to individuals peptides dataframes and iterate to get the peptide description
    for peptide, peptide_df in final_peps_df.groupby(level=0):
        n+=1
        description = ''
        peptide_df = peptide_df.reset_index(level = 0, drop = True)
        
        #break peptide df to individual source sequences dfs and iterate
        for seq_id, seq_for_pep_df in peptide_df.groupby(level=0):
            description += ' +++ seq_id:'+seq_id+' '
            seq_for_pep_df = seq_for_pep_df.reset_index(level = 0, drop = True)
            
            #ource sequences df to individual coordinate in source dfs and iterate
            for coo, row in seq_for_pep_df.iterrows():
                description+= '| coor:' + coo + ' editing_combs:'+str(row['editing_combinations_relative_to_sense_orf_base0']) + ' infered_editing_comb:'+str(row['inferred_editing_combination_relative_to_sense_orf_base0'])+' '
                
        # create a SeqRecord and write to all-peptides file
        record = SeqRecord(Seq(peptide, IUPAC.protein), id = 'pep_'+str(n), description = description)
        all_writer.write_record(record)
        
        #if peptide is generated from one sequence and from one location within that sequence - then its informative
        if len(peptide_df.index.to_series()) == 1:
            inf+=1
            if write_infotmative_peptides_fasta:
                inf_writer.write_record(record)
    
    all_writer.write_footer()        
    if write_infotmative_peptides_fasta:
        inf_writer.write_footer()         
    
    return n, inf    

    
"""
find regex (regex - pre compiled regex) in header
"""
def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return []
    

"""
find intersection of all editing combinations represented by peptide
"""
def infer_editing_comb(editing_combis):
    return (list(set.intersection(*map(set,[x[0] for x in editing_combis]))),
            list(set.intersection(*map(set,[x[1] for x in editing_combis]))))
    
    
def determine_edge(coor, seq_len):
    edge = None
    if '_0_' in coor:
        edge = 5
    elif '_'+str(seq_len) in coor:
        edge = 3
    return edge
    

def edit_rna_as_peptide(sequence, a_to_g_list, c_to_t_list):
    splited_ag = [sequence[i:j] for i,j in zip([0] + sorted([x+1 for x in a_to_g_list]), sorted([x for x in a_to_g_list]) + [None])]
    ag_swaped_seq = 'G'.join(splited_ag)
    splited_ct = [ag_swaped_seq[i:j] for i,j in zip([0] + sorted([x+1 for x in c_to_t_list]), sorted([x for x in c_to_t_list]) + [None])]
    new_seq = 'T'.join(splited_ct)
    return new_seq

def create_fully_edited_proteins_fasta(input_path, fasta_input, output_path, final_peps_df):
    
    #create a seq-id:sequence dictionary from input fasta file
    sequences_dict = {}
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        sequences_dict.update({record.id:record.seq})
    
    writer =  FastaWriter(open(output_path + 'proteins_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    #creating a dataframe of all editing cominations per protein
    comps_editing_combs = final_peps_df.groupby('seq_id').agg({'editing_combinations_relative_to_sense_orf_base0':lambda x: sorted([comb for sublist in list(x) for comb in sublist])})
    
    #for each seq_id, iterate over all editing combinations and creat edited peptides
    final_peps_df = final_peps_df.drop_duplicates(subset = 'seq_id', keep = 'first') #removing duplicates as only data in seq_id level is now needed 
    final_peps_df.set_index('seq_id', inplace = True)
    for index, combs_list in comps_editing_combs.iterrows():
        written_combs = []
        n=1
        start_nuc = final_peps_df.loc[index,'prot_start_nuc_base1']
        end_nuc = final_peps_df.loc[index,'prot_end_nuc_base1']
        strand = final_peps_df.loc[index,'strand']
        protein_basic_description = 'protein_start_nuc_base1:' + str(start_nuc) + '| protein_start_nuc_base1:' + str(end_nuc) + '| strand: ' + strand
        
        for comb in combs_list[0]:
            if comb not in written_combs:
                if comb == ([], []): #the original sequence
                    comb_id = index + '_original'
                    protein = sequences_dict[index].translate()
                    protein_description = protein_basic_description
                else:
                    comb_id = index + '_edited_' + str(n)
                    protein_description = protein_basic_description + '| a2t_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list(comb,strand, start_nuc, end_nuc, editing_type = 0))
                    protein_description = protein_description + '| c2t_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list(comb,strand, start_nuc, end_nuc, editing_type = 1))
                    protein = Seq(edit_rna_as_peptide(str(sequences_dict[index]),comb[0],comb[1]), generic_dna).translate()
                written_combs.append(comb)
                writer.write_record(SeqRecord(protein, id = comb_id, description = protein_description))
                n+=1
    
    writer.write_footer()
