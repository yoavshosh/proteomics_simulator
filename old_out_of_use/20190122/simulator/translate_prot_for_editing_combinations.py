import sys
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from peps_stats import get_inferred_editings_as_list


def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'

def edit_rna_as_peptide(sequence, a_to_g_list, c_to_t_list):
    splited_ag = [sequence[i:j] for i,j in zip([0] + sorted([x+1 for x in a_to_g_list]), sorted([x for x in a_to_g_list]) + [None])]
    ag_swaped_seq = 'G'.join(splited_ag)
    splited_ct = [ag_swaped_seq[i:j] for i,j in zip([0] + sorted([x+1 for x in c_to_t_list]), sorted([x for x in c_to_t_list]) + [None])]
    new_seq = 'T'.join(splited_ct)
    return new_seq



def create_edited_proteins_all_represented_combinations(input_path, fasta_input, output_path, final_peps_df, max_edits_per_pep = None, allow_change_in_cleavage_sites = False):
    
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
        protein_basic_description = 'protein_start_nuc_base1:' + str(start_nuc) + '| protein_ends_nuc_base1:' + str(end_nuc) + '| strand: ' + strand
        
        for comb in combs_list[0]:
            
            #flag editing combination for print\dont print in proteins file
            edit_prot = True
            if max_edits_per_pep != None:
                if len([site for edit_type in comb for site in edit_type]) > max_edits_per_pep:
                    edit_prot = False
            if not allow_change_in_cleavage_sites and edit_prot:
                if final_peps_df.loc[index,'N_terminus'] != 'no_change' or final_peps_df.loc[index,'C_terminus'] != 'no_change' or final_peps_df.loc[index,'cancelled_cs_in_pep']:
                    edit_prot = False
            
            
            #editing proteins and writing to file if combination not already writen and combination do not exceed editing events
            if comb not in written_combs and edit_prot:
                if comb == ([], []): #the original sequence
                    comb_id = index + '_original'
                    protein = sequences_dict[index].translate()
                    protein_description = protein_basic_description
                else:
                    comb_id = index + '_edited_' + str(n)
                    protein_description = protein_basic_description + '| a2g_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list(comb,strand, start_nuc, end_nuc, editing_type = 0))
                    protein_description = protein_description + '| c2t_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list(comb,strand, start_nuc, end_nuc, editing_type = 1))
                    protein = Seq(edit_rna_as_peptide(str(sequences_dict[index]),comb[0],comb[1]), generic_dna).translate()
                    n+=1
                written_combs.append(comb)
                writer.write_record(SeqRecord(protein, id = comb_id, description = protein_description))
    
    writer.write_footer()


def create_fully_edited_proteins_fasta(input_path, fasta_input,output_path):
    
    strand_regex = re.compile('(?<=strand:\s)[^\s]+')
    prot_start_nuc_regex = re.compile('(?<=prot_start_nuc:\s)[^\s]+')
    prot_end_nuc_regex = re.compile('(?<=prot_end_nuc:\s)[^\s]+')
    prot_start_regex = re.compile('(?<=prot_start:\s)[^\s]+')
    prot_end_regex = re.compile('(?<=prot_end:\s)[^\s]+')
    original_orf_start_regex = re.compile('(?<=original_orf_start:\s)[^\s]+')
    original_orf_end_regex = re.compile('(?<=original_orf_end:\s)[^\s]+')
    a2g_header_regex = re.compile(r'(?<=a2g_base0:\s).*?]')
    c2t_header_regex = re.compile(r'(?<=c2t_base0:\s).*?]')
    
    
    writer = FastaWriter(open(output_path + 'fully_edited_and_native_proteins_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        
        descrp = record.description
        a2g = eval(find_by_regex_in_header(descrp,a2g_header_regex))
        c2t = eval(find_by_regex_in_header(descrp,c2t_header_regex))
        
        protein_basic_description = ''
        #translate native protein
        seq_id = record.id + '_original'
        protein = record.seq.translate()
        writer.write_record(SeqRecord(protein, id = seq_id, description = protein_basic_description))
    
        if len(a2g+c2t):
            seq_id = record.id + '_fully_edited'
            protein_description = protein_basic_description + '| a2g_editing_events_relative_to_coding_sequence_base0: ' + str(a2g)
            protein_description = protein_description + '| c2t_editing_events_relative_to_coding_sequence_base0: ' + str(c2t)
            edited_seq = Seq(edit_rna_as_peptide(str(record.seq),a2g,c2t), generic_dna)
            protein = edited_seq.translate()
            writer.write_record(SeqRecord(protein, id = seq_id, description = protein_description))
            if len(edited_seq)%3:
                print(record.id)
                print(len(record.seq))
                print(len(edited_seq))
        
    writer.write_footer()



if __name__ == '__main__':
    
#    input_path = sys.argv[1]
#    fasta_input = sys.argv[2]
#    df_path = sys.argv[3]
#    df_file = sys.argv[4]
#    max_edits_per_pep = eval(sys.argv[5])
#    allow_cs_changes = eval(sys.argv[6])
#    fully_edited = eval(sys.argv[7])
    
    input_path = 'E:/RNA_editing_Large_files/proteomics_simulation/'
    fasta_input = 'in_frame_rna_rec_only_from_shahar_squ_orfs.fasta'
    
#    input_path = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/'
#    fasta_input = 'in_frame_rna_from_human_coding_sequences.fasta'
    df_path = input_path + 'results_from_in_frame_rna_rec_only_from_shahar_squ_orfs_3mc_7minl_4600maxm_20maxes/'
    df_file = 'peps_from_in_frame_rna_rec_only_from_shahar_squ_orfs_fasta.pickle'
    max_edits_per_pep = None
    allow_cs_changes = False
    fully_edited = True
    
    print('Creating proteins fasta from sequences in ' + fasta_input)

    if fully_edited:
        create_fully_edited_proteins_fasta(input_path, fasta_input, input_path)
    else:
        print('Reading peptide DataFrame from ' + df_file)
        final_peps_df = pd.read_pickle(df_path + df_file)
        create_edited_proteins_all_represented_combinations(input_path, fasta_input, df_path, final_peps_df, max_edits_per_pep = max_edits_per_pep, allow_change_in_cleavage_sites = allow_cs_changes)
    
    print("Translation of peptides is finished")
    
    
# =============================================================================
# def create_fully_edited_proteins_fasta_old(input_path, fasta_input, output_path, final_peps_df, max_edits_per_pep = None, allow_change_in_cleavage_sites = True):
#     
#     #create a seq-id:sequence dictionary from input fasta file
#     sequences_dict = {}
#     for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
#         sequences_dict.update({record.id:record.seq})
#     
#     writer =  FastaWriter(open(output_path + 'fully_edited_and_native_proteins_considering_params_from_' + fasta_input , 'w'), wrap=None)
#     writer.write_header()
#     
#     final_peps_df_no_dup = final_peps_df.drop_duplicates(subset = 'seq_id', keep = 'first').copy()
#     final_peps_df_no_dup.set_index('seq_id', inplace = True)
#     final_peps_df.set_index('seq_id', inplace = True)
#     
#     if not allow_change_in_cleavage_sites:
#         final_peps_df = final_peps_df[final_peps_df['N_terminus']=='no_change']
#         final_peps_df = final_peps_df[final_peps_df['C_terminus']=='no_change']
#     
#     #for each seq_id, iterate over all editing combinations and creat edited peptides
#     for comp in list(sequences_dict.keys()):
#         try:
#             start_nuc = str(final_peps_df_no_dup.loc[comp,'prot_start_nuc_base1'])
#             end_nuc = str(final_peps_df_no_dup.loc[comp,'prot_end_nuc_base1'])
#             strand = str(final_peps_df_no_dup.loc[comp,'strand'])
#             protein_basic_description = 'protein_start_nuc_base1:' + str(start_nuc) + '| protein_start_nuc_base1:' + str(end_nuc) + '| strand: ' + strand
#             
#             #translate native protein
#             seq_id = comp + '_original'
#             protein = sequences_dict[comp].translate()
#             writer.write_record(SeqRecord(protein, id = seq_id, description = protein_basic_description))
#     
#             protein_mini_df = final_peps_df.loc[[comp],:]
#             editing_sites_for_prot = []
#                         
#             for index, row in protein_mini_df.iterrows():
#                 #flag editing combination for print\dont print in proteins file
#                 consider_editing_sites = True
#                 if max_edits_per_pep != None:
#                     if row['sites_in_peptide_range'] > max_edits_per_pep:
#                         consider_editing_sites = False
#                 
#                 if consider_editing_sites:
#                     editing_sites_for_prot.append(row['inferred_editing_combination_relative_to_sense_orf_base0'])
#                      
#             a2g_sites = sorted(list(set([site for comb in editing_sites_for_prot for site in comb[0]])))
#             c2t_sites = sorted(list(set([site for comb in editing_sites_for_prot for site in comb[1]])))
#             
# #            print(comp)
# #            print(a2g_sites)
#             #write_fullt_edited_prot
#             if len(a2g_sites+c2t_sites):
#                 seq_id = comp + '_fully_edited'
#                 protein_description = protein_basic_description + '| a2g_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list((a2g_sites,c2t_sites),strand, start_nuc, end_nuc, editing_type = 0))
#                 protein_description = protein_description + '| c2t_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list((a2g_sites,c2t_sites),strand, start_nuc, end_nuc, editing_type = 1))
#                 protein = Seq(edit_rna_as_peptide(str(sequences_dict[index]),a2g_sites,c2t_sites), generic_dna).translate()
#                 writer.write_record(SeqRecord(protein, id = seq_id, description = protein_description))
#         
#         except KeyError:
#             print(comp + ' has no valid peptides - not added to protein file')
#         except AttributeError:
#             print(comp + ' do not fit code')
#         
#     writer.write_footer()
# 
# =============================================================================
