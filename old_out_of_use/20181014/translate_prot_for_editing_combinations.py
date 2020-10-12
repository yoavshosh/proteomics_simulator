import re
import sys
import timeit
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from peps_stats import get_inferred_editings_as_list

def edit_rna_as_peptide(sequence, a_to_g_list, c_to_t_list):
    splited_ag = [sequence[i:j] for i,j in zip([0] + sorted([x+1 for x in a_to_g_list]), sorted([x for x in a_to_g_list]) + [None])]
    ag_swaped_seq = 'G'.join(splited_ag)
    splited_ct = [ag_swaped_seq[i:j] for i,j in zip([0] + sorted([x+1 for x in c_to_t_list]), sorted([x for x in c_to_t_list]) + [None])]
    new_seq = 'T'.join(splited_ct)
    return new_seq

def create_fully_edited_proteins_fasta(input_path, fasta_input, output_path, final_peps_df, max_edits_per_pep = None):
    
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
            
            #flag editing combination for print\dont print in proteins file
            edit_prot = True
            if max_edits_per_pep != None:
                if len([site for edit_type in comb for site in edit_type]) > max_edits_per_pep:
                    edit_prot = False
            
            #editing proteins and writing to file if combination not already writen and combination do not exceed editing events
            if comb not in written_combs and edit_prot:
                if comb == ([], []): #the original sequence
                    comb_id = index + '_original'
                    protein = sequences_dict[index].translate()
                    protein_description = protein_basic_description
                else:
                    comb_id = index + '_edited_' + str(n)
                    protein_description = protein_basic_description + '| a2t_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list(comb,strand, start_nuc, end_nuc, editing_type = 0))
                    protein_description = protein_description + '| c2t_editing_events_relative_to_original_sequence_base1: ' + str(get_inferred_editings_as_list(comb,strand, start_nuc, end_nuc, editing_type = 1))
                    protein = Seq(edit_rna_as_peptide(str(sequences_dict[index]),comb[0],comb[1]), generic_dna).translate()
                    n+=1
                written_combs.append(comb)
                writer.write_record(SeqRecord(protein, id = comb_id, description = protein_description))
    
    writer.write_footer()

if __name__ == '__main__':
    
    input_path = sys.argv[1]
    fasta_input = sys.argv[2]
    df_path = sys.argv[3]
    df_file = sys.argv[4]
    max_edits_per_pep = eval(sys.argv[5])
    
    print('Reading peptide DataFrame from ' + df_file)
    final_peps_df = pd.read_pickle(df_path + df_file)
    print('Creating proteins fasta from sequences in ' + fasta_input)
    create_fully_edited_proteins_fasta(input_path, fasta_input, input_path, final_peps_df, max_edits_per_pep = max_edits_per_pep)
    
    
    
    
    