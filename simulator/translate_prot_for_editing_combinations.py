import sys
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from edited_peptides_from_seq import edit_rna_as_peptide, all_mm


def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'


def create_proteins_for_each_peptide(input_path, fasta_input, output_path, final_peptides, allow_change_in_cleavage_sites = False):
    """
    for each sequence create the native protein
    and create a version of thath protein for each peptide
    """
    
    final_edited_peptides = final_peptides[final_peptides['edited']]
    
    #create a seq-id:sequence dictionary from input fasta file
    sequences_dict = {}
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        sequences_dict.update({record.id:record.seq})    
    
    writer =  FastaWriter(open(output_path + 'proteins_per_peptide_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
        
    for key, mrna_sequence, in sequences_dict.items():
        
        #first print the native protein
        comb_id = key + '|original'
        protein = mrna_sequence.translate()
        writer.write_record(SeqRecord(protein, id = comb_id, description = ''))
        
        edited_peptides = final_edited_peptides[final_edited_peptides['seq_id'] == key]
        
        n=1
        for index, row in edited_peptides.iterrows():
                
            #flag editing combination for print\dont print in proteins file
            edit_prot = True
            if not allow_change_in_cleavage_sites and edit_prot:
                if final_peps_df.loc[index,'N_terminus'] != 'no_change' or final_peps_df.loc[index,'C_terminus'] != 'no_change' or final_peps_df.loc[index,'cancelled_cs_in_pep']:
                    edit_prot = False
                
            if edit_prot:
                permutation_coor = tuple(int(x) for x in row['permutation_coor_base0'].split('_') if x != '')
                protein = mrna_sequence[:permutation_coor[0]].translate() + row['biological_extended_peptide'] + mrna_sequence[permutation_coor[1]+1:]
                comb_id = key + '|edited_'+str(n) + '\t' + str(row['editing_combinations_relative_to_coding_seq_base0'])
                writer.write_record(SeqRecord(protein, id = comb_id, description = ''))
                n+=1
    
    writer.write_footer()
        
                
def create_edited_proteins_all_represented_combinations(input_path, fasta_input, output_path, final_peps_df, max_edits_per_pep = None, allow_change_in_cleavage_sites = False):
    
    """
    for each sequence create the native protein
    and create a version of that protein for each editing combination represented by that each edited peptide
    """
    
    #create a seq-id:sequence dictionary from input fasta file
    sequences_dict = {}
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        sequences_dict.update({record.id:record.seq})
    
    writer =  FastaWriter(open(output_path + 'proteins_per_combination_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    #creating a dataframe of all editing cominations per protein
#    comps_editing_combs = final_peps_df.groupby('seq_id').agg({'editing_combinations_relative_to_sense_orf_base0':lambda x: sorted([comb for sublist in list(x) for comb in sublist])})
    comps_editing_combs = final_peps_df.groupby('seq_id')['editing_combinations_relative_to_coding_seq_base0'].aggregate(lambda x: list(x))
    #for each seq_id, iterate over all editing combinations and creat edited peptides
    final_peps_df = final_peps_df.drop_duplicates(subset = 'seq_id', keep = 'first') #removing duplicates as only data in seq_id level is now needed 
    final_peps_df.set_index('seq_id', inplace = True)
    
    for index, combs_nested_list in comps_editing_combs.iteritems():
    
        written_combs = []
        n=1
        protein_basic_description = ''
        length = len(sequences_dict[index])
        flattened_comb_list = [c for l in combs_nested_list for c in l]
        
        for comb in flattened_comb_list:
            
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
                if comb == ([],[],[],[],[],[],[],[],[],[],[],[]): #the original sequence
                    comb_id = index + '_original'
                    protein = sequences_dict[index].translate()
                    protein_description = protein_basic_description
                else:
                    comb_id = index + '_edited_' + str(n)
                    protein_description = protein_basic_description + '| editing_combinations_base0_wrt_to_coding_sequence: ' + str(comb)
                    protein = Seq(edit_rna_as_peptide(str(sequences_dict[index]),(0,length-1),comb), generic_dna).translate()
                    n+=1
                written_combs.append(comb)
                writer.write_record(SeqRecord(protein, id = comb_id, description = protein_description))
    
    writer.write_footer()


def create_fully_edited_proteins_fasta(input_path, fasta_input,output_path):  
    """
    for each sequence create a native protein version and a fully edited version
    """    

    mm_headers = {}
    [mm_headers.update({mm:re.compile(r'(?<='+mm+'_base0:\s).*?]')}) for mm in all_mm]
    
    writer = FastaWriter(open(output_path + 'fully_edited_and_native_proteins_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        
        sites_dict = {}
        [sites_dict.update({mm:sorted(eval(find_by_regex_in_header(record.description,mm_headers[mm])))}) for mm in all_mm]
        sites_number = sum([len(sites_dict[mm]) for mm in all_mm])
        length = len(record.seq)
        comb = tuple([sites_dict[mm] for mm in all_mm])
        
        protein_basic_description = ''
        #translate native protein
        seq_id = record.id + '_original'
        protein = record.seq.translate()
        writer.write_record(SeqRecord(protein, id = seq_id, description = protein_basic_description))
    
        if sites_number:
            seq_id = record.id + '_fully_edited'
            protein_description = protein_basic_description + '| editing_combinations_base0_wrt_to_coding_sequence: ' + str(comb)
            edited_seq = Seq(edit_rna_as_peptide(str(record.seq),(0,length-1),comb), generic_dna)
            protein = edited_seq.translate()
            writer.write_record(SeqRecord(protein, id = seq_id, description = protein_description))
            if len(edited_seq)%3:
                print(record.id)
                print(len(record.seq))
                print(len(edited_seq))
        
    writer.write_footer()
    
        


def create_peptides_fasta(input_path,fasta_input,peps_df, extention=15):
    
    writer =  FastaWriter(open(input_path + 'peptides_extanded_by'+str(extention)+'_from' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        
        prot = record.seq.translate()
        for i,row in peps_df[peps_df['seq_id']==record.id].iterrows():
            rna_pep_coor = row['in_frame_coordinates_base0'].split('_')
            pep_start = int(rna_pep_coor[1])/3
            pep_end = int(rna_pep_coor[2])/3
            seq_start = max(0,pep_start-extention)
            seq_end = min(pep_end+extention,len(prot))
            extented_pep = prot[seq_start:pep_start] + row['biological_peptide'] + prot[min(pep_end+1,len(prot)):seq_end]
            if not row['edited']:
                seq_id = record.id + '_original_' + str(seq_start*3)+'_'+str(seq_end*3) + '_pep_id_' + str(i)
            else:
                seq_id = record.id + '_' + str(seq_start*3) + '_' + str(seq_end*3) + '_editing_range' + row['permutation_coor_base0'] + '_pep_id_' + str(i)
            writer.write_record(SeqRecord(extented_pep, id = seq_id, description = ''))
    
    writer.write_footer()
        
        
    

if __name__ == '__main__':
    
    input_path = '/private7/projects/yoav/in_silico_proteomics/squ/new_fixed_transcriptome/'
    fasta_input = 'mrna_coding_sequences_from_orfs_squ.fa'
    peps_df_path = '/private7/projects/yoav/in_silico_proteomics/squ/new_fixed_transcriptome/proteomics_simulation/results_from_ag_sites_no_stops_for_regular_inf_2mc_7minl_4700maxm_20maxes/peps_from_ag_sites_no_stops_for_regular_inf.txt'
    
    print('Reading peptides dataframe')
    peps_df = pd.read_csv(peps_df_path, sep = '\t', index_col=False)
    print('Creating extended peptides fasta')
    create_peptides_fasta(input_path,fasta_input,peps_df, extention=15)
    print('Finished')
    
#    input_path = sys.argv[1]
#    fasta_input = sys.argv[2]
#    df_path = sys.argv[3]
#    df_file = sys.argv[4]
#    max_edits_per_pep = eval(sys.argv[5])
#    allow_cs_changes = eval(sys.argv[6])
#    fully_edited = eval(sys.argv[7])
    
#    input_path = 'E:/RNA_editing_Large_files/proteomics_simulation/'
#    fasta_input = 'in_frame_rna_rec_only_from_shahar_squ_orfs.fasta'
    
#    input_path = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/'
#    fasta_input = 'in_frame_rna_from_human_coding_sequences.fasta'
#    df_path = input_path + 'results_from_in_frame_rna_rec_only_from_shahar_squ_orfs_3mc_7minl_4600maxm_20maxes/'
#    df_file = 'peps_from_in_frame_rna_rec_only_from_shahar_squ_orfs_fasta.pickle'
#    max_edits_per_pep = None
#    allow_cs_changes = False
#    fully_edited = True
    
#    print('Creating proteins fasta from sequences in ' + fasta_input)
#
#    if fully_edited:
#        create_fully_edited_proteins_fasta(input_path, fasta_input, input_path)
#    else:
#        print('Reading peptide DataFrame from ' + df_file)
#        final_peps_df = pd.read_pickle(df_path + df_file)
#        create_edited_proteins_all_represented_combinations(input_path, fasta_input, df_path, final_peps_df, max_edits_per_pep = max_edits_per_pep, allow_change_in_cleavage_sites = allow_cs_changes)
#    
#    print("Translation of peptides is finished")
    
    
    
    
  
