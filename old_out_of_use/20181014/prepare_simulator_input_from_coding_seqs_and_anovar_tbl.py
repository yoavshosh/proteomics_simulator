import re
import os
import sys
import pandas as pd
import numbers
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis


"""
find regex (regex - pre compiled regex) in string
"""
def find_by_regex_in_str(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'


"""
read dataframe from anovar output
"""
def read_anovar_data(input_path, input_file, only_recoding = True):
     
    anovar_gene_change_regex = re.compile('(.+?)\:(.+)\:(.+).?\:c.(.+?)\:p.(.+)')
    
    columns = ['gene_name', 'mm_type', 'mrna_mm_position','nuc_mm','prot_mm_position','aa_mm','genomic_coordinates']
    data = []
    
    with open(input_path + input_file, "r") as f:
        content = f.readlines()
        for line in content:
            fields = line.split("\t")
            mm_type = fields[1]
            gene_name_and_changes = fields[2].split(',')
            genomic_coordinates = fields[3]+'_'+fields[4]+'_'+fields[8]
            for mm in gene_name_and_changes:
                if mm != '':
                    mm_data = find_by_regex_in_str(mm,anovar_gene_change_regex)
                    gene_name = mm_data[1]
                    mrna_mm_position = int(mm_data[3][1:-1])
                    nuc_mm = mm_data[3][0]+mm_data[3][-1]
                    prot_mm_position = int(mm_data[4][1:-1])
                    aa_mm = mm_data[4][0]+mm_data[4][-1]
                    data.append([gene_name, mm_type, mrna_mm_position, nuc_mm, prot_mm_position, aa_mm, genomic_coordinates])
    
    anovar_df = pd.DataFrame(data = data, columns = columns)
    if only_recoding:
        anovar_df = anovar_df[anovar_df['mm_type'] == 'nonsynonymous SNV']
    
    
    anovar_df_a2g = anovar_df[anovar_df['nuc_mm'] == 'AG']
    anovar_df_a2g.set_index('gene_name', inplace = True)
    anovar_df_a2g = anovar_df_a2g.groupby(anovar_df_a2g.index).agg({list(anovar_df_a2g.columns)[1]:lambda x: list(set(list(x)))})

    anovar_df_c2t = anovar_df[anovar_df['nuc_mm'] == 'CT']
    anovar_df_c2t.set_index('gene_name', inplace = True)
    anovar_df_c2t = anovar_df_c2t.groupby(anovar_df_c2t.index).agg({list(anovar_df_c2t.columns)[1]:lambda x: list(set(list(x)))})
        
    return anovar_df, anovar_df_a2g, anovar_df_c2t


def create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs(fasta_path,fasta_file,out_path,anovar_df_a2g,anovar_df_c2t,only_recoding=True):
    
    if only_recoding:
        out_str = 'rec_only'
    else:
        out_str = 'all_sites'
    
    writer =  FastaWriter(open(out_path + 'in_frame_rna_' + out_str + '_from_' + fasta_file , 'w'), wrap=None)
    writer_bad = FastaWriter(open(out_path + 'bad_seq_from_' + fasta_file , 'w'), wrap=None)
    writer.write_header()
    writer_bad.write_header()
    for record in SeqIO.parse(open(fasta_path + fasta_file, "r"), "fasta"):
     
        final_sequence = record.seq[0:-3]
        try:
            a2g_list = [k-1 for k in anovar_df_a2g.loc[record.id][0]]
        except KeyError:
            a2g_list = []
        try:
            c2t_list = [k-1 for k in anovar_df_c2t.loc[record.id][0]]
        except KeyError:
            c2t_list = []        
        prot_start_nuc = 1
        prot_end_nuc = len(final_sequence) - 3
        prot_start = 'first_met_in_original_orf'
        prot_end = 'original_sense_strand_orf_end'
        strand = '+'
        orf_start = 1
        orf_end = len(record.seq) - 3
        description_str = '| a2g: ' + str(a2g_list) + ' | c2t: ' + str(c2t_list) + '| prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + ' | strand: ' + strand + ' | prot_start_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc) + ' | original_orf_start: ' + str(orf_start) + ' | original_orf_end: ' + str(orf_end)
    
        if (Seq(str(record.seq[0:3]), generic_dna).translate() != 'M' and Seq(str(record.seq[-3:len(record.seq)]), generic_dna).translate() != '*') or '*' in Seq(str(final_sequence), generic_dna).translate():
            writer_bad.write_record(record)
        else:
            current_record = SeqRecord(final_sequence, id = record.id, description = description_str)
            writer.write_record(current_record)
            
    writer_bad.write_footer()
    writer.write_footer()
    
    
if __name__ == '__main__':
        
#    input_path = 'C:/Users/user/Google_Drive/RNA_Editing/files/Yeast/results_from_yishai/'
#    fasta_file = 'SC_refGeneMrna.fa'
#    anovar_file = 'RES_final_result_editLvl_ann.exonic_variant_function'
#    only_recoding = True

    input_path = sys.argv[1]
    fasta_file = sys.argv[2]
    anovar_file = sys.argv[3]
    only_recoding = True  
    
    if not os.path.exists(input_path + 'proteomics_simulation/'):
        os.makedirs(input_path + 'proteomics_simulation/')
    out_path = input_path + 'proteomics_simulation/'
    
    print('Reading and processing Anovar tabel')
    anovar_df, anovar_df_a2g, anovar_df_c2t = read_anovar_data(input_path,anovar_file,only_recoding = only_recoding)
    print('Writing fasta file in the proteomics simulator format')
    create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs(input_path,fasta_file,out_path,anovar_df_a2g,anovar_df_c2t,only_recoding = only_recoding)
    
    