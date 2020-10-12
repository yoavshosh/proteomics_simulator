import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqUtils.ProtParam import ProteinAnalysis



def comps_string_to_list(row, substr_to_del):
    """
    from maxquant column of protein groups, creat a list of all protein groups of child peptide
    """
    return [x.replace(substr_to_del,"") for x in row.split(";")]




def read_peptides_tabel(tabel_path, fasta_file_name = '' ,sequence_col = 'Sequence'):    
    """
    read a table in txt file to dataframe
    """
    
    print('\nCreating peptides dataframe from ' + tabel_path)
    
    data = []
    with open(tabel_path, "r") as f:
        content = f.readlines()
        columns = content[0].split('\t')
        for i,line in enumerate(content[1:]):
            line_arr = line.split('\t')
            data.append(line_arr)
            
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(pd.to_numeric, errors='ignore')
    df = df.replace(np.nan, '', regex=True)
    df['proteins_list'] = df.apply(lambda row: comps_string_to_list(row['Proteins'], fasta_file_name + '|'), axis = 1)
    df['protein_sources'] = df.apply(lambda row: len(row.proteins_list), axis = 1)
    df['isobaric_peptide'] = df.apply(lambda row: row[sequence_col].replace('L','X').replace('I','X'), axis = 1)
    df['mq_length'] = df.apply(lambda row: len(row[sequence_col]), axis = 1)
    df['mq_molecular_weight'] = df.apply(lambda row: ProteinAnalysis(row[sequence_col]).molecular_weight(), axis = 1)
    return df


def original_version_detected(row, discovered_peptides):
   
    seq_id = row['seq_id']
    coor = r
    
     


if __name__ == '__main__':
   
    mq_table = 'E:/RNA_editing_Large_files/human_editom/maxquant_files/PXD004930_brain_cancer_PF_all_variants_fdr5.txt'
    
    edited_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/edited_peptides.pickle'
    all_peptides = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/peps_from_all_variants_c2t_a2g_fasta.pickle'
      
#    mq_table = sys.argv[1]
#    edited_peptides_from_simulation = sys.argv[2] 
    
    simulation_df = pd.read_pickle(all_peptides)
    edited_peptides = pd.read_pickle(edited_peptides_from_simulation)
    mq_df = read_peptides_tabel(mq_table)
    
    discovered_peptides = simulation_df.merge(mq_df, left_on  = 'peptide',right_on = 'isobaric_peptide',how = 'inner',suffixes = ('_comparison','_simulation'))
    discovered_peptides.to_pickle(edited_peptides_from_simulation.replace('.pickle','_discovered_by_maxquant.pickle'))
    discovered_peptides.to_excel(edited_peptides_from_simulation.replace('.pickle','_discovered_by_maxquant.xlsx'))
    
    discovered_edited_peptides = edited_peptides.merge(mq_df, left_on  = 'peptide',right_on = 'isobaric_peptide',how = 'inner',suffixes = ('_comparison','_simulation'))
    
    discovered_editing_sites_df = discovered_edited_peptides[discovered_edited_peptides['genomic_informative']].copy()
    a2g_discovered_editing_sites_list = list(set([site for group in list(discovered_editing_sites_df['genomic_a2g_sites_keys']) for site in group]))
    c2t_discovered_editing_sites_list = list(set([site for group in list(discovered_editing_sites_df['genomic_c2t_sites_keys']) for site in group]))
    
    print(str(len(a2g_discovered_editing_sites_list)) + ' A2G Editing sites discovered')
    print(str(len(c2t_discovered_editing_sites_list)) + ' C2T Editing sites discovered')
    