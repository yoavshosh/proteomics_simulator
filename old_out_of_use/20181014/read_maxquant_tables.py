import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

tabel_path = 'C:/Users/user/Desktop/shahar/'
simulation_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/'


"""
from maxquant column of protein groups, creat a list of all protein groups of child peptide
"""
def comps_string_to_list(row, substr_to_del):
    return [x.replace(substr_to_del,"") for x in row.split(";")]


"""
read a table in txt file to dataframe
"""
def read_peptides_tabel(tabel_path, tabel_name = 'peptides.txt', fasta_file_name = 'squ'):
    
    data = []
    with open(tabel_path + tabel_name, "r") as f:
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
    return df


def get_detected_sources(row, maxquant_df):
    if row['peptide'] in maxquant_df.index:
        return maxquant_df.loc[row['peptide'],'proteins_list']
    else:
        return '-'
    
def check_detected_peptides(row, maxquant_df):
    if row['peptides'] in maxquant_df.index:
        return True
    else:
        return False
    

def compare_maxquant_and_simulation_results(simulation_df, maxquant_df):
    
    simulation_df['detected'] = simulation_df.apply(lambda row: check_detected_peptides(row, maxquant_df), axis = 1)
#    simulation_df['max_quant_sources'] = simulation_df.apply(lambda row: get_detected_sources(row, maxquant_df), axis = 1)
#    simulation_df['detected_proteins'] = simulation_df.apply(lambda row: len(row.max_quant_sources), axis = 1)
    
    #printing all sites to file
    return simulation_df


def remove_fasta_descriptions(input_path, input_fasta):
    from Bio import SeqIO
    writer =  FastaWriter(open(input_path + 'no_description_' + input_fasta , 'w'), wrap=None)
    writer.write_header()
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        writer.write_record(SeqRecord(record.seq, id = record.id,description = ''))
    writer.write_footer()

    
    
    
    
    
    
    
