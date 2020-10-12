import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqUtils.ProtParam import ProteinAnalysis


"""
from maxquant column of protein groups, creat a list of all protein groups of child peptide
"""
def comps_string_to_list(row, substr_to_del):
    return [x.replace(substr_to_del,"") for x in row.split(";")]


"""
read a table in txt file to dataframe
"""
def read_peptides_tabel(tabel_path, table_name = 'peptides.txt', fasta_file_name = 'squ',sequence_col = 'Sequence'):
    
    print('\nCreating peptides dataframe from ' + tabel_path + table_name)
    
    data = []
    with open(tabel_path + table_name, "r") as f:
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
    writer =  FastaWriter(open(input_path + 'no_description_' + input_fasta , 'w'), wrap=None)
    writer.write_header()
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        writer.write_record(SeqRecord(record.seq, id = record.id,description = ''))
    writer.write_footer()

    
def read_simulation_table_from_pickle(simulation_path):
    print('\nCreating simulation table from ' + simulation_path)
    return pd.read_pickle(simulation_path)


def get_peptides_from_xls_tbl(xlsx_ms_tabel, peptides_col = 'peptide'):
    print('\nCreating ms/ms results table from ' + xlsx_ms_tabel)
    df = pd.read_excel(xlsx_ms_tabel)
    df['isobaric_peptide'] = df.apply(lambda x: x[peptides_col].replace('I','X').replace('L','X'), axis = 1)
    df['length'] = df.apply(lambda x: len(x[peptides_col]), axis = 1)
    df['molecular_weight'] = df.apply(lambda x: ProteinAnalysis(x[peptides_col]).molecular_weight(), axis = 1)
    
    return df
    

if __name__ == '__main__':
    
    mq_tables_path = 'E:/RNA_editing_Large_files/ms_analysis/native_and_fully_edited/'
    mq_table_names_list = ['mq_fdr005.txt','mq_fdr01.txt','mq_fdr1.txt','mq_fdr5.txt']
    simulation_path = 'E:/RNA_editing_Large_files/proteomics_simulation/results_from_in_frame_rna_rec_only_from_shahar_squ_orfs_3mc_7minl_4600maxm_20maxes/peps_from_in_frame_rna_rec_only_from_shahar_squ_orfs_fasta.pickle'   
    comparison_msms_tabel = 'E:/RNA_editing_Large_files/ms_analysis/native_and_fully_edited/tech_ms_no_dup.xlsx'
    min_length = 7
    max_weight = 4600

    comparison_msms_df = get_peptides_from_xls_tbl(comparison_msms_tabel)
    simulation_df = read_simulation_table_from_pickle(simulation_path)
    comparison_msms_processed_resutlts = comparison_msms_df.merge(simulation_df, left_on  = 'isobaric_peptide',right_on = 'peptide',how = 'left',suffixes = ('_comparison','_simulation'))
    comparison_msms_processed_resutlts.to_excel(comparison_msms_tabel.replace('.xlsx','_processed') + '.xlsx')

# =============================================================================
    
    #Analyzing comparison_msms table
    print('\nResult for comparison msms analysis:')
    unique_comparison_peptides = set(list(comparison_msms_processed_resutlts['peptide_comparison']))
    print('Unique peptides: ' + str(len(unique_comparison_peptides)))
    
    comparison_msms_processed_fit_simulation_inf = comparison_msms_processed_resutlts[comparison_msms_processed_resutlts[['seq_id']].notnull()['seq_id']]
    comparison_msms_processed_fit_simulation_inf = comparison_msms_processed_fit_simulation_inf[comparison_msms_processed_fit_simulation_inf['informative_peptide']]
    comparison_msms_processed_fit_simulation_inf = comparison_msms_processed_fit_simulation_inf[comparison_msms_processed_fit_simulation_inf['length']>=min_length]
    comparison_msms_processed_fit_simulation_inf = comparison_msms_processed_fit_simulation_inf[comparison_msms_processed_fit_simulation_inf['molecular_weight']<=max_weight]
    unique_comparison_inf_peptides_fit_simulation_params = list(comparison_msms_processed_fit_simulation_inf['peptide_comparison'])
    print(f'Informative peptides for simulation parameters (min_len={min_length}, max_mw={max_weight}): {len(unique_comparison_inf_peptides_fit_simulation_params)}')
     
    edited_unique_comparison_inf_peptides_fit_simulation_params = list(comparison_msms_processed_fit_simulation_inf[comparison_msms_processed_fit_simulation_inf['editing_sites_inferred']>0]['peptide_comparison'])
    print(f'Edited Peptides: {len(edited_unique_comparison_inf_peptides_fit_simulation_params)}')
    
    one_zero_site_no_cs_changes = comparison_msms_processed_fit_simulation_inf[np.logical_and(comparison_msms_processed_fit_simulation_inf['N_terminus']=='no_change',comparison_msms_processed_fit_simulation_inf['C_terminus']=='no_change')]
    one_zero_site_no_cs_changes = one_zero_site_no_cs_changes[np.logical_or(one_zero_site_no_cs_changes['sites_in_peptide_range']==1,one_zero_site_no_cs_changes['sites_in_peptide_range']==0)]
    edited_unique_comparison_inf_peptides_fit_simulation_params_one_site_no_cs_change = list(one_zero_site_no_cs_changes[one_zero_site_no_cs_changes['editing_sites_inferred']==1]['peptide_comparison'])
    print(f'Edited Peptides - No Cleavage Sites Change and only one site in range: {len(edited_unique_comparison_inf_peptides_fit_simulation_params_one_site_no_cs_change)}')
# =============================================================================

    for tbl in mq_table_names_list:
        mq_df = read_peptides_tabel(mq_tables_path,table_name = tbl)    
        mq_processed_resutlts = mq_df.loc[:,['Sequence','isobaric_peptide','proteins_list','protein_sources','mq_length','mq_molecular_weight']].merge(comparison_msms_df, left_on  = 'Sequence',right_on = 'peptide',how = 'left',suffixes = ('_mq','_comparison'))
        mq_processed_resutlts = mq_processed_resutlts.merge(simulation_df, left_on = 'isobaric_peptide_mq',right_on = 'peptide',how = 'left',suffixes = ('_msms','_simulation'))
        mq_processed_resutlts.to_excel(mq_tables_path + tbl.replace('txt','xlsx'))
        
# =============================================================================
        #Analyzing MaxQuant results fot tbl and compairing with comparison msms results
        
        print('\nResults for ' + str(tbl)+':')
        unique_mq_peptides = set(list(mq_processed_resutlts['Sequence']))
        unique_mq_comparison_intersection = [x for x in unique_mq_peptides if x in unique_comparison_peptides]
        print(f'Unique peptides: {len(unique_mq_peptides)}, {len(unique_mq_comparison_intersection)} intersct') 
     
        mq_processed_fit_simulation_inf = mq_processed_resutlts[mq_processed_resutlts[['seq_id']].notnull()['seq_id']]
        mq_processed_fit_simulation_inf = mq_processed_fit_simulation_inf[mq_processed_fit_simulation_inf['informative_peptide']]
        mq_processed_fit_simulation_inf = mq_processed_fit_simulation_inf[mq_processed_fit_simulation_inf['mq_length']>=min_length]
        mq_processed_fit_simulation_inf = mq_processed_fit_simulation_inf[mq_processed_fit_simulation_inf['mq_molecular_weight']<=max_weight]
        unique_mq_inf_peptides_fit_simulation_params = list(mq_processed_fit_simulation_inf['Sequence'])
        unique_mq_comparison_inf_peptides_fit_simulation_params = [x for x in unique_mq_inf_peptides_fit_simulation_params if x in unique_comparison_inf_peptides_fit_simulation_params]
        print(f'Informative peptides for simulation parameters (min_len={min_length}, max_mw={max_weight}): {len(unique_mq_inf_peptides_fit_simulation_params)}, {len(unique_mq_comparison_inf_peptides_fit_simulation_params)} intersect')
    
        edited_unique_mq_inf_peptides_fit_simulation_params = list(mq_processed_fit_simulation_inf[mq_processed_fit_simulation_inf['editing_sites_inferred']>0]['Sequence'])
        edited_unique_mq_comparison_inf_peptides_fit_simulation_params = [x for x in edited_unique_mq_inf_peptides_fit_simulation_params if x in edited_unique_comparison_inf_peptides_fit_simulation_params]
        print(f'Edited Peptides: {len(edited_unique_mq_inf_peptides_fit_simulation_params)}, {len(edited_unique_mq_comparison_inf_peptides_fit_simulation_params)} intersect')
        
        mq_one_zero_site_no_cs_changes = mq_processed_fit_simulation_inf[np.logical_and(mq_processed_fit_simulation_inf['N_terminus']=='no_change',mq_processed_fit_simulation_inf['C_terminus']=='no_change')]
        mq_one_zero_site_no_cs_changes = mq_one_zero_site_no_cs_changes[np.logical_or(mq_one_zero_site_no_cs_changes['sites_in_peptide_range']==1,mq_one_zero_site_no_cs_changes['sites_in_peptide_range']==0)]
        edited_unique_mq_inf_peptides_fit_simulation_params_one_site_no_cs_change = list(mq_one_zero_site_no_cs_changes[mq_one_zero_site_no_cs_changes['editing_sites_inferred']==1]['Sequence'])
        edited_unique_mq_comparison_inf_peptides_fit_simulation_params_one_site_no_cs_change = [x for x in edited_unique_mq_inf_peptides_fit_simulation_params_one_site_no_cs_change if x in edited_unique_comparison_inf_peptides_fit_simulation_params_one_site_no_cs_change]
        print(f'Edited Peptides - No Cleavage Sites Change and only one site in range: {len(edited_unique_mq_inf_peptides_fit_simulation_params_one_site_no_cs_change)}, {len(edited_unique_mq_comparison_inf_peptides_fit_simulation_params_one_site_no_cs_change)} intersect')
        
# =============================================================================
