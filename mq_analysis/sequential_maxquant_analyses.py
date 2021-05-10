import sys
import re
import datetime
import argparse
import multiprocessing as mp
import pandas as pd
import numpy as np
import traceback
from Bio.SeqUtils.ProtParam import ProteinAnalysis

all_mm = ['AG','AC','AT','CA','CG','CT','GA','GC','GT','TA','TG','TC']


def read_piptides_from_csv_create_df_with_correct_types(file):
    
#    columns = ['seq_id','in_frame_coordinates_base0','peptide','editing_combinations_relative_to_coding_seq_base0',
#               'editing_combinations_intersection_base0','biological_peptide','biological_extended_peptide','permutation_coor_base0',
#               'N_terminus','C_terminus','sites_in_permutation_range','cancelled_cs_in_pep','editing_combinations_intersection_base1',
#               'peptide_count','informative_peptide','edited','informative_original_version','editing_sites_inferred',
#               'ucsc_id','HGNC_protein_id','genomic_keys','genomic_informative']
    data = []
    with open(file, "r") as f:
        content = f.readlines()
        columns = content[0].split("\t")
        for i in range(len(columns)):
            columns[i] = columns[i].replace('\n','').replace('\r','')
        
        for j, line in enumerate(content[1:]):
            fields = line.split("\t")
            for i in range(len(fields)):
                fields[i] = fields[i].replace('\n','').replace('\r','')
            fields[3] = eval(fields[3])
            fields[4] = eval(fields[4])
            fields[10] = int(fields[10])
            fields[11] = eval(fields[11])
            fields[12] = eval(fields[12])
            fields[13] = int(fields[13])
            fields[14] = eval(fields[14])
            fields[15] = eval(fields[15])
            try:  #could be either True or False, but aslo a string "no_original_version"
                fields[16] = eval(fields[16])  
            except NameError: #
                pass
            fields[17] = int(fields[17])
            fields[20] = eval(fields[20])
            fields[21] = eval(fields[21])
            data.append(fields)
    
    f.close()
    return pd.DataFrame(data = data, columns = columns)
            
def read_editing_sites_wrt_coding_seqs(anovar_file):
    """
    read dataframe from anovar output
    read file with 4 columns as in columns and translate return dataframe
    """
    def split_gene_and_prot_names(row,col_name):
        names = row[col_name].split(';')
        row['ucsc_id'] = names[0]
        row['HGNC_protein_id'] = names[1]
        return row

    columns = ['gene_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand', 'aa_change','editing_level']
    data = []
    
    with open(anovar_file, "r") as f:
        content = f.readlines()
        for j, line in enumerate(content):
            fields = line.split("\t")
            for i in range(len(fields)):
                fields[i] = fields[i].replace('\n','').replace('\r','')
            data.append(fields)
    
    f.close()
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(lambda row: split_gene_and_prot_names(row, 'gene_name'), axis = 1)
    df['genomic_key_base1'] = df.apply(lambda row: row['HGNC_protein_id'] + ';' + row['chromosome']+'_'+row['strand']+'_'+row['genomic_position_base1'], axis = 1)
    df['coding_key_base1'] = df.apply(lambda row: row['ucsc_id']+'_'+row['position_base1'],axis = 1)
    return df


def read_peptides_tabel(tabel_path ,sequence_col = 'Sequence'): 
    """
    read a table in txt file to dataframe and add some columns
    """
    print('\nCreating peptides dataframe from ' + tabel_path)
    
    def comps_string_to_list(proteins_string, substr_to_del='',spliter=';',spliter_expected_per_name=1):
        """
        from maxquant column of protein groups, creat a list of all protein groups of child peptide
        """
        regex = re.compile(spliter)
        indeces_to_split = [-1]
        for i,g in enumerate(regex.finditer(proteins_string)):
            if (i+1)%(spliter_expected_per_name+1)==0:
                indeces_to_split.append(g.start())
        parts = [proteins_string[i+1:j] for i,j in zip(indeces_to_split, indeces_to_split[1:]+[None])]
        parts = [p.replace(substr_to_del,'') for p in parts]
        return parts
    
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
    df['proteins_list'] = df.apply(lambda row: comps_string_to_list(row['Proteins']), axis = 1)
    df['protein_sources'] = df.apply(lambda row: len(row.proteins_list), axis = 1)
    df['isobaric_peptide'] = df.apply(lambda row: row[sequence_col].replace('L','X').replace('I','X'), axis = 1)
    df['mq_length'] = df.apply(lambda row: len(row[sequence_col]), axis = 1)
    df['mq_molecular_weight'] = df.apply(lambda row: ProteinAnalysis(row[sequence_col]).molecular_weight(), axis = 1)
    return df


def read_samples_df(maxquant_files_list, sub_path_for_mq_dirs):
    
    def add_peptides_table_path_field(row, sub_path):
        return row['raw_path']+sub_path
    
    columns = ['sample_name', 'sub_name', 'quantification', 'phosphoproteom','tissue', 'raw_path', 'organized_sample', 'notes']
    data = []
    with open(maxquant_files_list, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            line_arr = line.split('\t')
            for i in range(len(line_arr)):
                line_arr[i] = line_arr[i].replace('\n','')
            data.append(line_arr)
    
    df = pd.DataFrame(data = data, columns = columns)
    df['mq_peptides_path'] = df.apply(lambda row: add_peptides_table_path_field(row, sub_path_for_mq_dirs), axis=1)
    return df


def search_edited_peptides(maxquant_file, edited_peptides_and_native_versions_df, sites_df, pos):   

    try:
        #read peptides discovered by maxquant
        mq_df = read_peptides_tabel(maxquant_file)    
        
        #merge dataframe to get all discovered native and edited peptides
        discovered_edited_peptides_and_native_versions = edited_peptides_and_native_versions_df.merge(mq_df, left_on  = 'peptide',right_on = 'isobaric_peptide',how = 'inner',suffixes = ('_comparison','_simulation'))
    
#            discovered_edited_peptides_and_native_versions['permutation_start_nuc'] = discovered_edited_peptides_and_native_versions.apply(lambda row: int(row['in_frame_coordinates_base0'].split('_')[1]), axis = 1)
#            discovered_edited_peptides_and_native_versions['permutation_end_nuc'] = discovered_edited_peptides_and_native_versions.apply(lambda row: int(row['in_frame_coordinates_base0'].split('_')[2]), axis = 1)
    
        if len(discovered_edited_peptides_and_native_versions):    
            #cleane_sites_df for sites in genes with proteomics evidence
            discovered_genes = list(set(list(discovered_edited_peptides_and_native_versions['HGNC_protein_id'])))
            sites_df = sites_df[sites_df['HGNC_protein_id'].isin(discovered_genes)]
            sites_edited_and_native_evidence = check_sites_discovery(sites_df, discovered_edited_peptides_and_native_versions)
            
            print(str(pos) + ': ' + maxquant_file + ' ' + str(len(mq_df)))    
            return pos,len(mq_df),sites_edited_and_native_evidence
        else:
            print(str(pos) + ': ' + maxquant_file + '    NO RELEVANT PEPTIDES')
            return pos,len(mq_df),None
        
    
    except Exception as e:
        print(str(pos) + ': ' + maxquant_file + '    PROBLEM USING FILE')
        print(e)
        return -1
    
    
def count_edited_peptides(discovered_peptides):
    if discovered_peptides is None:
        return 0
    else:
        edited_peptides = discovered_peptides[discovered_peptides['edited']]
        return len(set(list(edited_peptides['peptide'])))
    


def check_sites_discovery(sites_df, discovered_peps_df):
    
    """
    !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!
    this function creats duplicates in cases where peptide support more than one site.
    find a way to rm dups and store data per each supported sites (like aa_change) in a list that corresponds to editing sites list represented by peptide
    Note that this could affect calculations in the plot_maxquant_results_for_multiple_runs.py script !!!!
    !!!!!!!!!!!!!!
    !!!!!!!!!!!!!!
    !!!!!!!!!!!!!!
    """
    
    
    def collect_edited_unedited_peptides_in_site(peps_from_variant, genomic_site ,coding_position_base0):
        """
        this sub function returns a dataframe of all peptides containing a passed site in its edited and unedited versions
        """
        edited_peptides = []
        unedited_peptides = []
        
        for index, row in peps_from_variant.iterrows():
            
            #check if peptide supports edited version edited versions    
            edited_site = False
            for i, mm_type in enumerate(row['genomic_keys']):
                for j, s in enumerate(mm_type):
                    if s == genomic_site:
                        edited_site = True
                        temp_df = peps_from_variant.loc[[index]]
                        temp_df['site_is_edited'] = temp_df.apply(lambda x: edited_site, axis = 1)
                        edited_peptides.append(temp_df)
                        break
                if edited_site:
                    break
            
            #if peptide does not support edited version, check if supports unedited version                    
            if not edited_site: 
                unedited_site = True
                for comb in row['editing_combinations_relative_to_coding_seq_base0']:
                    for i, mm_type in enumerate(comb):
                        for j, s in enumerate(mm_type):
                            if s == coding_position_base0:
                                unedited_site = False
                if unedited_site:
                    temp_df = peps_from_variant.loc[[index]]
                    temp_df['site_is_edited'] = temp_df.apply(lambda x: edited_site, axis = 1)
                    unedited_peptides.append(temp_df)
        
        if len(edited_peptides) + len(unedited_peptides):
            return pd.concat(edited_peptides + unedited_peptides)
        else:
            return None
    
    all_supporting_peptides = []
    for s_index, s_row in sites_df.iterrows():
        variant = s_row['ucsc_id']
        coding_position_base0 = int(s_row['position_base0'])
        genomic_site = s_row['genomic_key_base1']
        site_mm_type = s_row['mm_type']
        site_aa_change = s_row['aa_change']
        site_editing_level = s_row['editing_level']
        
        #all discovered (genomic informative) peptides from variant that contain the potential editing
        peps_from_variant = discovered_peps_df[discovered_peps_df['ucsc_id']==variant]
        peps_from_variant = peps_from_variant[np.logical_and(peps_from_variant['permutation_start_nuc'] <= coding_position_base0, peps_from_variant['permutation_end_nuc'] >= coding_position_base0)]
                
        supporting_peptides = collect_edited_unedited_peptides_in_site(peps_from_variant, genomic_site ,coding_position_base0)              
        
        if supporting_peptides is not None:
            supporting_peptides['position_base0'] = supporting_peptides.apply(lambda x: coding_position_base0, axis = 1)
            supporting_peptides['genomic_site'] = supporting_peptides.apply(lambda x: genomic_site, axis = 1)
            supporting_peptides['mm_type'] = supporting_peptides.apply(lambda x: site_mm_type, axis = 1)
            supporting_peptides['aa_change'] = supporting_peptides.apply(lambda x: site_aa_change, axis = 1)
            supporting_peptides['editing_level'] = supporting_peptides.apply(lambda x: site_editing_level, axis = 1)
            all_supporting_peptides.append(supporting_peptides)
    
    if len(all_supporting_peptides):    
        return pd.concat(all_supporting_peptides)
    else:
        return None

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Analysis of MaxQuant proteomics results wrt peptides from proteomics simulation')
    run_parser = parser.add_argument_group('Run Analysis')
    run_parser.add_argument('-samples_list', dest='sampels_list', action='store', required = True, help='Path to file with list of samples and their data (fields: samples_name, sample_sub_name, quantification_method, tissue, phosphoproteom(True/False) path_to_raw_files_dir, notes1, notes2')
    run_parser.add_argument('-peptides', dest='peptides', action='store', required = True, help='path to edited peptides and their native versions from proteomics simulator')
    run_parser.add_argument('-sites', dest='sites', action='store', required = True, help='path to editing sites list used in proteomics simulation')
    run_parser.add_argument('-combined_output_dir_name', dest='combined_output_dir_name', action='store', default = 'combined', help='in raw files paths in samples list - the sub path in which MQ peptides results resides')
    run_parser.add_argument('-threads', dest='n_workers', action='store', default = '30', help='number of chiled processes for mq results analysis')
    run_parser.add_argument('-o', dest='output_suffix', action='store', default = datetime.datetime.now().strftime("%Y%m%d%H%M%S"), help='output suffix for output files')
    arguments = parser.parse_args()
    
    maxquant_files_list = arguments.sampels_list
    edited_peptides_and_native_versions_file = arguments.peptides
    all_editing_sites_file = arguments.sites
    combined_output_dir_name = arguments.combined_output_dir_name
    n_workers = int(arguments.n_workers)
    output_suffix = arguments.output_suffix

#    maxquant_files_list = 'E:/RNA_editing_Large_files/human_editom/maxquant_files/maxquant_files_list_mycomp.txt'
#    edited_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/edited_peptides.pickle'
#    all_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/peps_from_all_variants_c2t_a2g_fasta.pickle'
#    edited_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_sim_2mc_6minl_4600maxm_20maxes/edited_peptides.pickle'
#    all_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_sim_2mc_6minl_4600maxm_20maxes/peps_from_sim_fa.pickle'
#    maxquant_files_list = 'E:/RNA_editing_large_files/proteomics_simulator_squeed/proteomics_simulation/mq_files_info.txt'
#    edited_peptides_and_native_versions_file = 'E:/RNA_editing_Large_files/proteomics_simulator_squeed/proteomics_simulation/results_from_ag_sites_2mc_7minl_4700maxm_20maxes__/edited_peptides_and_native_versions.pickle'
#    all_editing_sites_file = 'E:/RNA_editing_Large_files/proteomics_simulator_squeed/squ_editing_site_info_new_format.txt'

    #path to output files
    l = maxquant_files_list.split('/')
    results_file = '/'.join(l[:-1]) + '/results_from_' + l[-1]+'_'+output_suffix
    sites_detection_file = '/'.join(l[:-1]) + '/sites_detection_from_' + l[-1] +'_'+output_suffix
    
    print('Reading Editing Sites table')
    sites_df = read_editing_sites_wrt_coding_seqs(all_editing_sites_file)
    
    print('\nReading samples table')
    maxquant_files_df = read_samples_df(maxquant_files_list, combined_output_dir_name+'/txt/peptides.txt')
    
    print('\nReading edited peptides table')
    if edited_peptides_and_native_versions_file[-7:] == '.pickle':
        edited_peptides_and_native_versions_df = pd.read_pickle(edited_peptides_and_native_versions_file, compression=None)
    elif edited_peptides_and_native_versions_file[-4:] == '.txt':
        edited_peptides_and_native_versions_df = read_piptides_from_csv_create_df_with_correct_types(edited_peptides_and_native_versions_file)
    else:
        raise Exception('peptides file should be txt or pickle format. file is ' + edited_peptides_and_native_versions_file)
        
    edited_peptides_and_native_versions_df = edited_peptides_and_native_versions_df[edited_peptides_and_native_versions_df['genomic_informative']]
    edited_peptides_and_native_versions_df['permutation_start_nuc'] = edited_peptides_and_native_versions_df.apply(lambda row: int(row['permutation_coor_base0'].split('_')[1]), axis = 1)
    edited_peptides_and_native_versions_df['permutation_end_nuc'] = edited_peptides_and_native_versions_df.apply(lambda row: int(row['permutation_coor_base0'].split('_')[2]), axis = 1)
    
    #start parallel processing - maxquant results tables bs the simulation peptides table
    print('\nGetting detection results for edited peptides (and native versions)')
    pool_edited = mp.Pool(n_workers)
    results_edited = [pool_edited.apply_async(search_edited_peptides, (maxquant_files_df.loc[i,'mq_peptides_path'], edited_peptides_and_native_versions_df.copy(), sites_df.copy(), i)) for i in range(len(maxquant_files_df))]
    pool_edited.close()
    pool_edited.join()
    output_edited_native = [r.get() for r in results_edited]
    n_problematic_files = sum([1 for i in range(len(output_edited_native)) if output_edited_native[i]==-1])
    n_good_files = len(output_edited_native)-n_problematic_files
    
    #adding general results to maxquant_fils_df
    maxquant_files_df['total_peptides'] = maxquant_files_df.apply(lambda row: -1 if output_edited_native[row.name] == -1 else output_edited_native[row.name][1], axis = 1)
    maxquant_files_df = maxquant_files_df[maxquant_files_df['total_peptides'] != -1]
    maxquant_files_df['edited_peptides'] = maxquant_files_df.apply(lambda row: count_edited_peptides(output_edited_native[row.name][2]), axis = 1)
    for i,mm in enumerate(all_mm):
        maxquant_files_df[mm+'_detected'] = maxquant_files_df.apply(lambda row: [] if output_edited_native[row.name][2] is None else list(set([site for group in list(output_edited_native[row.name][2]['genomic_keys']) for site in group[i]])), axis = 1)
    maxquant_files_df['total_editing_sitse'] = maxquant_files_df.apply(lambda row: sum([len(row[mm+'_detected']) for mm in all_mm]), axis = 1)
    
    print('\n\n' + str(n_good_files) + ' Good files')
    print(str(n_problematic_files) + ' Problematic files\n')
    print('\nWriting General Analyses results')
    maxquant_files_df.to_csv(results_file, sep='\t', index=None)
    
    #calculating discovered editing sites per mm type
    all_discovered_sites = [[],[],[],[],[],[],[],[],[],[],[],[]]
    for i, row in maxquant_files_df.iterrows():
        if row['total_peptides'] > 0:
            sites_detected_n = ''
            for t,mm in enumerate(all_mm):
                sites_detected_n += str(len(row[mm+'_detected'])) + mm + ' '
                all_discovered_sites[t] += row[mm+'_detected']
            print(row['sample_name']+'_'+row['sub_name'] + ': ' + sites_detected_n)

    print('Total Sites Discovered:')
    all_sites_detected_n = ''
    for t,mm in enumerate(all_mm):
        all_sites_detected_n+=str(len(set(all_discovered_sites[t]))) + mm + ' '
    print(all_sites_detected_n)
    print('Total: ' + str(sum([len(set(all_discovered_sites[t])) for t in range(len(all_mm))])))
       
    print('\nWriting sites nave and edited versions detection')
    all_detections_dfs_list = []
    for out in output_edited_native:
        if out != -1:
            pos = out[0]
            temp_df = out[2]
            if temp_df is not None:
                for col in ['sample_name', 'sub_name', 'quantification', 'tissue', 'mq_peptides_path']:
                    temp_df[col] = temp_df.apply(lambda x: maxquant_files_df.loc[pos,col], axis = 1)
                all_detections_dfs_list.append(temp_df)
    
    sites_results_df = pd.concat(all_detections_dfs_list, sort = False)
    sites_results_df.to_csv(sites_detection_file, sep='\t', index=None)
    
    exit(0)
    
