import os
import sys
import subprocess
import multiprocessing as mp
import pandas as pd
from analize_maxquant_results import *


def read_samples_df(maxquant_files_list):
    
    columns = ['sample_name', 'path']
    data = []
    with open(maxquant_files_list, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            line_arr = line.split('\t')
            for i in range(len(line_arr)):
                line_arr[i] = line_arr[i].replace('\n','')
            data.append(line_arr)
    
    df = pd.DataFrame(data = data, columns = columns)
    return df
    


def search_edited_peptides(maxquant_file, proteomics_simulation_df, pos):   
    mq_df = read_peptides_tabel(maxquant_file)    
    discovered_peptides = proteomics_simulation_df.merge(mq_df, left_on  = 'peptide',right_on = 'isobaric_peptide',how = 'inner',suffixes = ('_comparison','_simulation'))
    return pos,len(mq_df),discovered_peptides

    
def func(pos,maxquant_file):
    return (pos,len(maxquant_file))

if __name__ == '__main__':
#    
#    maxquant_files_list = 'E:/RNA_editing_Large_files/human_editom/maxquant_files/maxquant_files_list_mycomp.txt'
#    edited_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/edited_peptides.pickle'
#    all_peptides_from_simulation = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/peps_from_all_variants_c2t_a2g_fasta.pickle'
    
    maxquant_files_list = sys.argv[1]
    edited_peptides_from_simulation = sys.argv[2]
    all_peptides_from_simulation = sys.argv[2]
#    
    #path to output file
    l = maxquant_files_list.split('/')
    results_file = '/'.join(l[:-1]) + '/results_from_' + l[-1]
    peptide_file = '/'.join(l[:-1]) + '/peptides_lists_from_' + l[-1]
    
    
    maxquant_files_df = read_samples_df(maxquant_files_list)
    n_files = len(maxquant_files_df)
    print('\nReading simulation tables')
    edited_peptides_df = pd.read_pickle(edited_peptides_from_simulation)
    all_peptides_df = pd.read_pickle(all_peptides_from_simulation)

    #start parallel processing - maxquant results tables bs the simulation peptides table
    print('\nGetting detection results for edited peptides')
    pool_edited = mp.Pool(n_files)
    results_edited = [pool_edited.apply_async(search_edited_peptides, (maxquant_files_df.loc[i,'path'], edited_peptides_df, i)) for i in range(n_files)]
    pool_edited.close()
    pool_edited.join()
    output_edited = [r.get() for r in results_edited]
    print('\nSorting results')
    output_edited = sorted(output_edited, key=lambda x: x[0])
    
    print('\nGetting detection results for all peptides')
    pool_all = mp.Pool(n_files)
    results_all_peps = [pool_all.apply_async(search_edited_peptides, (maxquant_files_df.loc[i,'path'], all_peptides_df, i)) for i in range(n_files)]
    pool_all.close()
    pool_all.join()
    output_all_peps = [r.get() for r in results_edited]
    print('\nSorting results')
    output_all_peps = sorted(output_all_peps, key=lambda x: x[0])
    
    

    maxquant_files_df['total_peptides'] = maxquant_files_df.apply(lambda row: output_edited[row.name][1], axis = 1)
    maxquant_files_df['edited_peptides'] = maxquant_files_df.apply(lambda row: len(output_edited[row.name][2]), axis = 1)
    maxquant_files_df['a2g_detected'] = maxquant_files_df.apply(lambda row: list(set([site for group in list(output_edited[row.name][2]['genomic_a2g_sites_keys']) for site in group])), axis = 1)
    maxquant_files_df['c2t_detected'] = maxquant_files_df.apply(lambda row: list(set([site for group in list(output_edited[row.name][2]['genomic_c2t_sites_keys']) for site in group])), axis = 1)

    print('\nWriting General Analyses results')
    maxquant_files_df.to_csv(results_file, sep='\t', index=None)
    
    for i, row in maxquant_files_df.iterrows():
        print(row['sample_name'] + ': ' + str(len(row['a2g_detected'])) + ' A2G sites,' + str(len(row['c2t_detected'])) + ' C2T sites')
    
    print('\nWriting peptides combonations')
    with open(peptide_file, "w") as pep_file:
        pep_file.write('sample_name\tseq_id\tpeptide\A2G_editing sites\tC2T_Editing_sites\tsites_in_peptide_range\tsites_in_peptide_range\toriginal_version_detected\n')
        for out in output_edited:
            for i, row in out[2].iterrows():
                all_detected_peps = output_all_peps[out[0]][2]
                other_versions = all_detected_peps.loc[np.logical_and(all_detected_peps['seq_id']==row['seq_id'],all_detected_peps['in_frame_coordinates_base0']==row['in_frame_coordinates_base0'])]
                if len(other_versions[other_versions['edited']==False]):
                    original_version_found = True
                else:
                    original_version_found = False
                pep_file.write(maxquant_files_df.loc[out[0], 'sample_name'] + '\t' + row['seq_id'] + '\t' + row['peptide'] + '\t' + str(row['genomic_a2g_sites_keys']) + '\t' + str(row['genomic_c2t_sites_keys']) + '\t' + str(row['sites_in_peptide_range']) + '\t' + str(original_version_found) + '\n')
    