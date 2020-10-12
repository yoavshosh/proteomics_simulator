import sys
import numpy as np
import pandas as pd
from Bio import SeqIO


def read_editing_sites_wrt_coding_seqs(anovar_file):

    def split_gene_and_prot_names(row):
        names = row['gene_name'].split(';')
        row['ucsc_id'] = names[0]
        row['HGNC_protein_id'] = names[1]
        return row
         
    columns = ['gene_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand']
    data = []
    
    with open(anovar_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            fields = line.split("\t")
            for i in range(len(fields)):
                fields[i] = fields[i].replace('\n','')
#            fields.append(str(fields[4]+'_'+fields[7]+'_'+fields[6]))
#            fields.append(str(fields[0].split(';')[0])+'_'+str(fields[2]))
            data.append(fields)
            
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(split_gene_and_prot_names, axis = 1)
    df['genomic_key_base1'] = df.apply(lambda row: row['chromosome']+'_'+row['strand']+'_'+row['genomic_position_base1'], axis = 1)
    df['coding_key_base1'] = df.apply(lambda row: row['ucsc_id']+'_'+row['position_base1'],axis = 1)
    return df


def genomic_keys_for_sites(row, sites_df):
    
    """
    for each site represented by peptide, get the genomic key from sites_df
    creates a column with a list of genomic keys corresponding to the sites in the list of real_inferred_a2g_sites_coor_base1 in peps df
    genomic keys are: crhomosome_strand_site-position-wrt-to-chromosome
    """

    protein_id = row['seq_id']
    protein_sites_keys = [protein_id + '_' + str(site) for site in row['real_inferred_a2g_sites_coor_base1']]
    try:
        genomic_sites_keys = [sites_df.loc[key, 'genomic_key_base1'] for key in protein_sites_keys]
        return genomic_sites_keys   
    except KeyError:
        print(protein_sites_keys)
        

def genomic_informative(row, groupd_genomic_keys_by_pep,all_peps):
    
    """
    operates on edited peptides dataframe
    
    the parameter - groupd_genomic_keys_by_pep - a sereis.
    for each peptide, the value is a nested list. each list within a list contain the genomic_keys of all sites represented by a particulae instance of a peptide
    the nested list will contain n lists where n is the number of peptide's instances
    
    if a peptide is not informative in the sense that there are identical peptides that comes from other versions of the same gene
    check if all lists within the nested list in groupd_genomic_keys_by_pep are identical - meaning if all instances of peptide will represent the smae site(s).
    if so - the peptide is genomic informative
    """
    def all_instances_edited(pep_instances):
        if len(pep_instances[pep_instances['editing_sites_inferred']==0]):
            return False
        else:
            return True

    def checkEqual1(iterator):
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == rest for rest in iterator)
    
    
    if all_instances_edited(all_peps[all_peps['peptide']==row['peptide']]):
        if not row['informative_peptide']:
            if checkEqual1([sorted(l) for l in groupd_genomic_keys_by_pep[row['peptide']]]):
                return True
            else:
                return False
        else:
            return True
    else:
        return False
        
        
    
if __name__ == '__main__':
    
#    peps_path = sys.argv[1]
#    peps_file = sys.argv[2]
#    sites_file = sys.argv[3]
    
    peps_path = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/A2G_only/results_from_all_variants_3mc_6minl_5500maxm_20maxes/'
    peps_file = 'peps_from_all_variants_fasta.pickle'
    sites_file = 'E:/RNA_editing_Large_files/human_editom/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor.txt'
    
    print('Reading sites table')
    sites_df = read_editing_sites_wrt_coding_seqs(sites_file)
    sites_df.set_index('coding_key_base1', inplace = True)
    all_sites_genomic_keys = set(list(sites_df[sites_df['mm_type']=='AG']['genomic_key_base1']))
    print('Total number of A2G recoding sites = ' + str(len(all_sites_genomic_keys)))
    
    print('Reading proteomics simulation results')
    peps_df = pd.read_pickle(peps_path + peps_file)
    print('Getting genomic A2G coordinates of sites per edited peptide')
    peps_df['genomic_a2g_sites_keys'] = peps_df.apply(lambda row: genomic_keys_for_sites(row,sites_df), axis = 1)
    
    groupd_genomic_keys_by_pep = peps_df[peps_df['editing_sites_inferred']>0].groupby('peptide')['genomic_a2g_sites_keys'].apply(list)
    
    print('Infer genomic-informative peptides')
    edited_peps = peps_df[peps_df['editing_sites_inferred']>0].copy()
    edited_peps['genomic_informative'] = edited_peps.apply(lambda row: genomic_informative(row,groupd_genomic_keys_by_pep,peps_df), axis = 1)
    
    print('Getting detectable sites')
    detectable_edited_peps = edited_peps[edited_peps['genomic_informative']]
    detectable_editing_sites = list(set([site for group in list(detectable_edited_peps['genomic_a2g_sites_keys']) for site in group]))
    print(str(len(detectable_editing_sites)) + ' sites are detectable in at least one version of a gene')
    print('Writing (genomic) informative edited peptides table')
    edited_peps.to_pickle(peps_path + 'edited_peptides.pickle')
    edited_peps.to_excel(peps_path + 'edited_peptides.xlsx')
    
