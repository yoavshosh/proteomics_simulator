import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from edited_peptides_from_seq import all_mm


def read_refGene_df(path_to_refGene_file):
    """    
    string name;        	"Name of gene (usually transcript_id from GTF)"
    string chrom;       	"Chromosome name"
    char[1] strand;     	"+ or - for strand"
    uint txStart;       	"Transcription start position"
    uint txEnd;         	"Transcription end position"
    uint cdsStart;      	"Coding region start"
    uint cdsEnd;        	"Coding region end"
    uint exonCount;     	"Number of exons"
    uint[exonCount] exonStarts; "Exon start positions"
    uint[exonCount] exonEnds;   "Exon end positions"
    int score;            	"Score"
    string name2;       	"Alternate name (e.g. gene_id from GTF)"
    string cdsStartStat; 	"enum('none','unk','incmpl','cmpl')"
    string cdsEndStat;   	"enum('none','unk','incmpl','cmpl')"
    lstring exonFrames; 	"Exon frame offsets {0,1,2}"
    """
    
    columns = ['unknown_number', 'variant', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
               'exonCount', 'exonStarts', 'exonEnds', 'score', 'gene', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    data = []    
    
    with open(path_to_refGene_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            fields = line.split("\t")
            for i in range(len(fields)):
                fields[i] = fields[i].replace('\n','')
            data.append(fields)
    
    ref_Gene_df = pd.DataFrame(data = data, columns = columns)
    return ref_Gene_df


def split_gene_and_prot_names(row,col_name):
    names = row[col_name].split(';')
    row['ucsc_id'] = names[0]
    row['HGNC_protein_id'] = names[1]
    return row


def read_editing_sites_wrt_coding_seqs(anovar_file):

    columns = ['gene_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand', 'aa_change', 'editing_level']
    data = []
    
    with open(anovar_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            fields = line.split("\t")
            for i in range(len(fields)):
                fields[i] = fields[i].replace('\n','')
            data.append(fields)
            
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(lambda row: split_gene_and_prot_names(row, 'gene_name'), axis = 1)
    df['genomic_key_base1'] = df.apply(lambda row: row['HGNC_protein_id'] + ';' + row['chromosome']+'_'+row['strand']+'_'+row['genomic_position_base1'], axis = 1)
    df['coding_key_base1'] = df.apply(lambda row: row['ucsc_id']+'_'+row['position_base1'],axis = 1)
    return df


def genomic_keys_for_sites(row, sites_df):
    
    """
    for each site represented by peptide, get the genomic key from sites_df
    creates a column with a list of genomic keys corresponding to the sites in the list of real_inferred_a2g_sites_coor_base1 in peps df
    genomic keys are: crhomosome_strand_site-position-wrt-to-chromosome
    """
    
    genomic_sites_keys = [[],[],[],[],[],[],[],[],[],[],[],[]]
    protein_id = row['ucsc_id']
#    protein_id = row['seq_id']
    
    for i,mm in enumerate(all_mm):
        protein_sites_keys = [protein_id + '_' + str(site) for site in row['editing_combinations_intersection_base1'][i]]        
        if len(protein_sites_keys):
            try:
                genomic_sites_keys[i] += [sites_df.loc[key, 'genomic_key_base1'] for key in protein_sites_keys]
            except KeyError:
                raise Exception('Some site was not found in sites df. ' + protein_id + ' ' + str(protein_sites_keys))
    
    return tuple(genomic_sites_keys)


def genomic_informative(row, groupd_genomic_keys_per_mm_dict, all_peps):
    
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
        """
        check if all instances of peptide are edited or some are from an original sequences
        """
        if len(pep_instances[pep_instances['editing_sites_inferred']==0]):
            return False
        else:
            return True

    def checkEqual(iterator):
        """
        check if all elements in an iterator are identical
        """
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == rest for rest in iterator)
    
    def sites_nested_list_for_mm(groupd_genomic_keys_per_mm_dict,peptide):
        """
        create a dictionary mm:nested_list
        where each element in nested_list is a list of all genomic_keys of the specific mm_type that are represented by one instance of the peptide
        """

        sites_nested_list_for_mm ={}
        for mm in all_mm:
            groupd_genomic_keys_for_mm = groupd_genomic_keys_per_mm_dict[mm]
            sites_nested_list_per_pep_per_mm = [sorted(l) for l in groupd_genomic_keys_for_mm[peptide]]
            sites_nested_list_for_mm.update({mm:sites_nested_list_per_pep_per_mm})
        return sites_nested_list_for_mm
    
    #for edited peps, check carefully that this informative peptides are exactly from the same locus, and not different locus on the same gene that yield similar peptides.
    if row['editing_sites_inferred']>0:
        if all_instances_edited(all_peps[all_peps['peptide']==row['peptide']]):
            if row['informative_peptide']:
                return True #if peptide is informative in the regular sence (having one source) it will also bo informative in the weaker sence of origin from one exon
            else:    
                sites_nested_list_for_mm = sites_nested_list_for_mm(groupd_genomic_keys_per_mm_dict,row['peptide'])
                if all(checkEqual(sites_nested_list_for_mm[mm]) for mm in all_mm):
                    return True #if all sites represented by a peptide sequences are of the same genomic location then the peptide is "genomic informative"
                else:
                    return False
        else:
            return False #if one of the peptide instances are not from an edited sequence, surly the peptide is not informative for the editing sites from the edited sequence
    
    #for unedited peptides - a less stricked condition - peptides instances are from the same gene
    else:
        all_genes_for_pep = set(list(all_peps[all_peps['peptide'] == row['peptide']]['HGNC_protein_id']))
        if len(all_genes_for_pep) == 1:
            return True
        else:
            return False

def grouped_genomic_keys_for_peptide_instance(edited_peps):
    """
    this function creates a pandas series in which peptides are keys and values are nested lists:
    each element in the list is a list of genomic keys that are inferred by one instance of that peptied
    the function dose so for each mm type and returns a dictionary of mm:series
    
    the input is a dataframe in which all peptides are edited and conclusive with regard to specific editing(s) (editing_sites_inferred > 0)
    """
    groupd_genomic_keys_per_mm_dict = {}
    
    for i,mm in enumerate(all_mm):
        temp_df = edited_peps.copy()
        temp_df['inferred_genomic_sites'] = temp_df.apply(lambda row: row['genomic_keys'][i], axis = 1)
        temp_df = temp_df.groupby('peptide')['inferred_genomic_sites'].apply(list)
        groupd_genomic_keys_per_mm_dict.update({mm:temp_df})
        
    return groupd_genomic_keys_per_mm_dict
    

    
    
if __name__ == '__main__':
  
#    ref_Gene_file = 'E:/RNA_editing_Large_files/human_editom/hg38_refGene.txt'
#    peps_path = 'C:/Users/shosh/OneDrive/Desktop/test/results_from_sim_2mc_6minl_4600maxm_20maxes/'
#    peps_file = 'peps_from_sim_fa.pickle'
#    sites_file = 'E:/RNA_editing_Large_files/human_editom/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor.txt'
#    path_to_refGene_file = 'E:/RNA_editing_Large_files/human_editom/hg38_refGene.txt'
    
    peps_path = sys.argv[1]
    peps_file = sys.argv[2]
    sites_file = sys.argv[3]
    check_genomic_inf = eval(sys.argv[4])
#    
# =============================================================================
#     path_to_refGene_file = sys.argv[4]    
#     print('Reading refGene table')
#     refGene_df = read_refGene_df(path_to_refGene_file)
#     gene_to_variants = refGene_df.groupby('gene')['variant'].apply(list)
# =============================================================================    
  
    print('Reading sites table')
    sites_df = read_editing_sites_wrt_coding_seqs(sites_file)
    sites_df.set_index('coding_key_base1', inplace = True)
    
    for mm in all_mm:
        print('Total number of ' + mm + ' sites = ' + str(len(set(list(sites_df[sites_df['mm_type']==mm]['genomic_key_base1'])))) )
    
    
    print('Reading proteomics simulation results')
    peps_df = pd.read_pickle(peps_path + peps_file)
#    peps_df = peps_df.apply(lambda row: split_gene_and_prot_names(row, 'seq_id'), axis = 1)
    print('Getting genomic coordinates of sites for each edited peptide')

    edited_peptides_and_native_versions = peps_df[peps_df['sites_in_permutation_range']>0]
    edited_peptides_and_native_versions = edited_peptides_and_native_versions.apply(lambda row: split_gene_and_prot_names(row, 'seq_id'), axis = 1)
    
    edited_peptides_and_native_versions['genomic_keys'] = edited_peptides_and_native_versions.apply(lambda row: genomic_keys_for_sites(row, sites_df), axis = 1)
#    edited_peps = edited_peptides_and_native_versions[edited_peptides_and_native_versions[']]
    edited_peps = edited_peptides_and_native_versions[edited_peptides_and_native_versions['editing_sites_inferred']>0].copy()
    #edited peptides that also conclusively indicate that specific editing(s) occured (peptide coulde be generated by few editing)
#    edited_peps = peps_df[peps_df['editing_sites_inferred']>0].copy()
#    edited_peps = edited_peps.join(sites_df[['ucsc_id','HGNC_protein_id']].drop_duplicates('ucsc_id').set_index('ucsc_id'), how = 'inner', on = 'seq_id')
    
    edited_peptides_and_native_versions_list = list(set(list(edited_peptides_and_native_versions['peptide'])))
    
    #a slice of all the peptides that are or identical to some edited peptides (just for better preformance when calculating for edited sites if genomic informative or not)
    edited_peptides_and_native_versions_duplicates = peps_df[peps_df['peptide'].isin(edited_peptides_and_native_versions_list)]
    edited_peptides_and_native_versions_duplicates = edited_peptides_and_native_versions_duplicates.apply(lambda row: split_gene_and_prot_names(row, 'seq_id'), axis = 1)

    #each element in the dict is series of nested lists for each peptide - each elemnt is a list of all genomic keys (for specific mm) represented by one instance of the peptide. 
    groupd_genomic_keys_per_mm_dict = grouped_genomic_keys_for_peptide_instance(edited_peps)
 
#    edited_peps['genomic_informative'] = edited_peps.apply(lambda row: genomic_informative(row,groupd_genomic_keys_per_mm_dict,edited_peps_duplicates), axis = 1)
    if check_genomic_inf:
        print('Infer genomic-informative peptides')    
        edited_peptides_and_native_versions['genomic_informative'] = edited_peptides_and_native_versions.apply(lambda row: genomic_informative(row,groupd_genomic_keys_per_mm_dict,edited_peptides_and_native_versions_duplicates), axis = 1)
    else:
        print('Genomic informative peptides are infered based on proteom informative status')
        edited_peptides_and_native_versions['genomic_informative'] = edited_peptides_and_native_versions.apply(lambda row: row['informative_peptide'], axis=1)
        
    #creating edited_peps dataframe - now with genomic informative field
    edited_peps = edited_peptides_and_native_versions[edited_peptides_and_native_versions['editing_sites_inferred']>0].copy()
    
    print('Getting detectable sites')
    detectable_edited_peps = edited_peps[edited_peps['genomic_informative']]
    
    detectable_sites_per_mm = {}
    for i,mm in enumerate(all_mm):
        det_sites = list(set([site for group in list(detectable_edited_peps['genomic_keys']) for site in group[i]]))
        detectable_sites_per_mm.update({mm:det_sites})
        print(str(len(det_sites)) + ' ' + mm + ' sites are detectable in at least one version of a gene')
        
    print('Writing (genomic) informative edited peptides table')
    edited_peptides_and_native_versions.to_pickle(peps_path + 'edited_peptides_and_native_versions.pickle')
    edited_peptides_and_native_versions.to_csv(peps_path + 'edited_peptides_and_native_versions.txt', sep = '\t', index = False)
    edited_peptides_and_native_versions.to_excel(peps_path + 'edited_peptides_and_native_versions.xlsx')
    print('Post processing finished :)')
#    all_peps.to_pickle(peps_path + 'new_' + peps_file)
 
