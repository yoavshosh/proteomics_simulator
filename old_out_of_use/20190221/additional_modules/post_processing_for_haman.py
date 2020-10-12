import sys
import numpy as np
import pandas as pd
from Bio import SeqIO


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


def genomic_keys_for_sites(row, sites_df, mm_type = 'AG'):
    
    """
    for each site represented by peptide, get the genomic key from sites_df
    creates a column with a list of genomic keys corresponding to the sites in the list of real_inferred_a2g_sites_coor_base1 in peps df
    genomic keys are: crhomosome_strand_site-position-wrt-to-chromosome
    """

    protein_id = row['seq_id']
    if mm_type == 'AG':
        protein_sites_keys = [protein_id + '_' + str(site) for site in row['real_inferred_a2g_sites_coor_base1']]
    elif mm_type == 'CT':
        protein_sites_keys = [protein_id + '_' + str(site) for site in row['real_inferred_c2t_sites_coor_base1']]
    try:
        genomic_sites_keys = [sites_df.loc[key, 'genomic_key_base1'] for key in protein_sites_keys]
        return genomic_sites_keys
    except KeyError:
        print(key + ' was not found on sites df. ' + protein_id + protein_sites_keys)
        

def genomic_informative(row, groupd_genomic_keys_by_pep_a2g, groupd_genomic_keys_by_pep_c2t, all_peps):
    
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
        check 
        """
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == rest for rest in iterator)
    
    
    if all_instances_edited(all_peps[all_peps['peptide']==row['peptide']]):
        if not row['informative_peptide']:
            all_peptides_sites_representations_a2g = [sorted(l) for l in groupd_genomic_keys_by_pep_a2g[row['peptide']]]  #sorting each list in a nested list - each elemnt is a list of all a2g genomic keys represented by some instance of the peptide
            all_peptides_sites_representations_c2t = [sorted(l) for l in groupd_genomic_keys_by_pep_c2t[row['peptide']]]  #same for c2t
            if checkEqual(all_peptides_sites_representations_a2g) and checkEqual(all_peptides_sites_representations_c2t):
                return True #if all sites represented by a peptide sequences are of the same genomic location then the peptide is "genomic informative"
            else:
                return False
        else:
            return True #if peptide is informative in the regular sence (having one source) it will also bo informative in the weaker sence of origin from one exon
    else:
        return False #this function is intended to operate on edited peptides. so if one of the peptide instances are not from an edited sequence, surly the peptide is not informative for the editing sites from the edited sequence
        


    
    
    
if __name__ == '__main__':
    
    peps_path = sys.argv[1]
    peps_file = sys.argv[2]
    sites_file = sys.argv[3]
#    path_to_refGene_file = sys.argv[4]
    
#    ref_Gene_file = 'E:/RNA_editing_Large_files/human_editom/hg38_refGene.txt'
#    peps_path = 'E:/RNA_editing_Large_files/human_editom/proteomics_simulation/results_from_all_variants_c2t_a2g_3mc_6minl_5500maxm_20maxes/'
#    peps_file = 'peps_from_all_variants_c2t_a2g_fasta.pickle'
#    sites_file = 'E:/RNA_editing_Large_files/human_editom/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor.txt'
#    path_to_refGene_file = 'E:/RNA_editing_Large_files/human_editom/hg38_refGene.txt'
    
# =============================================================================
#     print('Reading refGene table')
#     refGene_df = read_refGene_df(path_to_refGene_file)
#     gene_to_variants = refGene_df.groupby('gene')['variant'].apply(list)
# =============================================================================
    
    
    print('Reading sites table')
    sites_df = read_editing_sites_wrt_coding_seqs(sites_file)
    sites_df.set_index('coding_key_base1', inplace = True)
    all_A2G_sites_genomic_keys = set(list(sites_df[sites_df['mm_type']=='AG']['genomic_key_base1']))
    print('Total number of A2G recoding sites = ' + str(len(all_A2G_sites_genomic_keys)))
    all_C2T_sites_genomic_keys = set(list(sites_df[sites_df['mm_type']=='CT']['genomic_key_base1']))
    print('Total number of C2T recoding sites = ' + str(len(all_C2T_sites_genomic_keys)))

    
    print('Reading proteomics simulation results')
    peps_df = pd.read_pickle(peps_path + peps_file)
    print('Getting genomic coordinates of sites for each edited peptide')
    peps_df['genomic_a2g_sites_keys'] = peps_df.apply(lambda row: genomic_keys_for_sites(row, sites_df, mm_type='AG'), axis = 1)
    peps_df['genomic_c2t_sites_keys'] = peps_df.apply(lambda row: genomic_keys_for_sites(row, sites_df, mm_type='CT'), axis = 1)
    peps_df['number_of_a2g_sites'] = peps_df.apply(lambda row: len(row['real_inferred_a2g_sites_coor_base1']), axis = 1)
    peps_df['number_of_c2t_sites'] = peps_df.apply(lambda row: len(row['real_inferred_c2t_sites_coor_base1']), axis = 1)
    
    #lists of genomic keys of all sites (A2G and C2T)
    a2g_edited_peps_list = list(set(list(peps_df[peps_df['number_of_a2g_sites']>0].copy()['peptide'])))
    c2t_edited_peps_list = list(set(list(peps_df[peps_df['number_of_c2t_sites']>0].copy()['peptide'])))
    #and all mm types sites together
    edited_peps_list = list(set(a2g_edited_peps_list+c2t_edited_peps_list))
    
    #a slice of all the peptides dataframe containing all peptides that are or identical to edited peptides (just for better preformance when calculating for edited sites if genomic informative or not)
    edited_peps_duplicates = peps_df[peps_df['peptide'].isin(edited_peps_list)]

    #series of nested lists for each peptide - each elemnt is a list of all a2g/c2t genomic keys represented by some instance of the peptide
    groupd_genomic_keys_by_pep_a2g = peps_df[peps_df['editing_sites_inferred']>0].groupby('peptide')['genomic_a2g_sites_keys'].apply(list)
    groupd_genomic_keys_by_pep_c2t = peps_df[peps_df['editing_sites_inferred']>0].groupby('peptide')['genomic_c2t_sites_keys'].apply(list)
 
    print('Infer genomic-informative peptides')
    edited_peps = peps_df[peps_df['editing_sites_inferred']>0].copy()
    edited_peps['genomic_informative'] = edited_peps.apply(lambda row: genomic_informative(row,groupd_genomic_keys_by_pep_a2g,groupd_genomic_keys_by_pep_c2t,edited_peps_duplicates), axis = 1)
    
    print('Getting detectable sites')
    detectable_edited_peps = edited_peps[edited_peps['genomic_informative']]
    detectable_a2g_editing_sites = list(set([site for group in list(detectable_edited_peps['genomic_a2g_sites_keys']) for site in group]))
    detectable_c2t_editing_sites = list(set([site for group in list(detectable_edited_peps['genomic_c2t_sites_keys']) for site in group]))
    print(str(len(detectable_a2g_editing_sites)) + ' A2G sites are detectable in at least one version of a gene')
    print(str(len(detectable_c2t_editing_sites)) + ' C2T sites are detectable in at least one version of a gene')

    print('Writing (genomic) informative edited peptides table')
    edited_peps.to_pickle(peps_path + 'edited_peptides.pickle')
    edited_peps.to_excel(peps_path + 'edited_peptides.xlsx')
    
