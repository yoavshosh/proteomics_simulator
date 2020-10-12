import pandas as pd
import sys


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

    columns = ['gene_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand']
    data = []
    
    with open(anovar_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            fields = line.split("\t")
            for i in range(len(fields)):
                fields[i] = fields[i].replace('\n','').replace('\r','')
            data.append(fields)
            
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(lambda row: split_gene_and_prot_names(row, 'gene_name'), axis = 1)
    df['genomic_key_base1'] = df.apply(lambda row: row['HGNC_protein_id'] + ';' + row['chromosome']+'_'+row['strand']+'_'+row['genomic_position_base1'], axis = 1)
    df['coding_key_base1'] = df.apply(lambda row: row['ucsc_id']+'_'+row['position_base1'],axis = 1)
    return df


def add_data_to_sites_tble(row, full_data):

    #adding aa change
    
    variants_changes_str = full_data.loc[full_data['genomic_key_base1_no_gene_name'] == row['genomic_key_base1'].split(';')[-1]][7].values[0]
    variants_changes_list =  variants_changes_str.split(':')
    for i in range(len(variants_changes_list)):
        if variants_changes_list[i] == row['ucsc_id']:
            aa_change = ''.join([j for j in variants_changes_list[i+3].split('p.')[1].split(',')[0] if not j.isdigit()])
            break
    
    #adding average editing level
    editing_level_data_str = full_data.loc[full_data['genomic_key_base1_no_gene_name'] == row['genomic_key_base1'].split(';')[-1]][4].values[0]
    editing_level_data_list =  editing_level_data_str.split(';')
    editing_level = float(editing_level_data_list[3])
    
    row['aa_change'] = aa_change
    row['editing_level'] = editing_level
    
    return row
            



if __name__ == '__main__':
    
    sites_file = sys.argv[1]
    full_data_file = sys.argv[2]
#    
#    sites_file = 'E:/RNA_editing_Large_files/human_editom/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor.txt'
#    full_data_file = 'E:/RNA_editing_Large_files/human_editom/allTissues_filtered_blat5bp_clust_100_wTissueNumber_wSums_summary_editLvl0.01inOneTissue_allFields_NoWxsLvl0.001_tissues_maxEditLvl_Min100reads.txt'
    
    sites_df = read_editing_sites_wrt_coding_seqs(sites_file)
    
    full_data = pd.read_csv(full_data_file, sep = '\t', header = None)
    full_data = full_data[full_data[8] == 'nonsynonymous SNV']
    
    full_data['genomic_key_base1_no_gene_name'] = full_data.apply(lambda row: row[0]+'_'+row[5]+'_'+str(row[2]), axis = 1)
    
    sites_df = sites_df.apply(lambda row: add_data_to_sites_tble(row, full_data), axis = 1)
    
    columns_to_write = ['gene_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand', 'aa_change', 'editing_level']
    
    sites_df.to_csv(sites_file[0:-4]+'_with_additional_data.txt', sep = '\t', header  = False, index = False, columns = columns_to_write, line_terminator = '\n')
