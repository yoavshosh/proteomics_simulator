import pandas as pd
import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
        

def split_gene_and_prot_names(row):
    names = row['variant_name'].split(';')
    row['ucsc_id'] = names[0]
    row['HGNC_protein_id'] = names[1]
    return row



def read_non_editing_sites_wrt_coding_seqs(anovar_file):
    """
    read dataframe from anovar output
    read file with 4 columns as in columns and translate return dataframe
    """
         
    columns = ['variant_name', 'position_base0', 'position_base1', 'exon_number', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand', 'genomic_key_base1','coding_key_base1']
    data = []
    
    with open(anovar_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            line = line.replace('\n','')
            fields = line.split("\t")
            fields.append(str(fields[4]+'_'+fields[7]+'_'+fields[6]))
            fields.append(str(fields[0])+'_'+str(fields[2]))
            data.append(fields)
            
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(split_gene_and_prot_names, axis = 1)
    return df


def sequences_dictionary(fasta_file):
    seqs_dict = {}
    for rec in SeqIO.parse(open(fasta_file, "r"), 'fasta'):
        seqs_dict.update({rec.id:str(rec.seq)})
    return seqs_dict
  
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Preparation of proteomics simulator input')
    run_parser = parser.add_argument_group('run random AG table')
    run_parser.add_argument('-fasta', dest='fasta_file', action='store', required = True, help='Path to fasta file')
    run_parser.add_argument('-sites', dest='anovar_file', action='store', required = True, help='Path to fasta file')
    run_parser.add_argument('-o', dest='output_name', action='store', default = 'random_ag_editings', help='output table name')
    run_parser.add_argument('-n', dest='number_of_sites', action='store',default = '974', help='number of AG mm to draw')
    arguments = parser.parse_args()
        
    fasta_file = arguments.fasta_file
    sites_file = arguments.anovar_file
    output_name = arguments.output_name
    number_of_sites = eval(arguments.number_of_sites)
    
#    sites_file = 'E:/RNA_editing_Large_files/human_editom/non_editing_a_sites_incomplete.txt'
#    fasta_file = 'E:/RNA_editing_Large_files/human_editom/ag_all_variants.fasta'
#    output_name = 'random_ag_editings'
#    number_of_sites = 974
    
    seqs_dict = sequences_dictionary(fasta_file)
    real_sites_df = read_non_editing_sites_wrt_coding_seqs(sites_file)
    genomic_sites_set = list(set(list(real_sites_df['genomic_key_base1'])))
    
    random_genomic_sites = []
    while len(random_genomic_sites) < number_of_sites:
        
        genomic_site = random.choice(genomic_sites_set)
        if genomic_site not in random_genomic_sites:
            
            recoding = False
            variants_for_gene =  real_sites_df[real_sites_df['genomic_key_base1'] == genomic_site].copy()
           
            original_codons_for_variants = []
            edited_codons_for_variants = []
            
            #running over all sequences containing chosen site to determine if recoding or not
            for index, row in variants_for_gene.iterrows():
                try:
                    seq = seqs_dict[row['ucsc_id']]
                except KeyError:
                    break
                
                base0_coding_position = int(row['position_base0'])
                
                #QA - break script if one of the sites is not an A (assumption is that the list contain only A's)
                if base0_coding_position>=len(seq):
                    print('site ' + row['genomic_key_base1'] + ' is out of recoding range for variant ' + row['ucsc_id'] + ' by ' +str(base0_coding_position-len(seq)))
                    break
                if seq[base0_coding_position] != 'A':
                    print('site ' + row['genomic_key_base1'] + ' is not an A site')
                    break

                #determine original and edited codons
                if base0_coding_position%3 == 0:
                    original_codon = seq[base0_coding_position:base0_coding_position+3]
                    edited_codon = 'G'+seq[base0_coding_position+1:base0_coding_position+3]
                elif base0_coding_position%3 == 1:
                    original_codon = seq[base0_coding_position-1:base0_coding_position+2]
                    edited_codon = seq[base0_coding_position-1] + 'G' + seq[base0_coding_position+1]
                elif base0_coding_position%3 == 2:
                    original_codon = seq[base0_coding_position-2:base0_coding_position+1]
                    edited_codon = seq[base0_coding_position-2:base0_coding_position] + 'G'
                    
                #orinial and edited aa based on codons and determine if chosen site is recoding in at list one variant
                original_aa = Seq(original_codon, generic_dna).translate()
                edited_aa = Seq(edited_codon, generic_dna).translate()
                
                if original_aa != edited_aa:
                    recoding = True
                    original_codons_for_variants.append(original_codon)
                    edited_codons_for_variants.append(edited_codon)
                else:
                    original_codons_for_variants.append(original_codon)
                    edited_codons_for_variants.append(edited_codon)
            
            if recoding:
                random_genomic_sites.append(genomic_site)
                if len(set(original_codons_for_variants))>1:
                    print(row['HGNC_protein_id'] + ' site:' + genomic_site + ' is in multiple original codons')
                if len(set(edited_codons_for_variants))>1:
                    print(row['HGNC_protein_id'] + ' site:' + genomic_site + ' is in multiple edited codons')
    
    final_random_ag_sites_df = real_sites_df[real_sites_df['genomic_key_base1'].isin(random_genomic_sites)].copy()
    final_random_ag_sites_df.drop('exon_number', axis = 1, inplace = True) 
    final_random_ag_sites_df.loc[:,'mm_type'] = 'AG'
    columns = ['variant_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand']
    final_random_ag_sites_df.to_csv('/'.join(sites_file.split('/')[0:-1]) + output_name + 'txt', columns=columns)

