import re
import os
import sys
import pandas as pd
import numbers
import operator
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis

all_mm = ['AG','AC','AT','CA','CG','CT','GA','GC','GT','TA','TG','TC']

def slice_anovar_site_per_mm(anovar_df, mm_list = ['AG','CT']):
    
    slices_dict = {}
    
    for mm in all_mm:
        if mm in mm_list:
            anovar_slice = anovar_df[anovar_df['mm_type']==mm]
            anovar_slice.set_index('ucsc_id', inplace = True)
            slices_dict.update({mm:anovar_slice})
        else:
            slices_dict.update({mm:None})
    
    return slices_dict
            

def get_genes_variants_length(refGene_file):
    
    """
    input - a refGene file containing data regarding variants of genes
    """
    def variant_length(row):
        starts = eval('['+row['exonStarts']+']')
        ends = eval('['+row['exonEnds']+']')
        l=0
        for st, en in zip(starts,ends):
            l+=en-st
        return l
        
    #see columns description in: http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=dm3&hgta_group=genes&hgta_track=refGene&hgta_table=refGene&hgta_doSchema=describe+table+schema
    columns = ['bin','variant_name', 'chromosome', 'strand',
               'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
               'exonCount',  'exonStarts', 'exonEnds',
               'score','gene_name','cdsStartStat','cdsEndStat','exonFrames']
    data = []
    
    with open(refGene_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            fields = line.split("\t")
            data.append(fields)
            
    df = pd.DataFrame(data = data, columns = columns)
    df['variant_length'] = df.apply(lambda row: variant_length(row), axis = 1)
    df['variant_length_tuple'] = df.apply(lambda row: (row['variant_name'],row['variant_length']),axis =1)
    
    #create a list of longest variants gene:[(variant1,length1),.....]
    variant_length_for_gene = df.groupby('gene_name')['variant_length_tuple'].apply(list) #Series index=gene, value=[(variant1,length1),.....]
    longest_variants = []
    for variants in variant_length_for_gene:
        variants.sort(key=operator.itemgetter(1))
        longest_variants.append(variants[-1][0]) #sort variants list for gene by ascending order and append last variant name
        
    return df,longest_variants


"""
find regex (regex - pre compiled regex) in string
"""
def find_by_regex_in_str(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'



def split_gene_and_prot_names(row):
    names = row['variant_name'].split(';')
    row['ucsc_id'] = names[0]
    row['HGNC_protein_id'] = names[1]
    return row



def read_editing_sites_wrt_coding_seqs(anovar_file):
    """
    read dataframe from anovar output
    read file with 4 columns as in columns and translate return dataframe
    """
         
    columns = ['variant_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand', 'genomic_key_base1','coding_key_base1']
    data = []
    
    with open(anovar_file, "r") as f:
        content = f.readlines()
        for i, line in enumerate(content):
            fields = line.split("\t")
            fields[3] = fields[3].replace('\n','')
            fields.append(str(fields[4]+'_'+fields[7]+'_'+fields[6]))
            fields.append(str(fields[0])+'_'+str(fields[2]))
            data.append(fields)
            
    df = pd.DataFrame(data = data, columns = columns)
    df = df.apply(split_gene_and_prot_names, axis = 1)
    return df

def create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs_final_sites_dfs(fasta_file,output_name,out_path,mm_df_dict,variants_to_use = []):
    """
    input - coding sequences as fasta file
            
            sites (wrt to coding sequence) dataframe - result of read_editing_sites_wrt_coding_seqs
            after ucsc_id column is set to index 
            different dataframes for different mm types
    
    output - fasta file in the format of proteomics simulator 
             some of the values in the header will be useless because the input includes that coding sequences
             so this function does not trim the sequences.
    """
    
    n_bad = 0
    n_good = 0
    sites_good = 0
    sites_bad = 0
    
    writer =  FastaWriter(open(out_path + output_name + '.fasta' , 'w'), wrap=None)
    writer_bad = FastaWriter(open(out_path + 'bad_seqs_' + output_name + '.fasta' , 'w'), wrap=None)
    writer.write_header()
    writer_bad.write_header()
    
    for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
        
        mm_loc_dict = {}
        
        rec_id = record.id.split(';')[0]
        use_variant = True
        
        if len(variants_to_use): #if a not-empty list is passed for variants_to_use, flag variants that are not in list so they will not be included in uotput
            if rec_id not in variants_to_use:
                use_variant = False
        
        if use_variant:
            
            final_sequence = str(record.seq[0:-3]).replace('a','A').replace('g','G').replace('t','T').replace('c','C')
            
            for mm in all_mm:
                if mm_df_dict[mm] is None:
                    mm_list = []
                else:
                    sites = mm_df_dict[mm]
                    try:
                        mm_list = [int(k)-1 for k in sites.loc[[rec_id]]['position_base1']]
                    except KeyError:
                        mm_list = []
                mm_loc_dict.update({mm:mm_list})

            
            prot_start_nuc = 1
            prot_end_nuc = len(final_sequence) - 3
            prot_start = 'first_met_in_original_orf'
            prot_end = 'original_sense_strand_orf_end'
            strand = '+'
            orf_start = 1
            orf_end = len(record.seq) - 3
             
            mm_str = ''
            for mm in mm_loc_dict:
                mm_str+= '| '+mm+'_base0: '+ str(mm_loc_dict[mm])
                    
            description_str = mm_str + ' | prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + ' | strand: ' + strand + ' | prot_start_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc) + ' | original_orf_start: ' + str(orf_start) + ' | original_orf_end: ' + str(orf_end)
        
            if (Seq(str(record.seq[0:3]), generic_dna).translate() != 'M' or Seq(str(record.seq[-3:len(record.seq)]), generic_dna).translate() != '*') or '*' in Seq(str(final_sequence), generic_dna).translate():
                writer_bad.write_record(record)
                n_bad+=1
                sites_bad+=sum([len(mm_loc_dict[mm]) for mm in all_mm])
            else:
                if len(final_sequence)%3:
                    final_sequence=final_sequence[0:-len(final_sequence)%3]
                current_record = SeqRecord(Seq(final_sequence,generic_dna), id = rec_id, description = description_str)
                writer.write_record(current_record)
                n_good+=1
                sites_good+=sum([len(mm_loc_dict[mm]) for mm in all_mm])
            
    writer_bad.write_footer()
    writer.write_footer()
    
    print(str(n_good) + ' good sequence with ' + str(sites_good) + 'sites')
    print(str(n_bad) + ' bad sequence with ' + str(sites_bad) + 'sites')
    


def read_fasta_dict(fasta_file):
    """
    create a dictionary of transcript_name:sequence
    from human fasta file 
    """    
    dict = {}
    for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
        dict.update({record.id.split(';')[0]:str(record.seq)})
    return dict
        

def check_frame(fasta_file):
    """
    print sequences id's for sequences that are not multiplet of 3
    """
    n=0
    for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
        if len(str(record.seq))%3:
            print(record.id)
        else:
            n+=1
    print(str(n))
            
    
def check_len(fasta_file,length):
    """
    print seuqneces id's for sequences that are shorter than specified length
    """
    for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
        if len(str(record.seq))<length:
            print(record.id)

    
if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Preparation of proteomics simulator input')
    run_parser = parser.add_argument_group('Run mqpar prep')
    run_parser.add_argument('-i_fasta', dest='fasta_file', action='store', required = True, help='Path to fasta file')
    run_parser.add_argument('-i_anovar', dest='anovar_file', action='store', required = True, help='Path to fasta file')
    run_parser.add_argument('-use_longest', dest='use_longest', action='store', default = 'False', help='True if output sould include only longest variant per each gene')
    run_parser.add_argument('-refGene_file', dest='refGene_file', action='store', default = '', help='refGene_file with data regarding all genes and thier variants')
    run_parser.add_argument('-o', dest='output_name', action='store', default = 'output.fasta', help='fasta output name')
    run_parser.add_argument('-mm_types', dest='mm_types', action='store', nargs = '+' ,default = 'AG CT', help='mm types from anovar file to include in fasta file')
    arguments = parser.parse_args()
    
    fasta_file = arguments.fasta_file
    anovar_file = arguments.anovar_file
    use_longest = eval(arguments.use_longest)
    refGene_file = arguments.refGene_file
    output_name = arguments.output_name
    mm_types = arguments.mm_types
    
    if use_longest:
        assert refGene_file != '', 'use_longest is the true but no refGene_file as input to determine longest variants per gene'
    
    
#    fasta_file = 'E:/RNA_editing_Large_files/human_editom/human_coding_sequences.fasta'
#    anovar_file = 'E:/RNA_editing_Large_files/human_editom/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor.txt'
#    output_name = 'test.fasta'
#    only_recoding = True

    input_path = '/'.join(fasta_file.split('/')[:-1]) + '/'

    if not os.path.exists(input_path + 'proteomics_simulation/'):
        os.makedirs(input_path + 'proteomics_simulation/')
    out_path = input_path + 'proteomics_simulation/'     

    print('Reading and processing Anovar tabel')
    sites = read_editing_sites_wrt_coding_seqs(anovar_file)
    all_mm_df_dict = slice_anovar_site_per_mm(sites, mm_list = mm_types)
    
    
    if use_longest:
        print('Preparing list of longest variants')
        variants_df, longest_variants = get_genes_variants_length(refGene_file)
        print('Writing fasta file in the proteomics simulator format')
        create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs_final_sites_dfs(fasta_file,output_name,out_path,all_mm_df_dict,variants_to_use=longest_variants)            
    else:
        print('Writing fasta file in the proteomics simulator format')
        create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs_final_sites_dfs(fasta_file,output_name,out_path,all_mm_df_dict)
        
