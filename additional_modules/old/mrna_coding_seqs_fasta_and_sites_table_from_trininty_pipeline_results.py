# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 15:11:08 2019

@author: shosh

From Editing sites table craeted by transcriptome, and fasta transcriptome conaining orfs,
create a table with relevant fields for proteomics simulation input fasta preperation
"""

import re
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, IUPAC 


def find_by_regex_in_line(line, regex):
    try:
        return regex.findall(line)[0]
    except:
        return []


def sites_coding_location(row, orfs_dictionaty):
    
    start = orfs_dictionaty[row['trinity']][0]
    end = orfs_dictionaty[row['trinity']][1]
    
    if row['strand'] == '+':
        editing_sites_in_coding_sequence = row['position_base1'] - start
    elif row['strand']== '-':
        editing_sites_in_coding_sequence = end - row['position_base1']
    
    return editing_sites_in_coding_sequence


def retrive_native_rna(record, sites_in_sequence):
    """
    this function create the native sequence
    given a sequence that is edited (in all strong sites, but anyway the function iterate over all sites) 
    and a dataframe of edited sites for the particular sequence (with mm_type column)
    """
    complementaries = {'A':'T','T':'A','G':'C','C':'G'}
    sequence = str(record.seq)
    
    for i, row in sites_in_sequence.iterrows():
        
        #determine native and edited nucs in input seq
        native_coding_nuc = row['mm'][0]
        if row['strand'] == '-':
            native_nuc_in_seq = complementaries[native_coding_nuc]
        elif row['strand'] == '+':
            native_nuc_in_seq = native_coding_nuc
    
        sequence = sequence[:row['position_base1']-1] + native_nuc_in_seq + sequence[row['position_base1']:]

    return sequence


def create_coding_mrna_fasta(input_path, input_fasta, sites_df):
    
    """
    create coding mrna from trinity fasta components file
    """    

    orfs_dictionaty = {}
    writer_transcripts = FastaWriter(open(input_path + 'mrna_coding_sequences_from_' + input_fasta , 'w'), wrap=None)
    writer_transcripts.write_header()

    for record in SeqIO.parse(open(input_path + input_fasta, "r"), "fasta"):
        
        rec_data = record.description.split('\t')
        strand = rec_data[6]
        start = int(rec_data[2])
        end = int(rec_data[4])
        protein = rec_data[-1].split('|')[1]
        sequence = str(record.seq)
        sites_for_sequence = sites_df[sites_df['trinity'] == record.id].copy()
        
        if len(sites_for_sequence):
            sequence = retrive_native_rna(record, sites_for_sequence)

        if strand == '+':
            coding_sequence = Seq(sequence[start-1:end], generic_dna)
        elif strand == '-':
            coding_sequence = Seq(sequence[start-1:end], generic_dna).reverse_complement()
        else:
            raise ValueError('No valid strand for record ' + str(record.id))
            break
        
        rec_id = record.id+';'+protein
        orfs_dictionaty.update({record.id:(start,end)})    
        writer_transcripts.write_record(SeqRecord(coding_sequence, id = rec_id+'\t', description=''))
    
    writer_transcripts.write_footer()
    
    return orfs_dictionaty
    

if __name__ == '__main__':

    input_path = sys.argv[1]
    input_fasta = sys.argv[2]
    sites_path = sys.argv[3]

    print('Reading sites file')
    editing_sites_columns = ['trinity','gene_name','position_base1','mm','dna_A_reads','dna_T_reads','dna_G_reads','dna_C_reads','rna_A_reads','rna_T_reads','rna_G_reads','rna_C_reads','nuc_in_transcriptome_file','total_rna_reads','total_dna_reads','p_val','aa_before','aa_after','analysis_type','orfs_length','editing_level','strand']
    sites_df = pd.read_csv(sites_path, names = editing_sites_columns, sep = '\t', index_col=False)
    
    print('creating coding mrna fasta file')
    orfs_dictionaty = create_coding_mrna_fasta(input_path, input_fasta, sites_df)
    
    print('creating sites file in new format')
    sites_df['coding_loc_base0'] = sites_df.apply(lambda row: sites_coding_location(row, orfs_dictionaty), axis = 1)
    sites_df['coding_loc_base1'] = sites_df.apply(lambda row: row['coding_loc_base0']+1, axis = 1)
    sites_df['position_base0'] = sites_df.apply(lambda row: row['position_base1']-1, axis = 1)
    sites_df['full_name'] = sites_df.apply(lambda row: row['trinity']+';'+row['gene_name'], axis = 1)
    sites_df['chromosome'] = sites_df.apply(lambda row: 'trinity', axis = 1)
    sites_df['aa_change'] = sites_df.apply(lambda row: row['aa_before']+row['aa_after'], axis = 1)
    
    new_sites_file_path = '/'.join(sites_path.split('.')[0].split('/')[:-1])+'/'
    new_sites_file_name = sites_path.split('/')[-1].split('.')[0] + '_new_format.txt'
    columns_to_wtire = ['full_name','coding_loc_base0','coding_loc_base1','mm','chromosome','position_base0','position_base1','strand','aa_change','editing_level']
    sites_df.to_csv(new_sites_file_path+new_sites_file_name,sep='\t',columns=columns_to_wtire,header=False,index=False)
