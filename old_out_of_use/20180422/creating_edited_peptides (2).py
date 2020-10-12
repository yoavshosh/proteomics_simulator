import re
import itertools as it
import sys
import timeit
import pandas as pd
import numpy as np
from all_edited_peptides import create_edited_rna_peptides
from rna_peptide import rna_peptide
from collections import deque
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_rna

#sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/' +'test.txt', 'w')

a2g_header_regex = re.compile(r'(?<=a2g:\s).*?]')
c2t_header_regex = re.compile(r'(?<=c2t:\s).*?]')
digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
look_ahead_for_editnig = 3
look_behind_for_editing = 6
input_path = 'C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/test_files/'
input_fasta = 'in_frame_rna_from_orfs_squ_edited_plusminus_whole_comp_ReverseComplimentMinustrand'

edited_peptides_from_rna_seq(input_path,input_fasta,a2g_header_regex,c2t_header_regex)

def edited_peptides_from_rna_seq(input_path,input_fasta,a2g_header_regex,c2t_header_regex):
    
    for record in SeqIO.parse(input_path + input_fasta + '.fasta', "fasta"):
        a2g_sites = find_by_regex_in_header(record.description,a2g_header_regex)
        c2t_sites = find_by_regex_in_header(record.description,c2t_header_regex)
        print(record.id)
        print(a2g_sites)
        print(c2t_sites)
        
        
        

    
"""
find regex (regex - pre compiled regex) in header
"""
def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'
    