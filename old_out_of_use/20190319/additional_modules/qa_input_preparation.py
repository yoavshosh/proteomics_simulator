import re
import sys
import pandas as pd
import numbers
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#i_path = r'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
#o_path = i_path
#fasta_input = 'orfs_squ.fasta'
#xls_path = r'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
#xls_file = 'editing_sites.xlsx'
#xls_columns_list = [0,3,4]
#xls_sheet = 'Editing_events_in_Squid'

#strand_regex = re.compile('(?<=Strand\s)[^\s]+')
#start_nuc_regex = re.compile('(?<=OrfStart\s)[^\s]+')
#end_nuc_regex = re.compile('(?<=OrfEnd\s)[^\s]+')

path1 = 'C:/Users/user/Google_Drive/RNA_Editing/files/old/'
path2 = 'C:/Users/user/Google_Drive/RNA_Editing/files/'
fasta1 = 'in_frame_rna_rec_only_from_orfs_squ.fasta'
fasta2 = 'in_frame_rna_rec_only_from_orfs_squ.fasta'
#a,b,c, = compare_two_fastas(path1,fasta1,path2,fasta2)
def compare_two_fastas(path1,fasta1,path2,fasta2):
    
    dict1 = {}
    dict2 = {}
    changed_recs = {}
    recs_not_in_dict2 = []
    recs_not_in_dict1 = []
    
    for record in SeqIO.parse(open(path1 + fasta1, "r"), "fasta"):
        dict1.update({record.id:str(record.seq)})

    for record in SeqIO.parse(open(path2 + fasta2, "r"), "fasta"):
        dict2.update({record.id:str(record.seq)})

    for rec in dict1:
        try:
            if dict1[rec] == dict2[rec]:
                pass
            elif dict1[rec] != dict2[rec]:
                changed_recs.update({rec:abs(len(dict2[rec])-len(dict1[rec]))})
        except(KeyError):
            recs_not_in_dict2.append(rec)
            
    for rec in dict2:
        try:
            if dict1[rec] == dict2[rec]:
                pass
            elif dict1[rec] != dict2[rec] and rec not in changed_recs:
                changed_recs.update({rec:abs(len(dict2[rec])-len(dict1[rec]))})
        except(KeyError):
            recs_not_in_dict1.append(rec)  
            
    return changed_recs, recs_not_in_dict1, recs_not_in_dict2
        
                
                

"""
import editng sites from xls to a pd framework
"""
#df = import_es_from_xls(xls_path,xls_file,xls_sheet = xls_sheet, xls_columns_list = xls_columns_list, only_recoding = True)
def import_es_from_xls(xls_path,xls_file,xls_sheet = 0, xls_columns_list = [0,3,4], only_recoding = True):
    df = pd.read_excel(xls_path + xls_file, index_col = 0, sheet_name = xls_sheet, usecols = xls_columns_list)
    if only_recoding:
        df = df.loc[df[df.columns[1]] != 'syn']
    df = df.groupby(df.index).agg({list(df.columns)[0]:lambda x: list(x)})
    return df



"""
find regex (regex - pre compiled regex) in header
"""
def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'
  
 
    
"""
find nucleotide location of first occurance of methionine within a given sequence
"""
def find_first_met(sequence):
    first_met = 'no_met'
    for i in range(0,len(sequence),3):
        if Seq(str(sequence[i:i+3]), generic_rna).translate() == 'M':
            first_met = i
            break     
    return first_met



"""
find first stop codon after first determined start codon (forst methionine or if not found - orf_start from original fasta)
"""
def find_first_stop(sequence):
    first_stop = 'no_stop'
    for i in range(0,len(sequence),3):
        if Seq(str(sequence[i:i+3]), generic_rna).translate() == '*':
            first_stop = i
            break     
    return first_stop


"""
check whether the codon after and before reading frame is stop
"""
def stop_beyond_edges(sequence, strand, start, end):
    
    if strand == '+':
        if len(sequence) > end+3:
            if Seq(str(sequence[end:end+3]), generic_rna).translate() == '*':
                seq_end = True
            else:
                seq_end = False
        else:
            seq_end = False
        
        if start > 3:
            if Seq(str(sequence[start-4:start]), generic_rna).translate() == '*':
                seq_beginning = True
            else:
                seq_beginning = False
        else:
            seq_beginning = False
                
            
    elif strand == '-':    
        if start > 3:
            if Seq(str(sequence[start-4:start-1]), generic_rna).reverse_complement().translate() == '*':
                seq_end = True
            else:
                seq_end = False
        else:
            seq_end = False
        
        if len(sequence) > end+3:
            if Seq(str(sequence[end:end+3]), generic_rna).reverse_complement().translate() == '*':
                seq_beginning = True
            else:
                seq_beginning = False
        else:
            seq_beginning = False
    
    return seq_beginning, seq_end
                

"""
determine real orf considering if codons beyond orf are stop
and considering stop codons and methionine within orf
let s be orf_start and e orf_end, M - methionin, * - stop_codon.  then:

sense strand    
...s...M(or several M)...e....                 -> no change in orf.                                                  prot_start=unknown     prot_end=unknown
...s.....................e....                 -> no change in orf.                                                  prot_start=unknown     prot_end=unknown
...s...*(or several *)...e....                 -> orf_end is first stop codon in sequence.                           prot_start=unknown     prot_end=first_stop_in_original_orf_that_is_not_befor_first_met
...s...M...*(or several * and M)....e...       -> orf_end is first stop codon in sequence after first M              prot_start=unknown     prot_end=first_stop_in_original_orf_that_is_not_befor_first_met
..*s...M(or several M)...e...                  -> orf_start is first M in sequence                                   prot_start=within_original_orf     prot_end=unknwon
..*s...M...*(or several * and M)....e...       -> orf_start is first M in sequence orf_end is first stop codon in sequence after first M, prot_start=within_original_orf, prot_end=first_stop_in_original_orf_that_is_not_befor_first_met
..*s...*(or several *)...e...                  -> orf_end is first stop codon in sequence. prot is flaged as bad     prot_start=outside_original_orf     prot_end=first_stop_in_original_orf_that_is_not_befor_first_met.         
...s.....................e*....                -> no change in orf.  prot is flaged as bad                           prot_start=outside_original_orf     prot_end=original_orf_end  
..*s.....................e*....                -> no change in orf.  prot is flaged as bad                           prot_start=outside_original_orf     prot_end=original_orf_end
..*s.....................e.....                -> no change in orf.  prot is flaged as bad                           prot_start=outside_original_orf     prot_end=unknown
all cases with stop codons within original orf, the present of a stop codon right after the original orf is not relevant

changes in orf are done within create_in_frame_rna_file function and consider the strand orientation.


"""
def rethink_reading_frame(orf_sequence, stop_at_seq_beginning, stop_at_seq_end):
    
    flag_bad = False
    subseq = orf_sequence
    first_met = find_first_met(orf_sequence)
    met_number = Seq(str(orf_sequence), generic_rna).translate().count('M')

    if stop_at_seq_beginning:
        if met_number == 1:
            prot_start = 'first_met_in_original_orf'
            subseq = orf_sequence[first_met:]
        elif met_number > 1:
            prot_start = 'within_original_orf'
            subseq = orf_sequence[first_met:]
        else:
            prot_start = 'outside_original_orf'
            flag_bad = True
    else:
        if first_met != 'no_met':
            prot_start = 'unknown'
        else:
            prot_start = 'outside_original_orf'
    
    first_stop = find_first_stop(subseq)
    
    if first_stop == 'no_stop':
        if stop_at_seq_end:
            prot_end = 'original_orf_end'
        else:
            prot_end = 'unknown'
    else:
        prot_end = 'first_stop_in_original_orf_that_is_not_befor_first_met'
        subseq = subseq[:first_stop]
        
    return subseq, first_met, first_stop, prot_start, prot_end, flag_bad
    

def reading_frame_rna(sequence,strand,start_nuc,end_nuc):
    
    in_farme_seq = Seq(str(sequence[start_nuc-1:end_nuc]),generic_dna)

    if strand == '+':
        final_rna = in_farme_seq
    elif strand == '-':
        final_rna = in_farme_seq.reverse_complement()
    else:
        final_rna = None
        
    return final_rna


def import_es_from_xls(xls_path,xls_file,xls_sheet = 0, xls_columns_list = [0,3,4], only_recoding = True):
    df = pd.read_excel(xls_path + xls_file, index_col = 0, sheet_name = xls_sheet, usecols = xls_columns_list)
    if only_recoding:
        df = df.loc[df[df.columns[1]] != 'syn']
    df = df.groupby(df.index).agg({list(df.columns)[0]:lambda x: list(x)})
    return df


def check_original_nuc(input_path, input_fasta, a2g_editing_sites_df):
    
    sites_not_a = {}
    
    strand_regex = re.compile('(?<=Strand\t)[^\t]+')
    
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        nucs_list = []
        strand = find_by_regex_in_header(record.description, strand_regex)
        try:
            for site in a2g_editing_sites_df.loc[record.id][0]:
                nucs_list.append(str(record.seq)[site-1])
            if strand == '-' and any(n in nucs_list for n in ['A','G']):
                sites_not_a.update({record.id:strand})
            if strand == '+' and any(n in nucs_list for n in ['T','C']):
                sites_not_a.update({record.id:strand})
        except KeyError:
            pass
    
    return sites_not_a
    
# =============================================================================
# strand_regex = re.compile('(?<=Strand\t)[^\t]+')
# start_nuc_regex = re.compile('(?<=OrfStart\t)[^\t]+')
# end_nuc_regex = re.compile('(?<=OrfEnd\t)[^\t]+')
# input_path = 'C:/Users/user/Google_Drive/RNA_Editing/files/'
# out_path = input_path
# only_recoding = True
# filter_sequences_with_stop_codons = True
# fasta_input = 'shahar_squ_orfs.fasta'
# 
# n_good = 0
# n_bad = 0
# 
# if only_recoding:
#     out_str = 'rec_only'
# else:
#     out_str = 'all_sites'
#  
# out_of_range_sites_file = open(input_path + 'out_of_frame_sites_' + out_str + '_from_' + fasta_input + '.txt', "w")
# writer =  FastaWriter(open(out_path + 'in_frame_rna_' + out_str + '_from_' + fasta_input , 'w'), wrap=None)
# writer_bad =  FastaWriter(open(out_path + 'bad_sequencese_' + out_str + '_from_' + fasta_input , 'w'), wrap=None)
# writer.write_header()
# writer_bad.write_header()
# 
# for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
#         
#     bad_seq = False
#     a2g_list = []
#     c2t_list = []
#     out_of_frame_a2g = []
#     out_of_frame_c2t = []
#     
#     strand = find_by_regex_in_header(record.description, strand_regex) #'+' if seq is read as is, '-' if reverse compliment is read
#     
#     #start_nuc and end_nuc are nucleotides locations within original sequence (regardless of strand) in input fasta stating the orfs 
#     orf_start = int(find_by_regex_in_header(record.description, start_nuc_regex))
#     orf_end = int(find_by_regex_in_header(record.description, end_nuc_regex))
#     
#     stop_at_seq_beginning, stop_at_seq_end = stop_beyond_edges(record.seq, strand, orf_start, orf_end)
# 
#     #get initial reading frame by the orfs stated in input fasta (input fasta orfs are blastx orfs prolonged until stop codons or until edges of aligned query sequences)
#     orf_sequence = reading_frame_rna(record.seq,strand,orf_start, orf_end)
#     if filter_sequences_with_stop_codons:
#         full_protein = Seq(str(orf_sequence), generic_rna).translate()
#         if str(full_protein).count('*'):
#             bad_seq = True
#     
#     #if stop cidins at seq edges, shorten reading frame and return new seq, first methionin location, and markers for wether or not seq edges are real or not
# #        final_sequence, first_met, prot_start, prot_end = rethink_reading_frame(orf_sequence)
#     final_sequence, first_met, first_stop, prot_start, prot_end, bad_orf = rethink_reading_frame(orf_sequence, stop_at_seq_beginning, stop_at_seq_end)
# 
#     prot_start_nuc = orf_start
#     prot_end_nuc = orf_end
#     #stop codon at beginning of sequence: sequence -> sequence from first methionin and updating in sense strand start_nuc
#     if stop_at_seq_beginning and not(prot_start == 'outside_original_orf'):
#         if strand == '+':
#             prot_start_nuc = prot_start_nuc+first_met
#         elif strand == '-':
#             prot_end_nuc = prot_end_nuc-first_met
#             
#         
#     #stop codon at end of sequence: sequence -> sequence untill stop codond and updating sense strand end_nuc
#     if first_stop != 'no_stop':
#         if strand == '+':
#             prot_end_nuc = prot_start_nuc+first_stop-1
#         elif strand == '-':
#             first_met_for_calc = first_met if stop_at_seq_beginning and not(prot_start == 'outside_original_orf') else 0
#             first_stop_for_clac = first_stop
#             prot_start_nuc = prot_end_nuc-first_met_for_calc-first_stop_for_clac+1
# 
#     #a2g editing site for sequence
#     if strand == '+':
#         try: #taking all in frame editing sites from sites dataframe
#             for k in sorted(a2g_editing_sites_df.at[record.id,a2g_editing_sites_df.columns[0]]):
#                 if prot_end_nuc>=k>=prot_start_nuc:
#                     a2g_list.append(k-prot_start_nuc)
#                 else: 
#                     out_of_frame_a2g.append(k)
#         except:
#             pass
#         
#         #c2t editing site for sequence
#         try:
#             for k in sorted(c2t_editing_sites_df.at[record.id,c2t_editing_sites_df.columns[0]]):
#                 if prot_end_nuc>=k>=prot_start_nuc:
#                     c2t_list.append(k-prot_start_nuc)
#                 else:
#                     out_of_frame_c2t.append(k)
#         except:
#             pass
#             
#         
#     elif strand == '-':
#         try:
#             for k in sorted(a2g_editing_sites_df.at[record.id,a2g_editing_sites_df.columns[0]]):
#                 if prot_end_nuc>=k>=prot_start_nuc:
#                     a2g_list.append(prot_end_nuc-k)
#                 else:
#                     out_of_frame_a2g.append(k)
#         except:
#             pass
#             
#         #c2t editing site for sequence
#         try:
#             for k in sorted(c2t_editing_sites_df.at[record.id,c2t_editing_sites_df.columns[0]]):
#                 if prot_end_nuc>=k>=prot_start_nuc:
#                     c2t_list.append(prot_end_nuc-k)
#                 else:
#                     out_of_frame_c2t.append(k)
#         except:
#             pass
#         
# 
#     if bad_seq:
#         n_bad+=1
#         current_record = SeqRecord(orf_sequence, id = record.id, 
#                                description = '| a2g: ' + str(a2g_list) + ' | c2t: ' + str(c2t_list) + 
#                                ' | prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + 
#                                ' | strand: ' + strand + ' | prot_srat_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc) + ' | original_orf_start: ' + str(orf_start) + ' | original_orf_end: ' + str(orf_end))    
#         writer_bad.write_record(current_record)
#         
#     else:
#         n_good+=1
#         current_record = SeqRecord(final_sequence, id = record.id, 
#                                description = '| a2g: ' + str(a2g_list) + ' | c2t: ' + str(c2t_list) + 
#                                ' | prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + 
#                                ' | strand: ' + strand + ' | prot_srat_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc) + ' | original_orf_start: ' + str(orf_start) + ' | original_orf_end: ' + str(orf_end))
#         
#         writer.write_record(current_record)
#         #write out of frame sites to designated file
#         if len(out_of_frame_c2t) + len(out_of_frame_a2g):
#             out_of_range_sites_file.write(record.id + ' |protein_frame (base1): ' + str(prot_start_nuc)+'-'+str(prot_end_nuc) + ' | strand: ' + strand + ' | a2g: ' + str(out_of_frame_a2g) + ' | c2t: ' + str(out_of_frame_c2t) + '\n')
# 
#             
# #        test = record.id +' | len: '+ str(len(final_sequence)) + ' | prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + ' | strand: ' + strand + ' | prot_start_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc)
# #        print(test)
#     if len(final_sequence)%3 or prot_end_nuc-prot_start_nuc+1 != len(final_sequence):# or bad_orf:# record.id == 'comp134935_c0_seq6':
#         break
#     if len(orf_sequence)%3:
#         break
#         
# out_of_range_sites_file.close()       
# writer.write_footer()
# writer_bad
# print(str(n_good) + ' sequences written as protein sense strand')
# print(str(n_bad) + ' bad sequences written to designated file')
# =============================================================================




