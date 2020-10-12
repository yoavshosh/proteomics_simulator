import re
import sys
import pandas as pd
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


"""
import editng sites from xls to a pd framework
"""
#df = import_es_from_xls(xls_path,xls_file,xls_sheet = xls_sheet, xls_columns_list = xls_columns_list, only_recoding = True)
def import_es_from_xls(xls_path,xls_file,xls_sheet = 0, xls_columns_list = xls_columns_list, only_recoding = False):
    if only_recoding:
        df = pd.read_excel(xls_path + xls_file, index_col= 0, sheetname = xls_sheet, parse_cols = xls_columns_list)
        df = df.loc[df['Codon changes'] != 'syn']
        df = df.groupby(df.index).agg({list(df.columns)[0]:lambda x: list(x)})
    else:
        df = pd.read_excel(xls_path + xls_file, index_col=0, sheetname = xls_sheet, parse_cols = xls_columns_list)
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
for a given sequence if beginning and end are not stop codons - than real orfs are actually undetermined
elsewise, determine orf from first methionine until stop stop codon
return new sequence, first methionon location, and markers for frame beginning and end (unknown if stop codons werent found)
"""
def rethink_reading_frame(sequence):
    
    first_met = 'no_met'
    first_aa = Seq(str(sequence[0:3]), generic_rna).translate()
    last_aa = Seq(str(sequence[-3:]), generic_rna).translate()

    if str(first_aa) != '*':
        frame_beginning = 'unknown'
    else:
        first_met = find_first_met(sequence)
        frame_beginning = 'within_original_orf'
        sequence = sequence[first_met:]
    
    if str(last_aa) == '*':
        print('end')
        sequence = sequence[:-3]
        frame_end = 'seq_end'
    else:
        frame_end = 'unknown'
        
    return sequence, first_met, frame_beginning, frame_end



def reading_frame_rna(sequence,strand,start_nuc,end_nuc):
    
    in_farme_seq = Seq(str(sequence[start_nuc-1:end_nuc]),generic_dna)

    if strand == '+':
        final_rna = in_farme_seq
    elif strand == '-':
        final_rna = in_farme_seq.reverse_complement()
    else:
        final_rna = None
        
    return final_rna
    

def get_aa_after_orf(seq,strand,start_nuc,end_nuc):
    if strand == '+':
        codon = seq[end_nuc:end_nuc+3]
#        print('+' + codon)
    elif strand == '-':
        codon = seq[start_nuc-4:start_nuc-1].reverse_complement()
#        print('-' + codon)
    aa = Seq(str(codon), generic_rna).translate()
    
    return aa
    
        

def check_codon_after_frame(input_path, input_fasta,strand_regex,start_nuc_regex,end_nuc_regex):
    
    aa_after_orf = {}
    
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        strand = find_by_regex_in_header(record.description, strand_regex)
        start_nuc = int(find_by_regex_in_header(record.description, start_nuc_regex))
        end_nuc = int(find_by_regex_in_header(record.description, end_nuc_regex))
        aa = get_aa_after_orf(record.seq,strand,start_nuc,end_nuc)
        aa_after_orf.update({record.id:str(aa)})
        
    return aa_after_orf
        


def create_in_frame_rna_file(out_path,input_path,fasta_input,a2g_editing_sites_df,c2t_editing_sites_df = [], only_recoding = False):
    
    if only_recoding:
        out_str = 'rec_only'
    else:
        out_str = 'all_sites'
     
    out_of_range_sites_file = open(input_path + 'out_of_frame_sites_' + out_str + '_from_' + fasta_input + '.txt', "w")
    writer =  FastaWriter(open(out_path + 'in_frame_rna_' + out_str + '_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    for record in SeqIO.parse(open(input_path + fasta_input, "r"), "fasta"):
        
        a2g_list = []
        c2t_list = []
        out_of_frame_a2g = []
        out_of_frame_c2t = []
        
        strand = find_by_regex_in_header(record.description, strand_regex) #'+' if seq is read as is, '-' if reverse compliment is read
        
        #start_nuc and end_nuc are nucleotides locations within original sequence (regardless of strand) in input fasta stating the orfs 
        start_nuc = int(find_by_regex_in_header(record.description, start_nuc_regex))
        end_nuc = int(find_by_regex_in_header(record.description, end_nuc_regex))

        #get initial reading frame by the orfs stated in input fasta (input fasta orfs are blastx orfs prolonged until stop codons or until edges of aligned query sequences)
        sequence = reading_frame_rna(record.seq,strand,start_nuc,end_nuc)
        
        #if stop cidins at seq edges, shorten reading frame and return new seq, first methionin location, and markers for wether or not seq edges are real or not
        sequence, first_met, frame_beginning, frame_end = rethink_reading_frame(sequence)
    
        #stop codon at beginning of sequence: sequence -> sequence from first methionin and updating in sense strand start_nuc
        if frame_beginning == 'within_original_orf':
            if strand == '+':
                start_nuc = start_nuc+first_met
            elif strand == '-':
                end_nuc = end_nuc-first_met
        
        #stop codon at end of sequence: sequence -> sequence untill stop codond and updating sense strand end_nuc
        if frame_end == 'seq_end':
            if strand == '+':
                end_nuc = end_nuc-3
            elif strand == '-':
                start_nuc = start_nuc+3

        #a2g editing site for sequence
        if strand == '+':
            try: #taking all in frame editing sites from sites dataframe
                for k in sorted(a2g_editing_sites_df.at[record.id,a2g_editing_sites_df.columns[0]]):
                    if end_nuc>=k>=start_nuc:
                        a2g_list.append(k-start_nuc)
                    else: 
                        out_of_frame_a2g.append(k-start_nuc)
            except:
                pass
            
            #c2t editing site for sequence
            try:
                for k in sorted(c2t_editing_sites_df.at[record.id,c2t_editing_sites_df.columns[0]]):
                    if end_nuc>=k>=start_nuc:
                        c2t_list.append(k-start_nuc)
                    else:
                        out_of_frame_c2t.append(k-start_nuc)
            except:
                pass
            
        
        elif strand == '-':
            try:
                for k in sorted(a2g_editing_sites_df.at[record.id,a2g_editing_sites_df.columns[0]]):
                    if end_nuc>=k>=start_nuc:
                        a2g_list.append(end_nuc-k)
                    else:
                        out_of_frame_a2g.append(end_nuc-k)
            except:
                pass
                
            #c2t editing site for sequence
            try:
                for k in sorted(c2t_editing_sites_df.at[record.id,c2t_editing_sites_df.columns[0]]):
                    if end_nuc>=k>=start_nuc:
                        c2t_list.append(end_nuc-k)
                    else:
                        out_of_frame_c2t.append(end_nuc-k)
            except:
                pass
            
        
        #write out of frame sites to designated file
        if len(out_of_frame_c2t) + len(out_of_frame_a2g):
            out_of_range_sites_file.write(record.id + ' |reading frame (base0): ' + str(start_nuc)+'-'+str(end_nuc) + ' | strand: ' + strand + ' | a2g: ' + str(out_of_frame_a2g) + ' | c2t: ' + str(out_of_frame_c2t))
 
        writer.write_record(SeqRecord(sequence, id = record.id, 
                                      description = '| a2g: ' + str(a2g_list) + ' | c2t: ' + str(c2t_list) + 
                                      ' | prot_beginning: ' + frame_beginning + ' | prot_end: ' + frame_end + 
                                      ' | strand: ' + strand + ' | start_nuc: ' + str(start_nuc) + ' | end_nuc: ' + str(end_nuc)))
    
    out_of_range_sites_file.close()       
    writer.write_footer()
    
#
#atg_editing_sites_df = import_es_from_xls(xls_path,xls_file,xls_sheet = xls_sheet, xls_columns_list = xls_columns_list)
#create_in_frame_rna_file(o_path,i_path,fasta_input,atg_editing_sites_df)            
#
#
#atg_editing_sites_df = import_es_from_xls(xls_path,xls_file,xls_sheet = xls_sheet, xls_columns_list = xls_columns_list, only_recoding = True)
#create_in_frame_rna_file(o_path,i_path,fasta_input,atg_editing_sites_df,only_recoding = True)            
                        
            
# =============================================================================
if __name__ == "__main__":
     
    i_path = sys.argv[1]
    xls_file = sys.argv[2]
    xls_sheet = sys.argv[3]
    only_recoding = sys.argv[4]
    fasta_input = sys.argv[5]
    o_path = i_path
    
    xls_columns_list = [0,3,4]
    strand_regex = re.compile('(?<=Strand\s)[^\s]+')
    start_nuc_regex = re.compile('(?<=OrfStart\s)[^\s]+')
    end_nuc_regex = re.compile('(?<=OrfEnd\s)[^\s]+')
     
    print('reading editing sites from Excel')
    atg_editing_sites_df = import_es_from_xls(i_path,xls_file,xls_sheet = xls_sheet, xls_columns_list = xls_columns_list)
    print('creating in-frame fasta from input fasta file (containing Editing sites list relative to frame - base0')
    create_in_frame_rna_file(o_path,i_path,fasta_input,atg_editing_sites_df,only_recoding = only_recoding)          
# =============================================================================

