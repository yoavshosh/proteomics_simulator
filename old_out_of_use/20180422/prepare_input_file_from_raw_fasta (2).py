import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_rna



seq = 'AACATGAAAGGTGACATGTTTTGTAGAAAGAAAAAACATAAAGAATGGCTAAAAATGGCTAAAGCTGGCCAATACACCTTAAAAGAGAAAATTCAAGAAGAGTTTTTGGAATGCAAAATCTGTTTTGAACCGTATGTAAAACCGAAGGCATTACCTTGCCTTCACTCATTCTGTGCCGAGTGTCTAAAGGACTACGTGCGAAAGAACCCCAACAAAAATGCAGTGCGTTTTTGCTGTCCAATCTGTCGCAAAGAAATCCCGATGCCTGCCGGTGGCATCGACGATTTCCAAGATAATTTTTGGTTATTGAGCTTGTCAAATTCCTTGGAAGAAGGAGAGGAAGACTGTCGTGTATCGTGCAATGGGAAGGCCGTGCGTAGCAATGCCTGGCCTACTCCAAAATGGAAACAACCCAAGAATGTTTCAACCTCATCTAAACCTTTGTATCCCCCATCATTCCAGGAACTTAAACTCAGGGGGTTTGAATGGTATTTTGGTAAAGTGAGTCGTAATGCTTCAGAAGAATGGCTTCTCCATCCTGGGCTACAGAAGGGAACATTCCTTATTCGACAAGGGGAGGCTCTTCCTGACACGTATACCTTGTCAGTTCGTGACTGTGACGAGCTGAGGGGTTACCTGGTGAAACATTACAAGATTCTAACCAAAAAAGCAAGTGATGGGGAGAAGGAAGTTTATTACATCACCCCAAAGCGAACTTTCCGTTCGCTTGAGGAACTGGTGAATCATTATTCAATTTCAGACGGCCTTTGCTGTAAACTAACACAAATATGTAACAAACCAAGGTCTTTACTCTGGGCAATGGAAAGAGGAAAACCGGATGATTTTATGACCACTAAGGACACTTTGCAATTGGTCAAGAAAATTGGCAGTGGCCAATTTGCTGAAGTTTATTATGCCAAATGGAACAACCAAGTAGAGGCAGCCGTAAAGATGCAAAAAAAGGATTGCGTGACAACCTCGGCTTTCTTAGATGAAGCTCAGATCTTAAAGACAATTCAACATGTAAATATAATCAAGTTGTTGGCTGTGTGCAGTGACGAACCTGTTTATTTGGTGACAGAATATATGCCTAATGGCCGACTCTCACAATACCTTCGAGAAGGAAAGGGGAAACAGCTCGGGGTTAACAGTCTCCTCTGGCTTGCGGCTCAGATCGCAGACGGAATGGCTTACATGGAAAAGGAGAATTTTGTTCATAGAAATTTGGGTGCCCGAAACATTCTGGTTGCTGACCAAAACAAAGTGAAAATTGCTGGTTTTGGGATGACAAAAGTGGCCGATGATCCTGATTACAATTTCAGAAAAGGTTTGAAAATGGCTGTAAAATGGATGGCCCCTGAAGTGTTGTTGTACAATAAATACAGCACAAAGGCTGACGTGTGGTCGTTTGGAATTGTCCTAATGGAAATATTTTCATATGGCAAAGAACCCTACGATGGTATGGGAAGCAAGGAAGCATTTGAAAATGTCCAGTCAGGTTACCGAATGCCTTGTCCTCACTGTTGCCCTGCCGAGGTATACAATGTGGCACTGACCTGTTGGAATATCAACACCCAGCGTCGGCCATCTTTTGACTTTCTTAACAGCTTCCTTCATGACTGGCACTATACGTCC'
i_path = r'C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/test_files/'
o_path = i_path
fasta_input = 'orfs_squ_edited_plusminus_whole_comp_ReverseComplimentMinustrand.fasta'
xls_path = r'C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/test_files/'
xls_file = 'editing_sites.xlsx'
xls_columns_list = [0,3]
xls_sheet = 'Editing_events_in_Squid'

strand_orientation_regex = re.compile('(?<=Strand\s)[^\s]+')
start_nuc_regex = re.compile('(?<=OrfStart\s)[^\s]+')
end_nuc_regex = re.compile('(?<=OrfEnd\s)[^\s]+')


"""
import editng sites from xls to a pd framework
"""
#df = import_es_from_xls(xls_path,xls_file,xls_sheets_list = 0, xls_columns_list = xls_columns_list)
def import_es_from_xls(xls_path,xls_file,xls_sheet = 0, xls_columns_list = xls_columns_list):
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
    

 
def reading_frame_rna(sequence,orientation,start_nuc,end_nuc):
    
    in_farme_seq = Seq(sequence[start_nuc-1:end_nuc],generic_rna)

    if orientation == '+':
        final_rna = in_farme_seq
    elif orientation == '-':
        final_rna = in_farme_seq.reverse_complement()
    else:
        final_rna = None
        
    return final_rna
    


def create_in_frame_rna_file(o_path,i_path,fasta_input,atg_editing_sites_df,c2t_editing_sites_df = []):
    
    writer =  FastaWriter(open(o_path + 'in_frame_rna_from_' + fasta_input , 'w'), wrap=None)
    writer.write_header()
    
    for record in SeqIO.parse(open(i_path + fasta_input, "r"), "fasta"):
        orientation = find_by_regex_in_header(record.description, strand_orientation_regex) #'+' if seq is read as is, '-' if reverse compliment is read
        start_nuc = int(find_by_regex_in_header(record.description, start_nuc_regex))
        end_nuc = int(find_by_regex_in_header(record.description, end_nuc_regex))

        sequence = reading_frame_rna(str(record.seq),orientation,start_nuc,end_nuc)
            
        #a2g editing site for sequence
        try:
            a2g_list = [k - start_nuc for k in sorted(atg_editing_sites_df.at[record.id,atg_editing_sites_df.columns[0]])]
        except:
            a2g_list = []
            
        #c2t editing site for sequence
        try:
            c2t_list = [k - start_nuc for k in sorted(atg_editing_sites_df.at[record.id,c2t_editing_sites_df.columns[0]])]
        except:
            c2t_list = []
            
        writer.write_record(SeqRecord(sequence, id = record.id, 
                                      description = '| a2g: ' + str(a2g_list) + ' | c2t: ' + str(c2t_list)))
            
    writer.write_footer()

#atg_editing_sites_df = import_es_from_xls(xls_path,xls_file,xls_sheet = xls_sheet, xls_columns_list = xls_columns_list)
#create_in_frame_rna_file(o_path,i_path,fasta_input,atg_editing_sites_df)

            