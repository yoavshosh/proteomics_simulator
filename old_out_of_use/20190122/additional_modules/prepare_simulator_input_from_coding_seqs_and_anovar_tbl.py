import re
import os
import sys
import pandas as pd
import numbers
import operator
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis


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

def create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs_final_sites_dfs(fasta_path,fasta_file,out_path,anovar_df_a2g,anovar_df_c2t=None,variants_to_use = []):
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
    
    writer =  FastaWriter(open(out_path + 'in_frame_rna_from_' + fasta_file , 'w'), wrap=None)
    writer_bad = FastaWriter(open(out_path + 'bad_seq_from_' + fasta_file , 'w'), wrap=None)
    writer.write_header()
    writer_bad.write_header()
    
    for record in SeqIO.parse(open(fasta_path + fasta_file, "r"), "fasta"):
        
        rec_id = record.id.split(';')[0]
        use_variant = True
        
        if len(variants_to_use): #if a not-empty list is passed for variants_to_use, flag variants that are not in list so they will not be included in uotput
            if rec_id not in variants_to_use:
                use_variant = False
        
        if use_variant:
            final_sequence = str(record.seq[0:-3]).replace('a','A').replace('g','G').replace('t','T').replace('c','C')
            try:
                a2g_list = [int(k)-1 for k in anovar_df_a2g.loc[[rec_id]]['position_base1']]
            except KeyError:
                a2g_list = []
                
            try:
                c2t_list = [int(k)-1 for k in anovar_df_c2t.loc[[rec_id]]['position_base1']]
            except KeyError:
                c2t_list = []   
            
            prot_start_nuc = 1
            prot_end_nuc = len(final_sequence) - 3
            prot_start = 'first_met_in_original_orf'
            prot_end = 'original_sense_strand_orf_end'
            strand = '+'
            orf_start = 1
            orf_end = len(record.seq) - 3
            description_str = '| a2g_base0: ' + str(a2g_list) + ' | c2t_base0: ' + str(c2t_list) + '| prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + ' | strand: ' + strand + ' | prot_start_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc) + ' | original_orf_start: ' + str(orf_start) + ' | original_orf_end: ' + str(orf_end)
        
            if (Seq(str(record.seq[0:3]), generic_dna).translate() != 'M' or Seq(str(record.seq[-3:len(record.seq)]), generic_dna).translate() != '*') or '*' in Seq(str(final_sequence), generic_dna).translate():
                writer_bad.write_record(record)
                n_bad+=1
                sites_bad+=len(a2g_list)+len(c2t_list)
            else:
                if len(final_sequence)%3:
                    final_sequence=final_sequence[0:-len(final_sequence)%3]
                current_record = SeqRecord(Seq(final_sequence,generic_dna), id = rec_id, description = description_str)
                writer.write_record(current_record)
                n_good+=1
                sites_good+=len(a2g_list)+len(c2t_list)
            
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



#==========================================================================================================
#Yeast Functions   
    
    

"""
read dataframe from anovar output
"""
def read_anovar_yeast_data(input_path, input_file, only_recoding = True):
     
    anovar_gene_change_regex = re.compile('(.+?)\:(.+)\:(.+).?\:c.(.+?)\:p.(.+)')
    
    columns = ['variant_name', 'mm_type', 'mrna_mm_position','nuc_mm','prot_mm_position','aa_mm','genomic_coordinates']
    data = []
    
    with open(input_path + input_file, "r") as f:
        content = f.readlines()
        for line in content:
            fields = line.split("\t")
            mm_type = fields[1]
            variant_name_and_changes = fields[2].split(',')
            genomic_coordinates = fields[3]+'_'+fields[4]+'_'+fields[8]
            for mm in variant_name_and_changes:
                if mm != '':
                    mm_data = find_by_regex_in_str(mm,anovar_gene_change_regex)
                    variant_name = mm_data[1]
                    mrna_mm_position = int(mm_data[3][1:-1])
                    nuc_mm = mm_data[3][0]+mm_data[3][-1]
                    prot_mm_position = int(mm_data[4][1:-1])
                    aa_mm = mm_data[4][0]+mm_data[4][-1]
                    data.append([variant_name, mm_type, mrna_mm_position, nuc_mm, prot_mm_position, aa_mm, genomic_coordinates])
    
    anovar_df = pd.DataFrame(data = data, columns = columns)
    if only_recoding:
        anovar_df = anovar_df[anovar_df['mm_type'] == 'nonsynonymous SNV']
    
    
    anovar_df_a2g = anovar_df[anovar_df['nuc_mm'] == 'AG']
    anovar_df_a2g.set_index('variant_name', inplace = True)
    anovar_df_a2g = anovar_df_a2g.groupby(anovar_df_a2g.index).agg({list(anovar_df_a2g.columns)[1]:lambda x: list(set(list(x)))})

    anovar_df_c2t = anovar_df[anovar_df['nuc_mm'] == 'CT']
    anovar_df_c2t.set_index('variant_name', inplace = True)
    anovar_df_c2t = anovar_df_c2t.groupby(anovar_df_c2t.index).agg({list(anovar_df_c2t.columns)[1]:lambda x: list(set(list(x)))})
        
    return anovar_df, anovar_df_a2g, anovar_df_c2t
 

def create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs(fasta_path,fasta_file,out_path,anovar_df_a2g,anovar_df_c2t,only_recoding=True):
    
    if only_recoding:
        out_str = 'rec_only'
    else:
        out_str = 'all_sites'
    
    writer =  FastaWriter(open(out_path + 'in_frame_rna_' + out_str + '_from_' + fasta_file , 'w'), wrap=None)
    writer_bad = FastaWriter(open(out_path + 'bad_seq_from_' + fasta_file , 'w'), wrap=None)
    writer.write_header()
    writer_bad.write_header()
    for record in SeqIO.parse(open(fasta_path + fasta_file, "r"), "fasta"):
     
        final_sequence = record.seq[0:-3]
        try:
            a2g_list = [k-1 for k in anovar_df_a2g.loc[record.id][0]]
        except KeyError:
            a2g_list = []
        try:
            c2t_list = [k-1 for k in anovar_df_c2t.loc[record.id][0]]
        except KeyError:
            c2t_list = []  
            
        prot_start_nuc = 1
        prot_end_nuc = len(final_sequence) - 3
        prot_start = 'first_met_in_original_orf'
        prot_end = 'original_sense_strand_orf_end'
        strand = '+'
        orf_start = 1
        orf_end = len(record.seq) - 3
        description_str = '| a2g: ' + str(a2g_list) + ' | c2t: ' + str(c2t_list) + '| prot_start: ' + str(prot_start) + ' | prot_end: ' + str(prot_end) + ' | strand: ' + strand + ' | prot_start_nuc: ' + str(prot_start_nuc) + ' | prot_end_nuc: ' + str(prot_end_nuc) + ' | original_orf_start: ' + str(orf_start) + ' | original_orf_end: ' + str(orf_end)
    
        if (Seq(str(record.seq[0:3]), generic_dna).translate() != 'M' and Seq(str(record.seq[-3:len(record.seq)]), generic_dna).translate() != '*') or '*' in Seq(str(final_sequence), generic_dna).translate():
            writer_bad.write_record(record)
        else:
            current_record = SeqRecord(final_sequence, id = record.id, description = description_str)
            writer.write_record(current_record)
            
    writer_bad.write_footer()
    writer.write_footer()
    
#==========================================================================================================

    
if __name__ == '__main__':
        
#    input_path = 'C:/Users/user/Google_Drive/RNA_Editing/files/Yeast/results_from_yishai/'
#    fasta_file = 'SC_refGeneMrna.fa'
#    anovar_file = 'RES_final_result_editLvl_ann.exonic_variant_function'
#    only_recoding = True


    animal = sys.argv[1]
    
    assert(animal in ['human','yeast']), animal + ' not a recognized animal'
    
    if animal == 'human':
        input_path = sys.argv[2]
        if not os.path.exists(input_path + 'proteomics_simulation/'):
            os.makedirs(input_path + 'proteomics_simulation/')
        out_path = input_path + 'proteomics_simulation/'     
        
        fasta_file = sys.argv[3]
        anovar_file = sys.argv[4]
        
        print('Reading and processing Anovar tabel')
        sites = read_editing_sites_wrt_coding_seqs(anovar_file)
        anovar_df_a2g = sites[sites['mm_type']=='AG']
        anovar_df_a2g.set_index('ucsc_id', inplace = True)
        
        anovar_df_c2t = sites[sites['mm_type']=='CT']
        anovar_df_c2t.set_index('ucsc_id', inplace = True)
        
        
        use_longest = eval(sys.argv[5])
        if use_longest:
            refGene_file = sys.argv[6]
            print('Preparing list of longest variants')
            variants_df, longest_variants = get_genes_variants_length(refGene_file)
            print('Writing fasta file in the proteomics simulator format')
            create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs_final_sites_dfs(input_path,fasta_file,out_path,anovar_df_a2g, anovar_df_c2t = anovar_df_c2t,variants_to_use=longest_variants)            
        else:
            print('Writing fasta file in the proteomics simulator format')
            create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs_final_sites_dfs(input_path,fasta_file,out_path,anovar_df_a2g, anovar_df_c2t = anovar_df_c2t)
        
    elif animal == 'yeast':
        input_path = sys.argv[2]
        fasta_file = sys.argv[3]
        anovar_file = sys.argv[4]
        only_recoding = True  
        
        if not os.path.exists(input_path + 'proteomics_simulation/'):
            os.makedirs(input_path + 'proteomics_simulation/')
        out_path = input_path + 'proteomics_simulation/'
        
        print('Reading and processing Anovar tabel')
        anovar_df, anovar_df_a2g, anovar_df_c2t = read_anovar_data(input_path,anovar_file,only_recoding = only_recoding)
        print('Writing fasta file in the proteomics simulator format')
        create_in_frame_rna_file_from_anovar_results_and_coding_mrna_seqs(input_path,fasta_file,out_path,anovar_df_a2g,anovar_df_c2t,only_recoding = only_recoding)
