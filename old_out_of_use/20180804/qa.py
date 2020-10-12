import re
import pandas as pd
import sys
#from in_frame_rna_fasta import import_es_from_xls, find_by_regex_in_header
from edited_peptides_from_seq import create_all_cleavage_sites, find_sites_in_window
from Bio import SeqIO 

#input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
#input_fasta = 'in_frame_rna_rec_only_from_orfs_squ.fasta'
#digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
#a2g_header_regex = re.compile(r'(?<=a2g:\s).*?]')
#c2t_header_regex = re.compile(r'(?<=c2t:\s).*?]')
#strand_orientation_regex = re.compile('(?<=Strand\s)[^\s]+')
#input_fasta = 'orfs_squ.fasta'

#input_fasta = 'in_frame_rna_recoding_sites_only_from_orfs_squ.fasta'
#df_sus = editing_sites_per_pep(input_path, input_fasta, a2g_header_regex, c2t_header_regex, digestion_rule)
#input_fasta = 'in_frame_rna_all_sites_from_orfs_squ.fasta'
#df_all = editing_sites_per_pep(input_path, input_fasta, a2g_header_regex, c2t_header_regex, digestion_rule)

#sites = editing_sites_per_pep(input_path, input_fasta, a2g_header_regex, c2t_header_regex, digestion_rule)

def editing_sites_per_pep(input_path, input_fasta, a2g_header_regex, c2t_header_regex,
                          digestion_rule, missed_cleavages = 0, 
                          look_ahead_for_editnig = 3, look_behind_for_editing = 6):
    n=0
    peps_df_list = []
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        n+=1
        print('record number ' + str(n))
#        print(str(record.seq))
        
        seq_length = len(str(record.seq))
        a2g_sites = sorted(eval(find_by_regex_in_header(record.description,a2g_header_regex)))
        c2t_sites = sorted(eval(find_by_regex_in_header(record.description,c2t_header_regex)))

        cleavage_sites, binary_for_fixed_cs = create_all_cleavage_sites(str(record.seq),digestion_rule,a2g_sites,c2t_sites,
                                                                        look_ahead_for_editnig,look_behind_for_editing)
        cleavage_sites_lengeth = len(cleavage_sites)
#        print(cleavage_sites)
#        print(binary_for_fixed_cs)
        
        for i in range(cleavage_sites_lengeth-1):
            
#            print(i)
            
            fixed_sites_cnt = 0
            j=1
            
            while (fixed_sites_cnt < missed_cleavages+1) and (cleavage_sites_lengeth-i-1>=j):
                
#                print([i,j])
            
                if binary_for_fixed_cs[i+j]: #count site if site ahead is fixed.
                    fixed_sites_cnt+=1
            
                #determining range for editing perutations:
                if binary_for_fixed_cs[i]: #site is fixed and therfore upstream editing dosent matter
                    edit_range_start = cleavage_sites[i]
                else:    
                    if cleavage_sites[i]-look_behind_for_editing < look_behind_for_editing: #sequence is close to beginning, looking for permutations from beginning
                        edit_range_start = 0
                    else:
                        edit_range_start = cleavage_sites[i]-look_behind_for_editing
            
                if binary_for_fixed_cs[i+j]: #site is fixed and therfore downstream editing dosent matter
                    edit_range_end = cleavage_sites[i+j]
                else:    
#                    print([i,j])
                    if cleavage_sites[i+j]+look_ahead_for_editnig > seq_length: #sequence is close to end, looking for permutations till end
                        edit_range_end = seq_length
                    else:
                        edit_range_end = cleavage_sites[i+j]+look_ahead_for_editnig
#                        print(str(edit_range_end))
                        
                perm_range = (edit_range_start,edit_range_end-1)
                perm_range_str = str(perm_range[0]) + '_' + str(perm_range[1])
                
                a2g_sites_in_range = find_sites_in_window(perm_range,a2g_sites)
                c2t_sites_in_range = find_sites_in_window(perm_range,c2t_sites)
                
                pep_df = pd.DataFrame({'compenent_id':[record.id],'pep_in_frame_extended_coor':[perm_range_str],'num_editing_sites':[len(a2g_sites_in_range) + len(c2t_sites_in_range)]})
                
                peps_df_list.append(pep_df)
                
                j+=1
    
    peps_df = pd.concat(peps_df_list)
    del(peps_df_list)
    peps_df.set_index(['compenent_id','pep_in_frame_extended_coor'], inplace = True)
    peps_df.to_pickle(input_path + 'editing_site_per_pep_from_' + input_fasta + '.pickle')
    peps_df.to_excel(input_path + 'editing_site_per_pep_from_' + input_fasta + '.xlsx')
    
    return peps_df
    
    

"""
return 2 lists of 'A' indeces and 'C' indeces' in sequence
"""
#a_list, c_list = c_a_list(sequence)
def c_a_list(sequence):
    
    a_list = []
    c_list = []
    
    for i in range(len(sequence)):
        if sequence[i:i+1] == 'A':
            a_list.append(i)
        elif sequence[i:i+1] == 'C':
            c_list.append(i)
            
    return a_list, c_list


def es_per_pep(es,cs):
    
    cs_per_pep = []
    for i in range(len(cs)-1):
        n = sum(1 for k in es if k>=cs[i] and k<=cs[i+1])
        if n == 19:
            print(str(i))
        cs_per_pep.append(n)
    return cs_per_pep



def get_mean_editing_sites(input_path,input_fasta,a2g_header_regex,c2t_header_regex):
    total_editing_sites = 0
    rec_num = 0
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        rec_num+=1
        a2g_sites = eval(find_by_regex_in_header(record.description,a2g_header_regex))
        c2t_sites = eval(find_by_regex_in_header(record.description,c2t_header_regex))
        editing_sits = len(a2g_sites)+len(c2t_sites)
        total_editing_sites+=editing_sits
    return total_editing_sites/rec_num



def get_mean_cleavage_sites(input_path,input_fasta,digestion_rule,a2g_header_regex,c2t_header_regex,look_ahead_for_editnig,look_behind_for_editing): 
    total_cleavage_sites = 0
    rec_num = 0
    for record in SeqIO.parse(input_path + input_fasta, "fasta"):
        rec_num+=1
        a2g_sites = eval(find_by_regex_in_header(record.description,a2g_header_regex))
        c2t_sites = eval(find_by_regex_in_header(record.description,c2t_header_regex))
        c,b = create_all_cleavage_sites(str(record.seq),digestion_rule,a2g_sites,c2t_sites,look_ahead_for_editnig,look_behind_for_editing)
        cleavage_sites = len(c)
        total_cleavage_sites+=cleavage_sites
    return total_cleavage_sites/rec_num
    
    

def editing_above_threshols(input_path,input_fasta,editing_threshod):
    n=0
    for record in SeqIO.parse(input_path + input_fasta, "fasta"): 
        
        a2g_sites = eval(find_by_regex_in_header(record.description,a2g_header_regex))
        c2t_sites = eval(find_by_regex_in_header(record.description,c2t_header_regex))
        if len(a2g_sites) + len(c2t_sites) > editing_threshod:
            print(record.id)
            n+=1            
    return n
    


#recs = getrecords('C:/Users/user/Google_Drive/RNA_Editing/yeast_proteomics/test_files/','in_frame_rna_from_orfs_squ.fasta',172)     
"""
get a list of record objects in requestes rows on fasta file
"""
def getrecords(input_path,input_fasta,*rec_nums):
    rec_list =[]
    n=0
    for record in SeqIO.parse(open(input_path + input_fasta, "r"), "fasta"):
        n+=1
        if n in rec_nums:
            rec_list.append(record)
            
    return rec_list
        

def get_rec_num(input_path,input_fasta):
    n=0
    for record in SeqIO.parse(open(input_path + input_fasta, "r"), "fasta"):
        n+=1
    return n



def get_comps_with_orientation(input_path, input_fasta, orientation, strand_orientation_regex):
    
    ids_list = []
    
    for record in SeqIO.parse(open(input_path + input_fasta, "r"), "fasta"):
        if find_by_regex_in_header(record.description, strand_orientation_regex) == orientation:
            ids_list.append(record.id)
            
    return ids_list
       

"""
find regex (regex - pre compiled regex) in header
"""
def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'
  


if __name__ == '__main__':
    
    digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
    a2g_header_regex = re.compile(r'(?<=a2g:\s).*?]')
    c2t_header_regex = re.compile(r'(?<=c2t:\s).*?]')
    input_path = sys.argv[1]
    input_fasta = sys.argv[2]
    
    sites = editing_sites_per_pep(input_path, input_fasta, a2g_header_regex, c2t_header_regex, digestion_rule)

