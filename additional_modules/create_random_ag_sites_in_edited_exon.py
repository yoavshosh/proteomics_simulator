# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 11:54:56 2019

@author: shosh

This scripts recives all transcriptome exons data and all sites data
for each edited exon it tries to find other similar codons within the exon that can be set as new random editing sites.
if there are not enough unedited codons that are similar to the edited ones, it seaches for codons in other exon of the same gene
"""
import pandas as pd
import numpy as np
import random
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



def flatten_sites_positions_in_exon(k,l):
    """
    This function unpack nested list of sites in exons
    (if two sites are in the same codon they will appear in l as a tuple)
    """
    all_sites_exon_position = []
    for s in l:
        if type(s) == tuple:
#            print(k)
            all_sites_exon_position = all_sites_exon_position + [k for k in s]
        else:
            all_sites_exon_position.append(s)
    return all_sites_exon_position


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
         
    columns = ['variant_name', 'position_base0', 'position_base1', 'mm_type', 'chromosome', 'genomic_position_base0', 'genomic_position_base1', 'strand', 'aa_change', 'editin_level', 'genomic_key_base1','coding_key_base1']
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


def find_all_edited_codons_in_exon(exon_row,exons_edited_codons_dict, genomic_keys_for_attended_sites):
    """
    All positions in this function are base0
    This funcition takes an exon_row from edited exons dataframe and retuens a nested dictionary of exons.
    each value is itself a dictionary of:
    1. a list of edited codons - each codon in list is a tuple of codon sequence and and all nucleotides in codon that are editing sites
    2. a list of position of the editing sites relative to the exon start (elements in this list are corresponding to elements in the codons list. 
       if some codon has more than one site in it, the corresponding value in that list will be a tuple of the editing sites positions relative to the exon start.
    3. corresponding codons positions in exon
    """
    
    def find_site_codon(exon, reading_frame, site_ex_position):
        """
        Given an exon sequence, its reading frame (0/1/2) and site position within the exon
        Returns the edited codon and the site position within the codon (0/1/2)
        """
        nucs_before_site = (site_ex_position-reading_frame)%3
        nucs_after_site = ((-site_ex_position-1)%3 - (-reading_frame)%3)%3
        codon = exon[site_ex_position-nucs_before_site:site_ex_position+nucs_after_site+1]
        return (codon,nucs_before_site), site_ex_position-nucs_before_site  #return codons and site location in codon (which in base0 eq to nucs before site) as tuple and codon location
    
    
    def merge_sites_in_same_codon(codons,reading_frame,sites_exon_positions):
        """
        This function merge elements in codons list if two elements are different site in the same codon
        """
        new_codons_list = [codons[0]]
        new_sites_exon_positions_list = [sites_exon_positions[0]]
        
        for i in range(1,len(sites_exon_positions)):
            #check if site souled be merged with precedent site
            if sites_exon_positions[i]-sites_exon_positions[i-1]<3 and (sites_exon_positions[i]-reading_frame)%3>(sites_exon_positions[i-1]-reading_frame)%3:
                new_codons_list[-1] = new_codons_list[-1] + (codons[i][1],)
                new_sites_exon_positions_list[-1] = (new_sites_exon_positions_list[-1],) + (sites_exon_positions[i],) 
            else:
                new_codons_list.append(codons[i])
                new_sites_exon_positions_list.append(sites_exon_positions[i])
        return new_codons_list, new_sites_exon_positions_list
    
    
    exon_key = exon_row['exon_key']
    exon_seq = exon_row['exon'].replace('a','A').replace('g','G').replace('c','C').replace('t','T')
    exon_coding_start = exon_row['exon_coding_start']
    reading_frame = (-exon_coding_start)%3
    sites_data = sorted([s.split(',') for s in exon_row['sites_data'].split(':')], key = lambda x: int(x[1]))
#    sites_exon_positions = [int(s[1]) - exon_coding_start for s in sites_data]
    sites_exon_positions = [int(s[1]) - exon_coding_start for s in sites_data if s[4]+'_'+s[7]+'_'+str(s[6]) not in genomic_keys_for_attended_sites]
    codons_and_positions = [find_site_codon(exon_seq, reading_frame, s) for s in sites_exon_positions]
    
    genomic_keys_for_attended_sites += [s[4]+'_'+s[7]+'_'+str(s[6]) for s in sites_data if s[4]+'_'+s[7]+'_'+str(s[6]) not in genomic_keys_for_attended_sites]
    
    if len([k for k in sites_exon_positions if k < reading_frame]):
        print(exon_key + ' containsites before reading frame starts')
    
    codons = []
    codons_positions = []
    codons_not_appanded = []
    for c,p in codons_and_positions:
       if len(c[0])<3:  #if codon is not complete due to vecinity to exon end, complete using preceeding exon
           print(exon_key + ' ' + str(p) + ' ' + str(c) + ' is short')
           variants = exons[exons['ucsc_id'] == exon_row['ucsc_id']]
           proceeding_exon = variants[variants['exon_number_in_var']==exon_row['exon_number_in_var']+1].squeeze()
           if len(proceeding_exon):
               complete_codon = c[0] + proceeding_exon['exon'][0:3-len(c[0])]
               complete_codon_and_nuc_edit = (complete_codon,c[1])
               codons.append(complete_codon_and_nuc_edit)
               codons_positions.append(p)
               print('codon was completed using proceeding exon in var') 
           else:
               print('no proceeding exon in var to complete codon')
               codons_not_appanded.append((c,p))
       else:
           codons.append(c)
           codons_positions.append(p)
    codons_positions = sorted(list(set([x for x in codons_positions])))
#    codons = [x[0] for x in codons_and_positions]
#    codons_positions = sorted(list(set([x[1] for x in codons_and_positions]))) # using a sorted(list(set())) shrinks the list to single unique values of codons positions corresponding to merge_sites_in_same_codon output
    
    if len(codons):
        codons_merged_sites, sites_ex_position_merged_locations = merge_sites_in_same_codon(codons, reading_frame, sites_exon_positions)
        exons_edited_codons_dict.update({exon_key:{'edited_codons':codons_merged_sites,'corresponding_sites_locations':sites_ex_position_merged_locations,'corresponding_codons_locations':codons_positions}})
    if len(codons_not_appanded):
        print('data was not added to exons_edited_codons_dict: exon ' + exon_key + ' codons ' + str(codons_not_appanded))
        
    return genomic_keys_for_attended_sites

def search_codon_in_sequence(codon, sequence, reading_frame):
    """
    given a sequence, its reading_frame and a codon
    this function will return a list of all start positions of the occurances of the codon in the sequence
    """
    similar_codons_in_exons = []
    for i in range(reading_frame,len(sequence),3):
        if sequence[i:i+3] == codon:
            similar_codons_in_exons.append(i)
    return similar_codons_in_exons


def search_for_codons_in_same_exon(exon_row, exons_edited_codons_dict, chosen_new_codons, codos_not_found_in_exon, genomic_locations_of_codons_already_found,sites_loc_wrt_genomic_locations_of_sites_already_found):
    """
    This function takes an exon_row from exons dataframe, and three dictionaries
    1. exons_edited_codons_dict - output of find_all_edited_codons_in_exon function
    2. chosen_new_codons - a nested dict similar to exons_edited_codons_dict (without the corresponding_sites_locations locations which could actually be constructed using the other 2 lists)
    3. codos_not_found_in_exon - a "slice" of exons_edited_codons_dict conataining all codons in exons_edited_codons_dict for which a similar codon to edit in exon could not be found (here as well with no corresponding_sites_locations list)
    
    for each real edited codon in exons_edited_codons_dict, this function tries to find all possible codons in exon in exon_row that could be set as new edited codons (do not already contain sites and were not already set as new edited codon) 
    and append them randomaly (for each real edited codon) as new sites to chosen_new_codons
    if there are no similar codons left to choose from, it appends the codon to codos_not_found_in_exon
    """ 
    exon_key = exon_row['exon_key']
    exon_seq = exon_row['exon'].replace('a','A').replace('g','G').replace('c','C').replace('t','T')
    exon_coding_start = exon_row['exon_coding_start']
    genomic_exon_start = exon_row['genomic_exon_start']
    strand = exon_row['strand']   
    chromosome = exon_row['chr']
    genomic_exon_end = exon_row['genomic_exon_end']

    reading_frame = (-exon_coding_start)%3
    all_real_codons_positions = exons_edited_codons_dict[exon_key]['corresponding_codons_locations']
       
    for j, cod in enumerate(exons_edited_codons_dict[exon_key]['edited_codons']):
        
        codon_location_in_exon = exons_edited_codons_dict[exon_key]['corresponding_codons_locations'][j]
        
        #calc codon genomic location
        if strand == '+':
            genomic_position_base0 = genomic_exon_start + codon_location_in_exon
            genomic_position_base1 = genomic_position_base0 + 1
        elif strand == '-':
            genomic_position_base0 = genomic_exon_end - codon_location_in_exon - 1
            genomic_position_base1 = genomic_position_base0 + 1
        genomic_codon_key_base1 = chromosome+'_'+strand+'_'+str(genomic_position_base1)
        
        genomic_codon_already_attended = False
        if genomic_codon_key_base1 in genomic_locations_of_codons_already_found:
            genomic_codon_already_attended = True
            print('codon ' + genomic_codon_key_base1 + ' was already taken cate of')
        
        if not genomic_codon_already_attended:
        
            possible_codons_to_edit = search_codon_in_sequence(cod[0], exon_seq, reading_frame)
            
            #get codons position of codons that were already chosen in prevevious iterations 
            if exon_key in chosen_new_codons.keys():
                chosen_new_codons_for_exon = chosen_new_codons[exon_key]['corresponding_codons_locations']
            else:
                chosen_new_codons_for_exon = []
            
            #all codons first nuc position in exon (if codon do not contain any real site or already drawn site)
            codons_to_draw = [i for i in possible_codons_to_edit if i not in all_real_codons_positions + chosen_new_codons_for_exon] 
            
            if len(codons_to_draw):
                if exon_key in chosen_new_codons.keys(): 
                    chosen_new_codons[exon_key]['corresponding_codons_locations'].append(random.choice(codons_to_draw))    
                    chosen_new_codons[exon_key]['edited_codons'].append(cod)    
                else:
                    chosen_new_codons.update({exon_key:{'edited_codons':[cod],'corresponding_codons_locations':[random.choice(codons_to_draw)]}})
                genomic_locations_of_codons_already_found.append(genomic_codon_key_base1)
                for p in cod[1:]:
                    sites_loc_wrt_genomic_locations_of_sites_already_found.append(genomic_codon_key_base1+'_'+str(p))
                
            else:
                if exon_key in codos_not_found_in_exon.keys():
                    codos_not_found_in_exon[exon_key]['edited_codons'].append(cod)
                    codos_not_found_in_exon[exon_key]['corresponding_codons_locations'].append(codon_location_in_exon) 
                else:
                    codos_not_found_in_exon.update({exon_key:{'edited_codons':[cod],'corresponding_codons_locations':[codon_location_in_exon]}})
    
    return genomic_locations_of_codons_already_found,sites_loc_wrt_genomic_locations_of_sites_already_found
    

def search_for_codons_in_different_exons(cod, original_exon_row, codon_location_in_original_exon, other_exons, exons_edited_codons_dict, chosen_new_codons, codons_not_found, genomic_keys_codons_already_added,sites_loc_wrt_genomic_locations_of_sites_already_found):
    """
    similar to search_for_codons_in_same_exon
    but instrad searching new codons possibilities for each codon in exon,
    it takes a single codon, and datagrame of other exons which did not contain the codon and tries to find codon possibilities within thos exons (randomaly).
    taking into consideration all forbiden codons in exons_edited_codons_dict and chosen_new_codons
    and append failiours (to find similar codons in exon) to codos_not_found_in_entire_gene
    """
    def check_exons_overlap(original_exon_row, other_row):
        """
        checks if original_exon overlap the other exons in which we try to find corresponding sites
        in order to avoid choosing one of the original sites in the original exons.
        Assumes original_exon_row and other_row are two exons of the same gene (check before calling this function)
        """
        overlap = True
        other_start = other_row['genomic_exon_start']
        other_end = other_row['genomic_exon_end']
        original_start = original_exon_row['genomic_exon_start']
        original_end = original_exon_row['genomic_exon_end']
        
        if (other_start<original_start and other_end<original_start) or (other_start>original_end and other_end>original_end):
            overlap = False

        return overlap

    original_exon_key = original_exon_row['exon_key']
    genomic_exon_start = original_exon_row['genomic_exon_start']
    strand = original_exon_row['strand']   
    chromosome = original_exon_row['chr']
    genomic_exon_end = original_exon_row['genomic_exon_end']
        
    #calc codon genomic location
    if strand == '+':
        genomic_position_base0 = genomic_exon_start + codon_location_in_original_exon
        genomic_position_base1 = genomic_position_base0 + 1
    elif strand == '-':
        genomic_position_base0 = genomic_exon_end - codon_location_in_original_exon - 1
        genomic_position_base1 = genomic_position_base0 + 1
    genomic_codon_key_base1 = chromosome+'_'+strand+'_'+str(genomic_position_base1)
    
    genomic_codon_already_attended = False
    if genomic_codon_key_base1 in genomic_keys_codons_already_added:
        genomic_codon_already_attended = True
        print('codon ' + genomic_codon_key_base1 + ' was already taken cate of')
    
    other_exons.set_index('exon_key', inplace = True)    
    exons_list = list(other_exons.index)
    np.random.shuffle(exons_list)       
 
    new_edited_codon_added = False   
    
    for ex in exons_list:
#        print(ex)       
        exon_row = other_exons.loc[ex,:]
        ovelapping_exons = False
        if original_exon_row['HGNC_protein_id'] == exon_row['HGNC_protein_id']: #if exons are of the same gene, dont search corresponding codons in overlapping exons
            ovelapping_exons = check_exons_overlap(original_exon_row, exon_row)
            
        if not genomic_codon_already_attended and not ovelapping_exons:

            exon_key = ex
            exon_seq = exon_row['exon'].replace('a','A').replace('g','G').replace('c','C').replace('t','T')
            exon_coding_start = exon_row['exon_coding_start']
            reading_frame = (-exon_coding_start)%3
            
            if exon_key in exons_edited_codons_dict.keys():
                all_real_codons_positions = exons_edited_codons_dict[exon_key]['corresponding_codons_locations']
            else:
                all_real_codons_positions = []
            
            possible_codons_to_edit = search_codon_in_sequence(cod[0], exon_seq, reading_frame)
            
            if exon_key in chosen_new_codons.keys():
                chosen_new_codons_for_exon = chosen_new_codons[exon_key]['corresponding_codons_locations']
            else:
                chosen_new_codons_for_exon = []
            
            codons_to_draw = [i for i in possible_codons_to_edit if i not in all_real_codons_positions + chosen_new_codons_for_exon]
            
            if len(codons_to_draw):
                new_edited_codon_added = True
                if exon_key in chosen_new_codons.keys(): 
                    chosen_new_codons[exon_key]['corresponding_codons_locations'].append(random.choice(codons_to_draw))    
                    chosen_new_codons[exon_key]['edited_codons'].append(cod)    
                else:
                    chosen_new_codons.update({exon_key:{'edited_codons':[cod],'corresponding_codons_locations':[random.choice(codons_to_draw)]}})
                
                genomic_keys_codons_already_added.append(genomic_codon_key_base1)
                for p in cod[1:]:
                    sites_loc_wrt_genomic_locations_of_sites_already_found.append(genomic_codon_key_base1+'_'+str(p))
                
                break

#    print(new_edited_codon_added)        
    #appending codons to codons_not_found dictionary in order to try and find a corresponding codon in a proceeding stage (or not, if codon was already attended due to overlapping exons containing it)
    if not new_edited_codon_added:
        if original_exon_key in codons_not_found.keys():
            codons_not_found[original_exon_key]['edited_codons'].append(cod)
            codons_not_found[original_exon_key]['corresponding_codons_locations'].append(codon_location_in_original_exon) 
        else:
            codons_not_found.update({original_exon_key:{'edited_codons':[cod],'corresponding_codons_locations':[codon_location_in_original_exon]}})
    
    return genomic_keys_codons_already_added,sites_loc_wrt_genomic_locations_of_sites_already_found




def retrive_site_row_for_proteomics_simulator_format(exon_row_for_calc, codon, codon_location_in_exon_for_calc, exon_occurance, mm_type = 'AG', qa_genomic_sites_list = [], genomic_keys_already_added = [], qa = False):
    """
    This function arguments are
    exon_row_for_calc = the row from the exons dataframe based on wihch the new site was created
    codon = a tuple, first element is the codon sequence, rest of it are locations in codon that rae edited
    codon_location_in_exon - base0 location in exon
    exon_occurance - a row from exons dataframe which represent some occurance of the codon
    mm_type - the editings type in the codon.
    
    it returns a list of tuples - each tuple is all the data needed to crate a row in the sites table format for the proteomics simulator
    (max returned list size is 3)
    """

    #values shared across all variants containing exon (mm_type is shared as well)  
    chromosome = exon_row_for_calc['chr']
    strand = exon_row_for_calc['strand']    
    editing_level = 0.
    HGNC_protein_id = exon_row_for_calc['HGNC_protein_id']
    
    codon_seq = codon[0]
    editings_in_codon = [i for i in codon[1:]]
    
    genomic_exon_start = exon_occurance['genomic_exon_start']
    genomic_exon_end = exon_occurance['genomic_exon_end']
    exons_gap_start = genomic_exon_start - exon_row_for_calc['genomic_exon_start']
    exons_gap_end = genomic_exon_end - exon_row_for_calc['genomic_exon_end']
    
    if strand == '+':
        codon_location_in_exon = codon_location_in_exon_for_calc - exons_gap_start
    elif strand == '-':
        codon_location_in_exon = codon_location_in_exon_for_calc + exons_gap_end
    
    #for qa - print if different codon thanexpected
    codon_seq_in_exon_occurance = exon_occurance['exon'][codon_location_in_exon:codon_location_in_exon+3].replace('a','A').replace('g','G').replace('c','C').replace('t','T')
    if len(codon_seq_in_exon_occurance)<3:
        if codon_seq_in_exon_occurance != exon_row_for_calc['exon'][codon_location_in_exon_for_calc:codon_location_in_exon_for_calc+3]:
            print(exon_occurance['exon_key'] +  ' codon:'+codon_seq_in_exon_occurance+str(codon_location_in_exon_for_calc) + ' has different codon than expected from ' + exon_row_for_calc['exon_key'] + ' codon:'+codon_seq+str(codon_location_in_exon_for_calc))
            print(str(exons_gap_end)+ ' ' + str(exons_gap_start))
    elif codon_seq_in_exon_occurance != codon_seq:
        print(exon_occurance['exon_key'] +  ' codon:'+codon_seq_in_exon_occurance+str(codon_location_in_exon_for_calc) + ' has different codon than expected from ' + exon_row_for_calc['exon_key'] + ' codon:'+codon_seq+str(codon_location_in_exon_for_calc))
        print(str(exons_gap_end)+ ' ' + str(exons_gap_start))
        
    sites_data = []
    if len(codon_seq)<3:
        print(exon_row_for_calc.name + ' ' + str(codon_location_in_exon) + ' is short')
    elif codon_seq in ['TAG','TAA','TGA']:
        print(exon_row_for_calc.name + ' ' + str(codon_location_in_exon) + ' is stop codon')
    else:
        #for each edit in exon, create a tuple of all data needed for sites table
        for edit in editings_in_codon:
            variant_name = exon_occurance['variant_name']
            ucsc_id = exon_occurance['ucsc_id']
            
            codon_as_list = list(codon_seq)
            codon_as_list[edit] = mm_type[1]
            edited_codon = ''.join(codon_as_list)
            aa_change = str(Seq(codon_seq,generic_dna).translate()) + str(Seq(edited_codon,generic_dna).translate()) 
            
            position_base0 = exon_occurance['exon_coding_start'] + codon_location_in_exon + edit
            position_base1 = position_base0 + 1
            coding_key_base1 = variant_name+'_'+str(position_base1)
    
            if strand == '+':
                genomic_position_base0 = genomic_exon_start + codon_location_in_exon + edit
                genomic_position_base1 = genomic_position_base0 + 1
            elif strand == '-':
                genomic_position_base0 = genomic_exon_end - codon_location_in_exon - edit - 1
                genomic_position_base1 = genomic_position_base0 + 1
            
            genomic_key_base1 = chromosome+'_'+strand+'_'+str(genomic_position_base1)    
            
            sites_data.append((variant_name,position_base0,position_base1,mm_type,chromosome,genomic_position_base0,genomic_position_base1,strand,aa_change,editing_level,genomic_key_base1,coding_key_base1,ucsc_id,HGNC_protein_id))
                
            #this is for qa. if this genomic site is suspecious and passed to function in qa_genomic_sites_list, print some data on it
            if qa:
                if genomic_key_base1 in qa_genomic_sites_list:
                    print('genomic key was specified wass passed to qa_genomic_sites_list')
                    print(exon_row_for_calc.name)
                    print(genomic_key_base1 + '  ' + coding_key_base1)
                    print(str(codon) + '  location inside exon:' + str(codon_location_in_exon))
                    
                if genomic_key_base1 in genomic_keys_already_added:
                    print('genomic key already added')
                    print(exon_row_for_calc.name)
                    print(genomic_key_base1 + '  ' + coding_key_base1)
                    print(str(codon) + '  location inside exon:' + str(codon_location_in_exon))
                
                if len(editings_in_codon)>1:
                    print('codon_contain_multiple editing site')
                    print(sites_data)
                    
                if len(sites_data) == 0:
                    print('sites data was not found for exon')
                    print(exon_row_for_calc['exon_key'])
                
    return sites_data


def exons_codons_dict_to_sites_df(exons_edited_codons_dict,all_exons_no_dup, all_exons, qa_genomic_sites_list = []):
    
    """
    This function takes
    exons_edited_codons_dict - exons to codons nested dictionaty
    all_exons_no_dup - a dataframe of exons based on which ethe xons_edited_sites dictionary was created
    exons - a dataframe of all exons occurancers (including duplicates which were filtered out to create all_exons_no_dup)
    qa_genomic_sites_list - a list of genomic sites keys for which to print some of the data calculated for the sites_df format - this is one of the retrive_site_row_for_proteomics_simulator_format parameters
    
    This function then creates a sites dataframe in simulator format
    sites in some exon will be associated to each of the variants containing the exon using all variants from exons dataframe
    """
    
    def find_all_exons_containing_codon(exon_row_for_calc, codon_location_in_exon_for_calc, all_exons_of_gene):
        strand = exon_row_for_calc['strand']
        genomic_exon_start = exon_row_for_calc['genomic_exon_start']
        genomic_exon_end = exon_row_for_calc['genomic_exon_end']
        if strand == '+':
            genomic_position_base0 = genomic_exon_start + codon_location_in_exon_for_calc
        elif strand == '-':
            genomic_position_base0 = genomic_exon_end - codon_location_in_exon_for_calc - 1
        all_exons_of_gene = all_exons_of_gene[all_exons_of_gene['genomic_exon_start'] <= genomic_position_base0]
        all_exons_of_gene = all_exons_of_gene[all_exons_of_gene['genomic_exon_end'] >= genomic_position_base0]
        
        return all_exons_of_gene
    
    
    genomic_keys_added = []
    sites_data = []
    for k, v in exons_edited_codons_dict.items():
        exon_row_for_calc = all_exons_no_dup.loc[k,:]
        all_exons_of_gene = all_exons[all_exons['HGNC_protein_id'] == exon_row_for_calc['HGNC_protein_id']]
#        all_exon_occurances = all_exons[all_exons['exon_key']==k]
        for j, codon in enumerate(v['edited_codons']):         
            all_exon_occurances = find_all_exons_containing_codon(exon_row_for_calc, v['corresponding_codons_locations'][j], all_exons_of_gene)   
            for index, exon_occurance in all_exon_occurances.iterrows():
                new_data = retrive_site_row_for_proteomics_simulator_format(exon_row_for_calc, codon, v['corresponding_codons_locations'][j],exon_occurance, qa_genomic_sites_list = qa_genomic_sites_list, genomic_keys_already_added = genomic_keys_added)
                sites_data += new_data
    new_sites_tbl_columns = ['variant_name','position_base0','position_base1','mm_type','chromosome','genomic_position_base0','genomic_position_base1','strand','aa_change','edigin_level','genomic_key_base1','coding_key_base1','ucsc_id','HGNC_protein_id']    
    sites_df = pd.DataFrame(data = sites_data, columns = new_sites_tbl_columns)
    print('Unique genomic sites created: ' + str(len(set(sites_df['genomic_key_base1']))))
    
    return sites_df
    

def check_new_codons_vs_original_ones(exons_edited_codons_dict, chosen_new_codons, all_exons_no_dup, exons):
    """
    This function recives two dictionaries of exons and edited codons, constructing a sites table from each
    and check overlap of those tables 
    
    """
    #reconstructing original sites dataframe from exons_edited_codons_dict dictionary
    original_sites = exons_codons_dict_to_sites_df(exons_edited_codons_dict,all_exons_no_dup,exons)
    original_genomic_sites = set(original_sites['genomic_key_base1'])
    #constructing new sites dataframe from exons_edited_codons_dict dictionary
    new_sites = exons_codons_dict_to_sites_df(chosen_new_codons,all_exons_no_dup,exons)
    new_genomic_sites = set(new_sites['genomic_key_base1'])    
    
    new_in_original = [i for i in new_genomic_sites if i in original_genomic_sites]  
    bad_new_sites = new_sites[new_sites['genomic_key_base1'].isin(new_in_original)]
    corresponding_original_sites = original_sites[original_sites['genomic_key_base1'].isin(new_in_original)]

    #creating sites table again. now printing data for problematic codons
    print('\nBad codons data from Original sites:')
    original_sites = exons_codons_dict_to_sites_df(exons_edited_codons_dict,all_exons_no_dup,exons, qa_genomic_sites_list=new_in_original)
    
    print('\nBad codons data from New sites:')
    new_sites = exons_codons_dict_to_sites_df(chosen_new_codons,all_exons_no_dup,exons, qa_genomic_sites_list=new_in_original)
    

if __name__ == '__main__':
    
#    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Preparation of random AG sites list for proteomics simulator input')
#    run_parser = parser.add_argument_group('generate random AG table')
#    run_parser.add_argument('-exons', dest='exons_file', action='store', required = True, help='path to exons file')
#    run_parser.add_argument('-sites', dest='sites_file', action='store', required = True, help='path to sites file')
#    run_parser.add_argument('-o', dest='output_name', action='store', default = 'random_ag_editings', help='output table name')
#    arguments = parser.parse_args()
#        
#    exons_file = arguments.exons_file
#    sites_file = arguments.sites_file
#    output_name = arguments.output_name
    
    sites_file = 'C:/Users/shosh/OneDrive/Desktop/test/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor_with_additional_data.txt'
    exons_file = 'C:/Users/shosh/OneDrive/Desktop/test/exons_and_ag_sites.txt'
    output_name = 'random_ag_editings'
    out_path = '/'.join(exons_file.split('/')[0:-1]) + '/'
    
    print("Reading sites data")
    sites = read_non_editing_sites_wrt_coding_seqs(sites_file)
    
    print("Reading all transcriptome's exons data")
    exons_columns = ['variant_name','exon_coding_start','exon_coding_end','exon_number_in_var','strand','chr','genomic_exon_start','genomic_exon_end','exon','AG_sites_num','sites_data']
    exons = pd.read_csv(exons_file, sep = '\t', names = exons_columns)
    exons = exons.apply(lambda row: split_gene_and_prot_names(row), axis = 1)
    exons['exon_key'] = exons.apply(lambda row: row['HGNC_protein_id']+'_'+row['chr']+'_'+row['strand']+'_'+str(row['genomic_exon_start'])+'_'+str(row['genomic_exon_end']), axis = 1)
    
    edited_exons = exons[exons['sites_data'].notnull()]  #first taking all edited exons - some exons contain edited site while appearing in some variant but do not contain any site in while in other variant - should check that it is not a bug in orshays commands for 
    edited_exons_no_dup = edited_exons.groupby('exon_key', group_keys=False).apply(lambda row: row.loc[row['AG_sites_num'].idxmax()])  #removing duplicates of exons keeping the one with largest number of sites
    unedited_exons_no_dup = exons[~exons['exon_key'].isin(edited_exons_no_dup['exon_key'])].drop_duplicates('exon_key')
    unedited_exons_no_dup.set_index('exon_key', inplace = True, drop = False)
    all_exons_no_dup = pd.concat([edited_exons_no_dup,unedited_exons_no_dup])
    
    codons_attended = []
    sites_attended = []
    #create exons_edited_codons_dict nested dictionary
    all_genomic_keys = []
    exons_edited_codons_dict = {}
    for index, exon_row in edited_exons_no_dup.iterrows():
        all_genomic_keys = find_all_edited_codons_in_exon(exon_row,exons_edited_codons_dict,all_genomic_keys)
    
    #senity check - number of unique sites from sites dataframe should be identical to uique sites from exons dataframe
    n_genomic_sites_from_exons_tbl = 0
    for k, v in exons_edited_codons_dict.items():  
        all_sites_in_exons = flatten_sites_positions_in_exon(k,v['corresponding_sites_locations'])
        n_genomic_sites_from_exons_tbl += len(all_sites_in_exons)
    print('Total sites from sites list: ' + str(len(set(sites[sites['mm_type'] == 'AG']['genomic_key_base1']))))
    print('Total sites from exons list: ' + str(n_genomic_sites_from_exons_tbl))    
    original_sites_df_from_exons = exons_codons_dict_to_sites_df(exons_edited_codons_dict,all_exons_no_dup,exons) 
    original_genomic_sites = set(original_sites_df_from_exons['genomic_key_base1'])
    print('Expecting ' + str(len(original_genomic_sites)) + ' New genomic sites')
    

    #create a dictionary of a similar structure of all chosen codons codons - new codons containing new editing sites
    #and a ductionaty of codons for which a similar codon in the exon could not be found    
    codons_not_found_in_exon = {}
    chosen_new_codons = {}
    for index, exon_row in edited_exons_no_dup.iterrows():
        if exon_row['exon_key'] in exons_edited_codons_dict.keys():
            codons_attended,sites_attended = search_for_codons_in_same_exon(exon_row, exons_edited_codons_dict, chosen_new_codons, codons_not_found_in_exon, codons_attended,sites_attended)
    new_codons_in_same_exons_with_overlapping_exons_count_twice = sum([len(v['corresponding_codons_locations']) for v in chosen_new_codons.values()])
    new_codons_in_same_exons_no_dup = len(set(codons_attended))
    new_sites_in_same_exons_no_dup = len(set(sites_attended))
    print('New codons in same exons: ' + str(new_codons_in_same_exons_no_dup) + ' ('+str(new_sites_in_same_exons_no_dup)+' sites)')
    
    #for codons_not_found_in_exon, try to find similar codons in other exons of the same gene
    codons_not_found_in_entire_gene = {}
    for k,v in codons_not_found_in_exon.items():
#        print('\noringinal exon ' + k)
        gene_name = k.split('_')[0]
        other_exons_in_gene = all_exons_no_dup[np.logical_and(all_exons_no_dup['HGNC_protein_id']==gene_name, all_exons_no_dup['exon_key'] != k)]
        original_exon_row = all_exons_no_dup.loc[k,:]
        if len(other_exons_in_gene):
            for i, codon in enumerate(v['edited_codons']):
                codons_attended,sites_attended = search_for_codons_in_different_exons(codon, original_exon_row, v['corresponding_codons_locations'][i], other_exons_in_gene.copy(), exons_edited_codons_dict, chosen_new_codons, codons_not_found_in_entire_gene, codons_attended,sites_attended)
        else:
            codons_not_found_in_entire_gene.update({k:v})        
    new_codons_in_other_exons_in_same_gene_no_dup = len([a for a in set(codons_attended) if a is not None]) - new_codons_in_same_exons_no_dup
    new_sites_in_other_wxons_in_same_gene_no_dup = len([a for a in set(sites_attended) if a is not None]) - new_sites_in_same_exons_no_dup
    print('New codons in other exons of the same gene: ' + str(new_codons_in_other_exons_in_same_gene_no_dup) + ' ('+str(new_sites_in_other_wxons_in_same_gene_no_dup)+' sites)')
     
    
    #for codons_not_found_in_entire_gene, try to find similar codons in other genes
    codons_not_found_in_entire_data_set = {}
    for k,v in codons_not_found_in_entire_gene.items():
#        print('\noringinal exon ' + k)
        gene_name = k.split('_')[0]
        other_genes_exons = all_exons_no_dup[all_exons_no_dup['HGNC_protein_id'] != gene_name]
        original_exon_row = all_exons_no_dup.loc[k,:]
        for i, codon in enumerate(v['edited_codons']):
            codons_attended,sites_attended = search_for_codons_in_different_exons(codon, original_exon_row, v['corresponding_codons_locations'][i], other_genes_exons.copy(), exons_edited_codons_dict, chosen_new_codons, codons_not_found_in_entire_data_set, codons_attended,sites_attended)
    new_codons_in_other_genes = len([a for a in set(codons_attended) if a is not None]) - new_codons_in_same_exons_no_dup - new_codons_in_other_exons_in_same_gene_no_dup
    new_sites_in_other_genes = len([a for a in set(sites_attended) if a is not None]) - new_sites_in_same_exons_no_dup - new_sites_in_other_wxons_in_same_gene_no_dup
    print('New codons in other genes: ' + str(new_codons_in_other_genes) + ' ('+str(new_sites_in_other_genes)+' sites)')
    
    original_edited_codons = sum([len(v['corresponding_codons_locations']) for v in exons_edited_codons_dict.values()])
    new_edited_codons = sum([len(v['corresponding_codons_locations']) for v in chosen_new_codons.values()])
    print('Original edited codons (containing at least one editing site): ' + str(original_edited_codons))
    print('New edited codons (containing at least one editing site): ' + str(new_edited_codons))

    #from chosen_new_codons - prepare a table similar to real sites table
    print('Preparing perliminary editing sites table format for proteomics simulator\n(containing only instances of sites for specific exons for which new codons were found, not including instances in overlapping exons)')
#    all_exons_no_dup.set_index('exon_key', inplace =  True)
    new_sites = exons_codons_dict_to_sites_df(chosen_new_codons,all_exons_no_dup,exons)
    
    #writing only columns expected by the simulator input preparation script
    columns_to_write = ['variant_name','position_base0','position_base1','mm_type','chromosome','genomic_position_base0','genomic_position_base1','strand','aa_change','edigin_level']
    new_sites[columns_to_write].to_csv(out_path + output_name+'.txt', sep = '\t', index = False, header = False)
    print('New editing sites file was written: ' + out_path + output_name+'.txt')

    

    