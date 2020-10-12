import os
import operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Charts import plot_expansion_polarization, plot_reads_per_seq
from CleaveFastaPeptidesFunctions import get_chain, get_clone, get_reads
from Bio.Alphabet import generic_protein
from Proteomics_Simulation_Params import *    



"""
create_clone_exp_pol_dict
from source file:
get a dictionary of all clones of a certain chain
for each - value is a dictionary of all sequences id's of the clone and for each sequence in doctionary, value is a list of all subset reads
"""
#d = create_clone_sequences_reads_dict(files_dir + 'source_sequences/','all_subset_unique_sequences_non_single.fasta','IGH',clone_identifier='CDR3: ')
def create_clone_sequences_reads_dict(path,file_name,chain,clone_identifier = 'clone: '):
    
#    print('\n')
#    print('collecting sequences for all clones of selected chain...')
    
    clone_sequences_reads = {}
        
    for record in SeqIO.parse(open(path + file_name, "r"), "fasta"):

        #assuming all sequences of certain cdr3 are of the same chain
        if get_chain(str(record.description)) == chain:
            
            clone = get_clone(str(record.description), clone_identifier)
            reads = get_reads(str(record.description))
            
            #update cdr and its sequences to clone in clone_sequences_reads
            if clone in clone_sequences_reads:
                clone_sequences_reads[clone].update({str(record.id):reads})
            else:
                clone_sequences_reads.update({clone:{str(record.id):reads}})
                
    return clone_sequences_reads


        
                
"""
data_size
for files containing unique peptides cleaved by generate_cleaved_peptide_files
return
n - number of unique peptides created
k - number of peptides created from sources before srhinking data base to unique peptides
nember of unique sources from which data base was created
"""
#x,y,z = data_size('C:/linux_share/Proteomics_project/files/cleaved_peptides/','trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta')
def data_size(path, file_name):
    
    n=0 #cnt number of unique peptides
    k=0 #cnt total number of peptides (NOT UNIQUE PEPTIDES!!)
    sources_list = [] #list countaining all sources id's
    
    for record in SeqIO.parse(open(path + file_name, "r"), "fasta"):
            
        n += 1
        
        if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
            des = record.description[len(record.id):]
            data = eval(des)
        else:
            data = eval(record.description)
            
        for source in data:
            sources_list.append(source['source_id'])
            k += len(source['overlaps']) #number of appearances of peptide for one source is number of overlaps with that source

    return n, k, len(set(sources_list))


             

"""
MolWeightList
input - path
        file - file with peptides (after creating for each if clone_informative or coverage informative), description format:
                {'informative': clone_informative, 'clone_coverage': 0, 'data':[sources_list]}
        inf - one of the groups: clove_informative\clone_informative\0
        
create a list of molecular weights for each peptide in group
""" 

#mw_list = MolWeightList(files_dir+"cleaved_peptides/",'inf_lib_heavy_chains_from_trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta')
def molecular_weight_list(path,file,groups = ['clone_informative','coverage_informative']):
    
    mw_list = []
    amb_records_list = []
    n=0
    k = 0
    
    print('\n') 
    print('analyzing molecular weight in file: ' + file)
    
    for record in SeqIO.parse(open(path + file, "r"), "fasta"):
            
        n += 1
            
        if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
            des = record.description[len(record.id):]
            rec_des = eval(des)
        else:
            rec_des = eval(record.description)
    
        if rec_des['informative'] in groups:
            k += 1
            if any(substring in str(record.seq) for substring in ['X','B','Z','J','_','*']):
                amb_records_list.append(str(record.id))
            else:
                mw_list.append(round(ProteinAnalysis(str(record.seq)).molecular_weight(),3))
    
               
    print('total records in file: ' + str(n))
    print('total informative peptides in file: ' + str(k))            
    print('informative peptides with defined molecular weight: ' + str(len(mw_list)))            
    print('informative peptides with ambiguous molecular weight: ' + str(len(amb_records_list)))
    
    return mw_list


