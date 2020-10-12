import os
import sys
import numpy as np
import operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ScanFile import create_clone_sequences_reads_dict
from Bio.Alphabet import generic_protein
from CleaveFastaPeptidesFunctions import get_clone, get_reads
from Charts import plot_expansion_polarization, plot_reads_per_seq, print_dict_to_xlsx, plot_expansion_polarization_rps_residuals
from Proteomics_Simulation_Params import *


"""
get_clones
use create_clone_sequences_reads_dict to get a clone-seqences dictionary
creates a dictionary of all clones from specified subset (have sequences that have reads in those subsets)
value is a dictionary of the keys:
    expansion - percentage of clone sequences from total sequences for specified chain and subset
    polarization - percentage of clone reads (sum of its sequences reads) from total reads of sequences of specified chain and subsets
"""
#dtest = get_clones(father_path,'final_fasta_mouse.fasta','IGH',clone_identifier = 'CDR3: ')
#dtest = get_clones(father_path,'final_fasta_mouse_simulation.fasta','IGH',346596,clone_identifier = 'CDR3: ')
#dtest = get_clones(father_path,'all_subset_unique_sequences_non_single.fasta','IGH',231,clone_identifier = 'CDR3: ')
def get_clones(path,file_name,chain,*subsets,clone_identifier = 'clone: ',plot = True):
    
    #define input folder
#    if not os.path.isdir(path + "source_sequences"):
#        print('source_sequences folder was not found') #print error if sources folder does not exist
#    input_path = path + "source_sequences/"
    input_path = path
    
    clone_sequences_reads = create_clone_sequences_reads_dict(input_path,file_name,chain,clone_identifier = clone_identifier)
    
    #create output path
    if not os.path.isdir(path + 'clones_analysis_for_' +  file_name[:-6]):
        os.makedirs(path + 'clones_analysis_for_' +  file_name[:-6])
    output_path = path + 'clones_analysis_for_' +  file_name[:-6] + '/'
    
    print('\n')
    print('collecting clones of selected subsets...')
    
    temp_clone_exp_pol_dict = {}
    clone_exp_pol_dict = {}    
    total_reads = 0
    total_seqs = 0
    
    
    #create a string contaning chosen subsets to be used in some files writing
    if len(subsets):
        subsets_str = 'subset(s) ' + str(subsets)
        subsets_str = subsets_str[:-1] + ')'
        print('\n')
        print('calculating expansion\polarization values for clones in ' + subsets_str)
    else:
        print('\n')
        print('no subset specified, calculating expansion\polarization values for clones in all subsets')
        subsets_str = 'all subsets'
        
    ex_pol_dist_name = 'Polarization-Expansion of ' + chain + ' clones from ' + subsets_str
    ex_pol_dist_res_name = 'Polarization-Expansion RPS Residuals of ' + chain + ' clones from ' + subsets_str
    rps_dist_name = 'Reads per Sequence of ' + chain + ' clones from ' + subsets_str
    
    #iterate over all clones in clone_sequences_reads   
    for clone in clone_sequences_reads:
        sequences_for_clone = 0
        reads_for_clone = 0
#        print(clone)
        
        #iterate over all sequences of certain clone
        for sequence in clone_sequences_reads[clone]:
            reads = 0
#                print(clone_sequences_reads[clone][query][sequence])
            if len(subsets):   
                reads = sum(subset_reads[1] for subset_reads in clone_sequences_reads[clone][sequence] if subset_reads[0] in subsets)
            else:
                reads = sum(subset_reads[1] for subset_reads in clone_sequences_reads[clone][sequence])
                          
            if reads:
                #sum all sequences for the clone that have reads in the specified subsets
                sequences_for_clone += 1
                reads_for_clone += reads
        
#        print(str(sequences_for_clone))
#        print(str(reads_for_clone))
        total_seqs += sequences_for_clone
        total_reads += reads_for_clone
            
        if sequences_for_clone:
            temp_clone_exp_pol_dict.update({clone:{'seqs':sequences_for_clone,'reads':reads_for_clone}})
          
    #create final clone_exp_pol_dict with expansion and polariztion values for each clone
    for clone in temp_clone_exp_pol_dict:
        clone_exp_pol_dict.update({clone:{'seqs':temp_clone_exp_pol_dict[clone]['seqs'],'reads':temp_clone_exp_pol_dict[clone]['reads'],'expansion':100*temp_clone_exp_pol_dict[clone]['seqs']/total_seqs,'polarization':100*temp_clone_exp_pol_dict[clone]['reads']/total_reads,'reads_per_seq':temp_clone_exp_pol_dict[clone]['reads']/temp_clone_exp_pol_dict[clone]['seqs']}})

    #run next lines only if colones from specified chain and subset were found 
    if len(clone_exp_pol_dict):
        
        population_reads_per_seq = total_reads/total_seqs
            
        if plot:
            
            print('total reads for subsets: ' + str(total_reads))
            print('total sequences for subsets: ' + str(total_seqs))
            print('total clones per subset (' + chain + ' chain): ' + str(len(clone_exp_pol_dict)))
        
            #create expantion polarization scatter graph    
            plot_expansion_polarization(output_path,clone_exp_pol_dict,ex_pol_dist_name,ex_pol_dist_name.replace('\n','').replace(' ','_').replace('-','_').replace(',',''))
            plot_expansion_polarization_rps_residuals(output_path,clone_exp_pol_dict,ex_pol_dist_res_name,ex_pol_dist_res_name.replace('\n','').replace(' ','_').replace('-','_').replace(',','').replace('__','_'))
        
            #create reads per sequence graph
            reads_per_seq_list = sorted([clone_exp_pol_dict[clone]['reads_per_seq'] for clone in list(clone_exp_pol_dict.keys())],reverse=True)
            plot_reads_per_seq(output_path,reads_per_seq_list,population_reads_per_seq,rps_dist_name,rps_dist_name.replace('\n','').replace(' ','_').replace('-','_').replace(',','').replace('__','_'))
    
            #creat a sorted list (by expantion rate) of (clone,expantion,polarization,reads_per_seq)
            sorted_clones_by_expantion = sorted([(clone,clone_exp_pol_dict[clone]['expansion'],clone_exp_pol_dict[clone]['polarization'],clone_exp_pol_dict[clone]['reads_per_seq']) for clone in clone_exp_pol_dict], key=operator.itemgetter(1), reverse = True)
    
            #print log file of clones expantion-polarization
            with open(output_path + ex_pol_dist_name.replace('\n','').replace('-','_').replace(' ','_').replace(',','') +'.txt' , "w") as handle:
                handle.write('chain: ' + chain + '\n')
                if len(subsets):    
                    handle.write('subsets: ' + subsets_str[10:] + '\n')
                else:
                    handle.write('subsets: ' + subsets_str + '\n')
                handle.write('total number of sequences: ' + str(total_seqs) + '\n')
                handle.write('total number of reads: ' + str(total_reads) + '\n')
                handle.write('\n')
                handle.write('clones expantion-polarization data:' + '\n')
                for clone in sorted_clones_by_expantion:    
                    handle.write('clone: ' + clone[0] + '  expansion: ' + str(round(clone[1],3)) + '%  polarization: ' + str(round(clone[2],3)) + '% reads_per_seq: ' + str(round(clone[3],3)) + '\n')
            handle.close()
            
            #create an excel file
            clones_list = list(clone_exp_pol_dict.keys())
            xls_seqs = [clone_exp_pol_dict[clone]['seqs'] for clone in clone_exp_pol_dict]
            xls_reads = [clone_exp_pol_dict[clone]['reads'] for clone in clone_exp_pol_dict]
            xls_expansion = [clone_exp_pol_dict[clone]['expansion'] for clone in clone_exp_pol_dict]
            xls_polarization = [clone_exp_pol_dict[clone]['polarization'] for clone in clone_exp_pol_dict]
            xls_reads_per_seq = [clone_exp_pol_dict[clone]['reads_per_seq'] for clone in clone_exp_pol_dict]
        
            dictionary = {'Clones':clones_list,'Sequences':xls_seqs,'Reads':xls_reads,'Expansion':xls_expansion,'Polarization':xls_polarization,'reads_per_seq':xls_reads_per_seq}
            print_dict_to_xlsx(output_path,ex_pol_dist_name.replace('\n','').replace('-','_').replace(' ','_'),dictionary)
            
            print('log file for all clones was created')

    
    else:
        print('no clones were found for chain-' + chain + ' subset- ' + subsets_str + ' in file: ' + file_name)
    
    return clone_exp_pol_dict


"""
get_clones_by_params
same as get_clones but final results include only clones with expansion\polarization\reads_per_seq values specified in optional arguments
"""
#dtestparams = get_clones_by_params(father_path,'final_fasta_mouse_simulation.fasta','IGH',346596, expansion=[4,8], polarization=[4,8], reads_per_seq=[], clone_identifier = 'CDR3: ')
def get_clones_by_params(path,file_name,chain,*subsets, expansion=[0,100], polarization=[0,100], reads_per_seq=[], clone_identifier = 'clone: ',plot = True):
    
    chosen_clones_dict = {}
    total_seqs = 0
    total_reads = 0
            
    #create output path
    if not os.path.isdir(path + 'clones_analysis_for_' +  file_name[:-6]):
        os.makedirs(path + 'clones_analysis_for_' +  file_name[:-6])
    output_path = path + 'clones_analysis_for_' +  file_name[:-6] + '/'
    
    #create dict of clones that meet conditions:
    exp_pol_dict =  get_clones(path,file_name,chain,*subsets,clone_identifier = clone_identifier,plot = False)
    total_reads_without_cond = int(np.ceil(max([exp_pol_dict[clone]['reads']/exp_pol_dict[clone]['seqs'] for clone in exp_pol_dict])))
    
    #create a string contaning chosen subsets to be used in some files writing
    if len(subsets):
        subsets_str = 'subset(s) ' + str(subsets)
        subsets_str = subsets_str[:-1] + ')'
    else:
        subsets_str = 'all subsets' 
            
    if len(reads_per_seq):
        ex_pol_dist_name = 'Polarization-Expansion of ' + chain + ' clones from ' + subsets_str + ' \n' +  'exp ' + str(expansion[0]) + '-' + str(expansion[1]) + ' pol ' + str(polarization[0]) + '-' + str(polarization[1]) + ' RPS ' + str(reads_per_seq[0]) + '-' + str(reads_per_seq[1])
        ex_pol_dist_res_name = 'Polarization-Expansion RPS Residuals of ' + chain + ' clones from ' + subsets_str + ' \n' +  'exp ' + str(expansion[0]) + '-' + str(expansion[1]) + ' pol ' + str(polarization[0]) + '-' + str(polarization[1]) + ' RPS ' + str(reads_per_seq[0]) + '-' + str(reads_per_seq[1])
        ex_pol_dist_file_name = '_'.join(ex_pol_dist_name.splitlines()).replace(' ','_').replace('-','_').replace('__','_')
        rps_dist_name = 'Reads per Sequence of ' + chain + ' clones from ' + subsets_str + ' \n' + 'exp ' + str(expansion[0]) + '-' + str(expansion[1]) + ' pol ' + str(polarization[0]) + '-' + str(polarization[1]) + ' RPS ' + str(reads_per_seq[0]) + '-' + str(reads_per_seq[1])
        rps_dist_file_name = '_'.join(rps_dist_name.splitlines()).replace(' ','_').replace('-','_').replace('__','_')
    else:
        ex_pol_dist_name = 'Polarization-Expansion of ' + chain + ' clones from ' + subsets_str + ' \n' +  'exp ' + str(expansion[0]) + '-' + str(expansion[1]) + ' pol ' + str(polarization[0]) + '-' + str(polarization[1])
        ex_pol_dist_res_name = 'Polarization-Expansion RPS Residuals of ' + chain + ' clones from ' + subsets_str + ' \n' +  'exp ' + str(expansion[0]) + '-' + str(expansion[1]) + ' pol ' + str(polarization[0]) + '-' + str(polarization[1])
        ex_pol_dist_file_name = '_'.join(ex_pol_dist_name.splitlines()).replace(' ','_').replace('-','_').replace('__','_')
        rps_dist_name = 'Reads per Sequence of ' + chain + ' clones from ' + subsets_str + ' \n' + 'exp ' + str(expansion[0]) + '-' + str(expansion[1]) + ' pol ' + str(polarization[0]) + '-' + str(polarization[1])
        rps_dist_file_name = '_'.join(rps_dist_name.splitlines()).replace(' ','_').replace('-','_').replace('__','_')
    
    if not len(reads_per_seq):
        reads_per_seq = [1,total_reads_without_cond]
    
    #creating_dict_of_desirable_parameters_values_only
    for clone in exp_pol_dict:
        if all([expansion[0] <= exp_pol_dict[clone]['expansion'] <= expansion[1],
                polarization[0] <= exp_pol_dict[clone]['polarization'] <= polarization[1],
                reads_per_seq[0] <= exp_pol_dict[clone]['reads_per_seq'] <= reads_per_seq[1]]):
            
            total_seqs += exp_pol_dict[clone]['seqs']
            total_reads += exp_pol_dict[clone]['reads']
            chosen_clones_dict.update({clone:exp_pol_dict[clone]})
    
    if total_seqs:
        population_reads_per_seq = total_reads/total_seqs
    
    if len(chosen_clones_dict):
        if plot:  
            print(str(len(chosen_clones_dict)) + ' clones of specified parameters')
            print(str(total_seqs) + ' sequences for clones')
            
            #create expantion polarization scatter graph    
            plot_expansion_polarization(output_path,chosen_clones_dict,ex_pol_dist_name,ex_pol_dist_file_name)
            plot_expansion_polarization_rps_residuals(output_path,chosen_clones_dict,ex_pol_dist_res_name,ex_pol_dist_res_name.replace('\n','').replace(' ','_').replace('-','_').replace(',',''))
        
            #create reads per sequence graph
            reads_per_seq_list = sorted([chosen_clones_dict[clone]['reads_per_seq'] for clone in chosen_clones_dict],reverse=True)
            plot_reads_per_seq(output_path,reads_per_seq_list,population_reads_per_seq,rps_dist_name,rps_dist_file_name)
    
            #creat a sorted list (by expantion rate) of (clone,expantion,polarization,reads_per_seq)
            sorted_clones_by_expantion = sorted([(clone,chosen_clones_dict[clone]['expansion'],chosen_clones_dict[clone]['polarization'],chosen_clones_dict[clone]['reads_per_seq']) for clone in chosen_clones_dict], key=operator.itemgetter(1), reverse = True)
    
            #print log file of clones expantion-polarization
            with open(output_path + ex_pol_dist_name.replace('\n','').replace('-','_').replace(' ','_') +'.txt' , "w") as handle:
                handle.write('chain: ' + chain + '\n')
                if len(subsets):    
                    handle.write('subsets: ' + subsets_str[10:] + '\n')
                else:
                    handle.write('subsets: ' + subsets_str + '\n')
                handle.write('total number of clones: ' + str(len(chosen_clones_dict)) + '\n')
                handle.write('total number of sequences: ' + str(total_seqs) + '\n')
                handle.write('total number of reads: ' + str(total_reads) + '\n')
                handle.write('\n')
                handle.write('clones expantion-polarization data:' + '\n')
                for clone in sorted_clones_by_expantion:    
                    handle.write('clone: ' + clone[0] + '  expansion: ' + str(round(clone[1],3)) + '%  polarization: ' + str(round(clone[2],3)) + '% reads_per_seq: ' + str(round(clone[3],3)) + '\n')
            handle.close()
        
            #create an excel file
            clones_list = list(chosen_clones_dict.keys())
            xls_seqs = [chosen_clones_dict[clone]['seqs'] for clone in chosen_clones_dict]
            xls_reads = [chosen_clones_dict[clone]['reads'] for clone in chosen_clones_dict]
            xls_expansion = [chosen_clones_dict[clone]['expansion'] for clone in chosen_clones_dict]
            xls_polarization = [chosen_clones_dict[clone]['polarization'] for clone in chosen_clones_dict]
            xls_reads_per_seq = [chosen_clones_dict[clone]['reads_per_seq'] for clone in chosen_clones_dict]
        
            dictionary = {'Clones':clones_list,'sequences':xls_seqs,'reads':xls_reads,'expansion':xls_expansion,'polarization':xls_polarization,'reads_per_seq':xls_reads_per_seq}
            print_dict_to_xlsx(output_path,ex_pol_dist_name.replace('\n','').replace('-','_').replace(' ','_'),dictionary)
    
#        print(output_path)
            print('log file for all clones of desirable parameters values was created')
    
    else:
        print('no clones of specified parameters were found')
    
    return chosen_clones_dict





"""
get_sequences_by_clones_params
"""

#get_sequences_by_clones_params(father_path,'final_fasta_mouse_simulation.fasta','IGH',346596, expansion=[4,8], polarization=[4,8], reads_per_seq=[], clone_identifier = 'CDR3: ')
def get_sequences_by_clones_params(path,file_name,chain,*subsets, expansion=[0,100], polarization=[0,100], reads_per_seq=[], clone_identifier = 'clone: '):

    #define input folder
#    if not os.path.isdir(path + "source_sequences"):
#        print('source_sequences folder was not found') #print error if sources folder does not exist
#    input_path = path + "source_sequences/"
    input_path = path    

    #create output path
    if not os.path.isdir(path + 'clones_analysis_for_' +  file_name[:-6] + '/clones'):
        os.makedirs(path + 'clones_analysis_for_' +  file_name[:-6] + '/clones')
    output_path = path + 'clones_analysis_for_' +  file_name[:-6] + '/clones/' 
    
    chosen_clones_dict = get_clones_by_params(path,file_name,chain,*subsets, expansion=expansion, polarization=polarization, reads_per_seq=reads_per_seq, clone_identifier = clone_identifier, plot = False)

    if len(subsets):
        subsets_str = '_subsets_' + '_'.join([str(subset) for subset in subsets])
    else:
        subsets_str = '_all_subsets'

    if len(reads_per_seq):
        output_file_name = 'expansion_' + str(expansion[0]) + '_' + str(expansion[1]) + '_polarization_' + str(polarization[0]) + '_' + str(polarization[1]) + '_reads_per_seq_' + str(reads_per_seq[0]) + '_' + str(reads_per_seq[1]) + subsets_str + '_clones_sequences_from_' + file_name
    else:
        output_file_name = 'expansion_' + str(expansion[0]) + '_' + str(expansion[1]) + '_polarization_' + str(polarization[0]) + '_' + str(polarization[1]) + subsets_str + '_clones_sequences_from_' + file_name        

    sources_cnt = 0 
    with open(output_path + output_file_name , "w") as handle:
        for record in SeqIO.parse(open(input_path + file_name, "r"), "fasta"):
            if get_clone(str(record.description), identifier = clone_identifier) in chosen_clones_dict and any([any([x in [sub[0] for sub in get_reads(str(record.description))] for x in subsets]),len(subsets)==0]):
                sources_cnt += 1
                rec = SeqRecord(Seq(str(record.seq),generic_protein), id = str(record.id), description = str(record.description))
                SeqIO.write(rec, handle, "fasta")  

    
    print(str(len(chosen_clones_dict)) + ' clones of specified parameters')
    print(str(sources_cnt) + ' sequences for clones')
    handle.close()
    
    if sources_cnt:
        print('clone sequences fasta file was created')
    else:
        os.remove(output_path + output_file_name)

    return chosen_clones_dict


"""
get_clone_sequences
"""
#get_clone_sequences(files_dir,'final_fasta.fasta','CARSVSAYNAFDIW')
#get_clone_sequences(files_dir,'final_fasta_mouse.fasta','CARITTATGAMDYW')
def get_clone_sequences(path,file_name,clone,clone_identifier = 'CDR3: '):
    
        #define input folder
#    if not os.path.isdir(path + "source_sequences"):
#        print('source_sequences folder was not found') #print error if sources folder does not exist
#    input_path = path + "source_sequences/"
    input_path = path    

    #create output path
    if not os.path.isdir(path + 'clones_analysis_for_' +  file_name[:-6] + '/clones'):
        os.makedirs(path + 'clones_analysis_for_' +  file_name[:-6] + '/clones')
    output_path = path + 'clones_analysis_for_' +  file_name[:-6] + '/clones/' 
    
    sources_cnt = 0 
    with open(output_path + clone + '_clone_sequences_from_' + file_name , "w") as handle:
        for record in SeqIO.parse(open(input_path + file_name, "r"), "fasta"):
            if get_clone(str(record.description), identifier = clone_identifier) ==  clone:
                sources_cnt += 1
                rec = SeqRecord(Seq(str(record.seq),generic_protein), id = str(record.id), description = str(record.description))
                SeqIO.write(rec, handle, "fasta")  
    handle.close()
    
    print('\n')
    print('clone: ' + clone)
    print(str(sources_cnt) + ' sources for clone')
    if sources_cnt:
        print('\n')
        print(clone + ' sources fasta file was created')
    else:
        os.remove(output_path + clone + '_clone_sequences_from_' + file_name)
    
        
    
if __name__ == '__main__':

        
    if len(sys.argv) == 2:
        function = sys.argv[1]
    else:
        function = ''
        
    if function == 'get_clone_sequences':
        print('\n')
        path = str(input('Insert antibodies file path: ')).replace('\\','/') + '/'
        print('\n')
        file_name = str(input('Insert antibodies file name: '))
        
        print('\n')
        cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
        while all([cdr3_clonotyping!='y',cdr3_clonotyping!='n']):
            cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
            
        if cdr3_clonotyping == 'y':
            clone_identifier = 'CDR3: '
        else:
            print('\n')
            clone_identifier = str(input('Insert clone identifier (a string in record description before clone id):'))
        
        print('\n')
        clone = str(input('Insert clone: '))
        
        while clone is not '':
            get_clone_sequences(path,file_name,clone,clone_identifier=clone_identifier)
            print('\n\n')
            clone = str(input('Insert another clone or press Enter to finish: '))
            
    
    elif function == 'get_sequences_by_clones_params':
        print('\n')
        path = str(input('Insert antibodies file path: ')).replace('\\','/') + '/'
        print('\n')
        file_name = str(input('Insert antibodies file name: '))
        print('\n')
        chain = str(input('Insert chain: '))
        
        print('\n')
        subsets = input('Insert subset(s) seperated by space (press enter for all subsets):')
        if subsets == '':
            subsets = []
        else:
            subsets = [int(x) for x in subsets.split(' ')]
        
        print('\n')
        expansion = eval(input('Insert expansion rate boundaries ([a,b] where a is lower bound): '))
        print('\n')
        polarization = eval(input('Insert polarization rate boundaries ([a,b] where a is lower bound): '))
        print('\n')
        reads_per_seq = eval(input('Insert reads per sequence boundaries ([a,b] where a is lower bound, [] if all values wanted): '))
        
        print('\n')
        cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
        while all([cdr3_clonotyping!='y',cdr3_clonotyping!='n']):
            cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
            
        if cdr3_clonotyping == 'y':
            get_sequences_by_clones_params(path,file_name,chain,*subsets, expansion=expansion, polarization=polarization, reads_per_seq=reads_per_seq, clone_identifier = 'CDR3: ')
        else:
            print('\n')
            clone_identifier = str(input('Insert clone identifier (a string in record description before clone id):'))
            get_sequences_by_clones_params(path,file_name,chain,*subsets, expansion=expansion, polarization=polarization, reads_per_seq=reads_per_seq, clone_identifier=clone_identifier)
        
        
    elif function == 'get_clones_by_params':
        print('\n')
        path = str(input('Insert antibodies file path: ')).replace('\\','/') + '/'
        print('\n')
        file_name = str(input('Insert antibodies file name: '))
        print('\n')
        chain = str(input('Insert chain: '))
        
        print('\n')
        subsets = input('Insert subset(s) seperated by space (press enter for all subsets):')
        if subsets == '':
            subsets = []
        else:
            subsets = [int(x) for x in subsets.split(' ')]
        
        print('\n')
        expansion = eval(input('Insert expansion rate boundaries ([a,b] where a is lower bound): '))
        print('\n')
        polarization = eval(input('Insert polarization rate boundaries ([a,b] where a is lower bound): '))
        print('\n')
        reads_per_seq = eval(input('Insert reads per sequence boundaries ([a,b] where a is lower bound, [] if all values wanted): '))
        
        print('\n')
        cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
        while all([cdr3_clonotyping!='y',cdr3_clonotyping!='n']):
            cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
            
        if cdr3_clonotyping == 'y':
            get_clones_by_params(path,file_name,chain,*subsets, expansion=expansion, polarization=polarization, reads_per_seq=reads_per_seq, clone_identifier = 'CDR3: ')
        else:
            print('\n')
            clone_identifier = str(input('Insert clone identifier (a string in record description before clone id):'))
            get_clones_by_params(path,file_name,chain,*subsets, expansion=expansion, polarization=polarization, reads_per_seq=reads_per_seq, clone_identifier=clone_identifier)
        
    elif function == 'get_clones':
        
        print('\n')
        path = str(input('Insert antibodies file path: ')).replace('\\','/') + '/'
        print('\n')
        file_name = str(input('Insert antibodies file name: '))
        print('\n')
        chain = str(input('Insert chain: '))
       
        print('\n')
        subsets = input('Insert subset(s) seperated by space (press enter for all subsets):')
        if subsets == '':
            subsets = []
        else:
            subsets = [int(x) for x in subsets.split(' ')]
        
        print('\n')
        cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
        while all([cdr3_clonotyping!='y',cdr3_clonotyping!='n']):
            cdr3_clonotyping = str(input('clone by CDR3 sequences?[y/n] '))
            
        if cdr3_clonotyping == 'y':
            get_clones(path,file_name,chain,*subsets,clone_identifier='CDR3: ')
        else:
            print('\n')
            clone_identifier = str(input('Insert clone identifier (a string in record description before clone id):'))
            get_clones(path,file_name,chain,*subsets,clone_identifier=clone_identifier)

    else:
        print('\n')
        print('function is not recognized, please Insert a valid function name: "python clones_analysis.py funcion_name"')
        
        
    
    

            
    
            
            
            