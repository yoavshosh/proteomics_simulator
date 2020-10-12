# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:42:53 2019

@author: shosh

full process of proteomimcs simulator
1. preperation of input (based on coding sequences fasta file and editing sites table (in required format))
2. digestion of transcriptome and creation of all peptides (native and edited) adtabase
3. post processing of peptides database for genomic locations of editing sites and marking of "genomic informative" peptides
"""

import subprocess
import argparse
import os
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Proteomics Simulator Arguments')
run_parser = parser.add_argument_group('Run Editom Digestion')
#parameters for input preparation stage
run_parser.add_argument('-i_fasta', dest='fasta_file', action='store', required = True, help='Path to fasta file')
run_parser.add_argument('-i_sites', dest='sites_file', action='store', required = True, help='Path to fasta file')
run_parser.add_argument('-use_longest', dest='use_longest', action='store', default = 'False', help='True if output sould include only longest variant per each gene')
run_parser.add_argument('-refGene_file', dest='refGene_file', action='store', default = 'no_file', help='refGene_file with data regarding all genes and thier variants')
run_parser.add_argument('-o', dest='output_name', action='store', default = 'output', help='fasta output name')
run_parser.add_argument('-mm_types', dest='mm_types', action='store', nargs = '+' ,default = 'AG', help='mm types from anovar file to include in fasta file')
run_parser.add_argument('-first_is_met', dest='met_as_good_records', action='store',default = 'True', help='keep only records that start with metionin')
run_parser.add_argument('-stop_as_bad_records', dest='stop_as_bad_records', action='store',default = 'True', help='keep only records that do not contain stop codons')
run_parser.add_argument('-last_is_stop', dest='last_is_stop', action='store',default = 'True', help='keep only records that end with stop codon')

#parameters for digestion stage
run_parser.add_argument('-mc', '--missed_cleavages', dest='missed_cleavages', action='store', default = '2', help='Number of allowed missed cleavages per peptide')
run_parser.add_argument('-min_aa', dest='min_aa_number', action='store', default = '7', help='Minimal number of amino acids per peptide')
run_parser.add_argument('-max_mass', dest='max_mass', action='store', default = '4600', help='Maximal peptides mass')
run_parser.add_argument('-max_sites', dest='max_sites_per_pep', action='store', default = '20', help='maximal number of editing sites per peptide')
run_parser.add_argument('-print_peps_xlsx', dest='print_peps_xlsx', action='store', default = 'False', help='printing peptides tabel in xlsx format (not recommended for large data sets)')

#parameters for post processing stage
run_parser.add_argument('-check_genomic_inf', dest='check_genomic_inf', action='store', default = 'True', help='for each peptide check if edited sites are all from the same genomic location - in cases where few variants of the same gene are present in proteom and thus are not marked as informative in 2 stage')

arguments = parser.parse_args()

#defining paths
input_coding_mrna_path = '/'.join(arguments.fasta_file.split('/')[:-1]) + '/'
med_path = input_coding_mrna_path + 'proteomics_simulation/' #this path is created in the preparation script
param_str = '_'+str(arguments.missed_cleavages)+'mc_'+str(arguments.min_aa_number)+'minl_'+str(arguments.max_mass)+'maxm_'+str(arguments.max_sites_per_pep)+'maxes'
#creating ouput folder
if not os.path.exists(med_path + 'results_from_' + arguments.output_name + param_str+  '/'):
    os.makedirs(med_path + 'results_from_' + arguments.output_name + param_str+ '/')
output_path = med_path + 'results_from_' + arguments.output_name + param_str + '/'

#log file to which all stdout is written
logfile = open(med_path + arguments.output_name+'_log.txt','w')


input_prparation_script = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') + '/' + 'prepare_simulator_input_from_coding_seqs_and_sites_tbl.py'
input_prparation_cmd = 'python ' + input_prparation_script + ' ' + ' '.join(['-i_fasta',arguments.fasta_file,
                                                                             '-i_sites',arguments.sites_file,
                                                                             '-use_longest',arguments.use_longest,
                                                                             '-refGene_file', arguments.refGene_file,
                                                                             '-o', arguments.output_name,
                                                                             '-mm_types',' '.join(arguments.mm_types),
                                                                             '-first_is_met',arguments.met_as_good_records,
                                                                             '-stop_as_bad_records', arguments.stop_as_bad_records,
                                                                             '-last_is_stop', arguments.last_is_stop])

print('\nStage 1 - Preparation of fasta file for digestion:\n')
logfile.write('\nStage 1 - Preparation of fasta file for digestion:\n')
try:
    p_prep = subprocess.Popen(input_prparation_cmd, shell = True, universal_newlines = True, stdout=subprocess.PIPE)
    while True:
        line = p_prep.stdout.readline()
        if line == '' and p_prep.poll() is not None:
            break
        if line:
            print(line.strip())
            logfile.write(line.strip()+'\n')
    rc_prep = p_prep.poll()
#    exit(0)
except subprocess.CalledProcessError as e:
    print(e.output.decode())
    logfile.write(e.output.decode()+'\n')
    exit(1)


while True:
    if rc_prep is not None:
        cleave_script = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') + '/' + 'cleave_rna_seqs_as_prots.py'
        cleave_cmd = 'python ' + cleave_script + ' ' + ' '.join([med_path, arguments.output_name+'.fasta', arguments.missed_cleavages, arguments.min_aa_number, arguments.max_mass, arguments.max_sites_per_pep, arguments.print_peps_xlsx])
        print('\n\nStage 2 - Digestion of mrna sequences and peptides DB creation:\n')
        logfile.write('\n\nStage 2 - Digestion of mrna sequences and peptides DB creation:\n')
        try:
            p_cleave = subprocess.Popen(cleave_cmd, shell = True, universal_newlines = True, stdout=subprocess.PIPE)
            while True:
                line = p_cleave.stdout.readline()
                if line == '' and p_cleave.poll() is not None:
                    break
                if line:
                    print(line.strip())
                    logfile.write(line.strip()+'\n')
            rc_cleave = p_cleave.poll()
#            exit(0)
        except subprocess.CalledProcessError as e:
            print(e.output.decode())
            logfile.write(e.output.decode()+'\n')
            exit(1)
        break
    
 
while True:
    if rc_cleave is not None:       
        post_script = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') + '/' + 'post_processing_for_genomic_locations.py'
        post_cmd = 'python ' + post_script + ' ' + ' '.join([output_path,'peps_from_'+arguments.output_name+'.pickle',arguments.sites_file,arguments.check_genomic_inf])
        print('\n\nStage 3 - Post processing of peptides DB for genomic locations and marking of "genomic" informative peptides:\n')
        logfile.write('\n\nStage 3 - Post processing of peptides DB for genomic locations and marking of "genomic" informative peptides:\n')
        try:
            p_post = subprocess.Popen(post_cmd, shell = True, universal_newlines = True, stdout=subprocess.PIPE)
            while True:
                line = p_post.stdout.readline()
                if line == '' and p_post.poll() is not None:
                    break
                if line:
                    print(line.strip())
                    logfile.write(line.strip()+'\n')
            rc = p_post.poll()
            exit(0)
        except subprocess.CalledProcessError as e:
            print(e.output.decode())
            logfile.write(e.output.decode()+'\n')
            exit(1)
        break

logfile.close()
