# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:03:46 2019

@author: shosh
This script runs MQ program for raw files in a directory.
It does so sequentialy - for each given proteome for which MQ should run against.
"""

import argparse
import subprocess
import time


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Preparing mqpar file and running maxquant for several proteoms sequentialy')
    run_parser = parser.add_argument_group('Run mqpar prep and run MQ')
    run_parser.add_argument('-raw_path', dest='raw_files_path', action='store', required = True, help='Path to raw input directory')
    run_parser.add_argument('-proteoms', dest='proteoms', action='store', required = True, help='list of proteomes to run maxquant against. if two or more files should be included in the same run, pass them separated by coma.')
    run_parser.add_argument('-fdr', dest='fdr', action='store', default = '0.01', help='FDR for peptides detection - for now it is the fdr of psmFdrCrosslink\proteins\peptides\sites detection')
    run_parser.add_argument('-min_aa', dest='min_pep_len', action='store', default = '7', help='Minimal number of amino acids per peptide')
    run_parser.add_argument('-max_mass', dest='max_pep_mass', action='store', default = '4600', help='Maximal peptides mass')
    run_parser.add_argument('-max_mc', dest='max_mc', action='store', default = '2', help='maximal number of missed cleavages per peptide')
    run_parser.add_argument('-quantification', dest='quantification', action='store', default = 'lfq', help='quantification method used in MS/MS')
    run_parser.add_argument('-phospho', dest='phospho', action='store', default = 'False', help='raw files are of phosphoproteomics')
    run_parser.add_argument('-enz', dest='enzymes', action='store', nargs='+', default = ['trypsin'], help='enzymes used for cleavage for peptides search. type unspecific for unspecific')
    run_parser.add_argument('-o', dest='output_names', action='store', default = 'ag_editings non_ag_editings random_ag_editings', help='names seperated by spaces: mqpar/log files/combined directory (MQ output) names for each proteom run. should be a list of size equal to number of proteoms passed. each name will be assigned to each proteom respectively')
    run_parser.add_argument('-run_mq', dest='run_maxquant', action='store', default = 'True', help='run MaxQuant right after mqpar preparation')
    arguments = parser.parse_args()
    
    raw_files_path = arguments.raw_files_path
    proteoms = arguments.proteoms.split(' ')
    fdr = arguments.fdr
    min_pep_len = arguments.min_pep_len
    max_pep_mass = arguments.max_pep_mass
    max_mc = arguments.max_mc
    output_names = arguments.output_names.split(' ')
    run_maxquant = arguments.run_maxquant
    enzymes = arguments.enzymes
    quantification = arguments.quantification
    phospho = arguments.phospho
    
    assert len(proteoms)==len(output_names), str(len(proteoms))+" proteoms passed. "+str(len(output_names))+" output names passed\neach proteom should have its own unique name"
    assert len(set(output_names))==len(output_names), str(len(proteoms))+"pass a unique output name for each run"

    for i,prot in enumerate(proteoms):
        
        zip_when_done = 'False'
        if i==len(proteoms)-1:
            zip_when_done = 'True'
        
        fasta_files = prot.split(',')
        fasta_files_str = ''
        for f in fasta_files:
            fasta_files_str += f+' '
        fasta_files_str = fasta_files_str.rstrip()
        
        enzymes_str = ''
        for e in enzymes:
            enzymes_str+= e+' '
        enzymes_str = enzymes_str.rstrip()
             
        cmd = 'python ~/scripts/proteomics_simulator/additional_modules/prepare_mqpar_file_for_maxquant_analysis.py ' + '-raw_path '+raw_files_path + ' -fasta_files '+fasta_files_str + ' -fdr '+fdr + ' -min_aa '+min_pep_len + ' -max_mass '+max_pep_mass + ' -max_mc '+max_mc + ' -enz '+enzymes_str + ' -quantification '+quantification + ' -phospho '+phospho + ' -o '+output_names[i] + ' -run_mq '+run_maxquant + ' -delet_intermediate_stages '+'True' +  ' -zip_when_done '+zip_when_done
        print('\n\n\n======================================================================================\n\n\n')
        print('Run name: ' + output_names[i])
        print(cmd)
        print('\nPreparing parameters and running MQ for ' + output_names[i])
        print('proteome fasta file are:')
        for f in fasta_files:
            print(f)
        print('\n')
        time.sleep(10)
        p = subprocess.Popen(cmd, shell = True, universal_newlines = True, stdout=subprocess.PIPE)
        while p.poll() is None:
            time.sleep(60*60)
#            print(output_names[i] + ' MQ analysis has not finished yet')
#        print(output_names[i] + ' MQ analysis has finished/terminated')
        while True:
            line = p.stdout.readline()
            if not line:
                break
            print (line.rstrip())
    
    print('\n\n\n======================================================================================\n\n\n')
    print('Sequential Running finished')
    exit(0)
            
        

    