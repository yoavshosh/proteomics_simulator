import subprocess
import argparse
import os
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Proteomics Simulator Arguments')
run_parser = parser.add_argument_group('Run Editom Digestion')
run_parser.add_argument('-i', '--input_path', dest='input_path', action='store', required = True, help='Path to fasta input directory')
run_parser.add_argument('-f', '--input_fasta', dest='input_fasta', action='store', required = True, help='Name of input fasta file')
run_parser.add_argument('-mc', '--missed_cleavages', dest='missed_cleavages', action='store', default = '2', help='Number of allowed missed cleavages per peptide')
run_parser.add_argument('-min_aa', dest='min_aa_number', action='store', default = '7', help='Minimal number of amino acids per peptide')
run_parser.add_argument('-max_mass', dest='max_mass', action='store', default = '4600', help='Maximal peptides mass')
run_parser.add_argument('-max_sites', dest='max_sites_per_pep', action='store', default = '20', help='maximal number of editing sites per peptide')
run_parser.add_argument('-print_peps_fasta', dest='print_peps_fasta', action='store', default = 'False', help='printing peptides tabel in fasta format if True (not recommended as files are usualy large and record description is not convenient to present peptides data')
run_parser.add_argument('-print_peps_xlsx', dest='print_peps_xlsx', action='store', default = 'False', help='printing peptides tabel in xlsx format (not recommended for large data sets)')
arguments = parser.parse_args()

#creating output folder
param_str = '_'+str(arguments.missed_cleavages)+'mc_'+str(arguments.min_aa_number)+'minl_'+str(arguments.max_mass)+'maxm_'+str(arguments.max_sites_per_pep)+'maxes'
#creating ouput folder
if not os.path.exists(arguments.input_path + 'results_from_' + arguments.input_fasta.split('.')[0] + param_str+  '/'):
    os.makedirs(arguments.input_path + 'results_from_' + arguments.input_fasta.split('.')[0] + param_str+ '/')
output_path = arguments.input_path + 'results_from_' + arguments.input_fasta.split('.')[0] + param_str + '/'


#assuming that main proteomics simulator script (cleave_rna_seqs_as_prots.py) is within the same directory as this file'
cleave_script = os.path.dirname(os.path.realpath(__file__)).replace('\\','/') + '/' + 'cleave_rna_seqs_as_prots.py'


#cleave_script = 'python C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/scripts/cleave_rna_seqs_as_prots.py'
cmd = 'python ' + cleave_script + ' ' + ' '.join([arguments.input_path, arguments.input_fasta, arguments.missed_cleavages, arguments.min_aa_number, arguments.max_mass, arguments.max_sites_per_pep, arguments.print_peps_fasta, arguments.print_peps_xlsx])
logfile = open(output_path + 'log.txt','w')

#run cleave_rna_seqs_as_prots.py with arguments
try:
    p = subprocess.Popen(cmd, shell = True, universal_newlines = True, stdout=subprocess.PIPE)
    while True:
        line = p.stdout.readline()
        if line == '' and p.poll() is not None:
            break
        if line:
            print(line.strip())
            logfile.write(line.strip()+'\n')
    rc = p.poll()
    exit(0)
except subprocess.CalledProcessError as e:
    print(e.output.decode())
    logfile.write(e.output.decode()+'\n')
    exit(1)

logfile.close()
