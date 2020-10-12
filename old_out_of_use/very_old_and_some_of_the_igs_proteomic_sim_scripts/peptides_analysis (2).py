import sys
import inspect
import importlib
import os
import subprocess
import logging
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import generic_protein
from CleaveFastaPeptidesFunctions import get_clone

# peptide('C:/Proteomics_project/files/source_sequences/','C:/Proteomics_project/files/source_sequences/cleaved_peptides_from_simulation/chain_sorted_files/final_peptides_files/trypsin_0mc/','simulation.fasta','informative_checked_from_heavy_chains_from_trypsin_0_miss_peptides.fasta','NTLFLQMNNLR')
class peptide:
    
    def __init__(self,peptides_path,peptides_file,peptide_seq,sources_path = 'not_specified',sources_file = 'not_specified'):
        
        self.pep = peptide_seq
        self.sources_path = sources_path
        self.sources_file = sources_file
        self.peptides_path = peptides_path
        self.peptides_file = peptides_file
        self.sources_list = []
        
        
    """
    get_peptide_sources
    for peptide generated
    get all sources from source file into a fasta file
    """
    def get_peptide_sources(self):
        
        #create output path
        if not os.path.isdir(self.peptides_path + '/chosen_peptides'):
            os.makedirs(self.peptides_path + '/chosen_peptides')
        output_path = self.peptides_path + 'chosen_peptides/' 
        
        for record in SeqIO.parse(open(self.peptides_path + self.peptides_file, "r"), "fasta"):
            if str(record.seq) == self.pep:
                if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
                    des = record.description[len(record.id):]
                    data = eval(des)
                else:
                    data = eval(record.description)
                
                for source in data['sources_data']:
                    self.sources_list.append(source['source_id'])
    
        with open(output_path + self.pep + '_peptide_sources.fasta' , "w") as handle:
            for record in SeqIO.parse(open(self.sources_path + self.sources_file, "r"), "fasta"):
                if record.id in self.sources_list:
                    rec = SeqRecord(Seq(str(record.seq),generic_protein), id = str(record.id), description = str(record.description))
                    SeqIO.write(rec, handle, "fasta")  
    
        handle.close()
        
        if len(self.sources_list):
            print('\nfile created: ' + self.pep + '_peptide_sources.fasta')
        else:
            os.remove(output_path + self.pep + '_peptide_sources.fasta')


    """
    read_peptide_data
    for a given peptide in a given output file
    prints all relevant data for this peptide including all peptide sources
    """
    def read_peptide_data(self):
        
        #create output path
        if not os.path.isdir(self.peptides_path + '/chosen_peptides'):
            os.makedirs(self.peptides_path + '/chosen_peptides')
        output_path = self.peptides_path + 'chosen_peptides/' 
    
            
        pep_found = 0
    
        print('\n')
        for record in SeqIO.parse(open(self.peptides_path + self.peptides_file, "r"), "fasta"):
            if str(record.seq) == self.pep:
                pep_found = 1
                
                fn = output_path + self.pep + '_peptide_data.txt'
                os.remove(fn) if os.path.exists(fn) else None
                
                level = logging.INFO
                format = '%(message)s'
                handlers = [logging.FileHandler(fn), logging.StreamHandler()]
                logging.basicConfig(level = level, format = format, handlers = handlers)
            
                if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
                    des = record.description[len(record.id):]
                    data = eval(des)
                else:
                    data = eval(record.description)
        
                logging.info('peptide id: ' + str(record.id))
                logging.info('length: ' + str(len(record.seq)))
                logging.info('molecular weight: ' + str(round(ProteinAnalysis(str(record.seq)).molecular_weight(),3)))
                try:
                    logging.info('number of missed cleavages: ' + str(data['cleavage_sites']))
                except KeyError:
                    pass
                
                logging.info('number of sources: ' + str(len(data['sources_data'])))
                try:
                    logging.info('informative: ' + str(data['informative']))
                    logging.info('clone coverage: ' + str(round(100*data['clone_coverage'],2)) + '%')
                except KeyError:
                    pass
    
                logging.info('\n')
                if data['isobaric_associations'] == 'none':
                    logging.info('isobaric associations: none' )
                else:
                    print('isobaric associations:' )
                    for chain in data['isobaric_associations']:
                        if data['isobaric_associations'][chain]:    
                            logging.info(chain + ':')
                            [logging.info(str(x)) for x in data['isobaric_associations'][chain]]
                        
                logging.info('\n')
                logging.info('sources data:')
                logging.info('-------------')
                for source in data['sources_data']:
                    logging.info('\n')
                    logging.info('source id: ' + str(source['source_id']) + ', ' + 'chain: ' + str(source['chain'])
                          + ', ' + 'isotype: ' + str(source['isotype']) + ', ' + 'reads: ' + str(source['reads'])
                          + ', ' + 'clone: ' + str(source['clone'])
                          + ', ' + 'CDR3-P: ' + str(source['query_seq']))
                    logging.info('overlaps of peptide with cdr3-p:')
                
                    is_overlap = 0
                    for overlap in source['overlaps']:
                        if overlap['is_overlap']:
                            is_overlap = 1
                            logging.info('overlap position relative to cdr3-p: ')
                            logging.info('beginnig: ' + str(overlap['beginning']) + ', ' + 'end: ' + str(overlap['end']))
                            logging.info('peptide coverage of cdr3-p: ' + str(round(100*min(overlap['coverage'],1.0),2)) + '%')
                    if not is_overlap:
                        logging.info('no overlaps of peptide and cdr3-p')
                        
                logging.shutdown()
        
        if not pep_found:
            print('peptide was not found')
            
    
def get_clone_peptides(peptides_path, peptides_file, clone):
    
    inf_list = []
    peptides = 0
    
    #create output path
    if not os.path.isdir(peptides_path + '/chosen_peptides'):
        os.makedirs(peptides_path + '/chosen_peptides')
    output_path = peptides_path + 'chosen_peptides/' 
    
    with open(output_path + clone + '_clone_peptides.fasta' , "w") as handle:
        for record in SeqIO.parse(open(peptides_path + peptides_file, "r"), "fasta"): 

            if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
                data = eval(record.description[len(record.id):])
            else:
                data = eval(record.description)
    
            clones_list = [source['clone'] for source in data['sources_data']]
            
            if clone in clones_list:
                peptides += 1
                inf_list.append(data['informative'])
                rec = SeqRecord(Seq(str(record.seq),generic_protein), id = str(record.id), description = str(record.description))
                SeqIO.write(rec, handle, "fasta")
    
    pep_cnt_dict = Counter(inf_list)
    
    print('\n')
    print(str(len(inf_list)) + ' unique peptides found for clone ' + clone + ' :')
    for key in pep_cnt_dict:
        print(str(pep_cnt_dict[key]) + ' ' + str(key).replace('_',' ').replace('0','uninformative') + ' peptides')
        
    print('\n')
    if peptides:
        print('file created in peptides path (in chosen peptides directory): ' + clone + '_clone_peptides.fasta')
    else:
        os.remove(output_path + clone + '_clone_peptides.fasta')
        print('no peptides were found for clone ' + clone)
            


if __name__ == '__main__':
    
    if len(sys.argv) == 2:
        function = sys.argv[1]
    else:
        function = ''
        
    if function == 'get_peptide_sources':
        print('\n')
        sources_path = str(input('Insert antibodies file path: ')).replace('\\','/') + '/'
        print('\n')
        sources_file = str(input('Insert antibodies file name: '))
        print('\n')
        peptides_path = str(input('Insert peptides file path: ')).replace('\\','/') + '/'
        print('\n')
        peptides_file = str(input('Insert peptides file name: '))
        print('\n')
        peptide_seq = str(input('Insert peptide sequences: '))
        print('\n')
        while peptide_seq is not '':
            x = peptide(peptides_path,peptides_file,peptide_seq,sources_path=sources_path,sources_file=sources_file)
            x.read_peptide_data()
            x.get_peptide_sources()
            print('\n\n')
            peptide_seq = str(input('Insert another peptide sequences or press Enter to finish: '))
        print('\n')
        
        
    elif function == 'peptide_data':
        print('\n')
        peptides_path = str(input('Insert peptides file path: ')).replace('\\','/') + '/'
        print('\n')
        peptides_file = str(input('Insert peptides file name: '))
        print('\n')
        peptide_seq = str(input('Insert peptide sequences: '))
        print('\n')
        while peptide_seq is not '':
            x = peptide(peptides_path,peptides_file,peptide_seq)
            x.read_peptide_data()
            print('\n\n')
            peptide_seq = str(input('Insert another peptide sequences or press Enter to finish: '))
        print('\n')
        
        
    elif function == 'get_clone_peptides':
        print('\n')
        peptides_path = str(input('Insert peptides file path: ')).replace('\\','/') + '/'
        print('\n')
        peptides_file = str(input('Insert peptides file name: '))
        print('\n')
        clone = str(input('Insert clone: '))
        print('\n')
        while clone is not '':
            get_clone_peptides(peptides_path,peptides_file,clone)
            print('\n\n')
            clone = str(input('Insert another clone ID or press Enter to finish: '))
        print('\n')
        
    else:
        print('\n')
        print('function is not recognized, please Insert a valid function name: "python peptydes_analysis.py funcion_name"')
    





