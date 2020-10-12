import os
import sys
import DataProcessing
import Charts
from time import gmtime, strftime, localtime
from ScanFile import molecular_weight_list
from CleaveFastaPeptides import generate_cleaved_peptide_files
from Proteomics_Simulation_Params import *


#e,f,g,h = cleave_Igs(files_dir,'simulation.fasta',rules_list_alt,mw_threshold=[350,6000],miss_cleavages_filtration=[0,1])
def cleave_Igs(father_path, input_file, rules_list, clone_identifier = clone_identifier, continuation = continue_CDR3_element, mw_threshold = [], miss_cleavages_filtration = [], create_filtered_database=False):
    
    heavy_list = []
    heavy_threshold_misses_list = []
    light_list = []
    light_threshold_misses_list = []
    labels_for_cov_chart = []
    
    filter_peps = False
    if len(mw_threshold) or len(miss_cleavages_filtration):
        filter_peps = True
    
    #concatenate all misscleavag numbers needed for inf\uninf examinning
    misses_string = ''
    if len(miss_cleavages_filtration): 
        misses_string = '_mc_' + '_'.join([str(i) for i in miss_cleavages_filtration])
    
    mw_string = ''
    if len(mw_threshold):
        mw_string = '_mw_threshold_' + str(mw_threshold[0]) + '_' + str(mw_threshold[1])
        
    filter_string = '\nMW threshold: ' + str(mw_threshold[0])+'-'+str(mw_threshold[1]) + ', ' + misses_string.replace('_mc_','miss cleavages: ').replace('_',', ')
    if filter_string.endswith(', '):
        filter_string = filter_string[:-2]
    
    for rules in rules_list:
        
        rules_string = ''
        for rule in rules:
            rules_string = rules_string + rule['name'] + '-' + str(rule['miss_cleavages']) + 'mc_' 
        rules_string = rules_string[:-1]
        labels_for_cov_chart.append(rules_string)
        
        
        #cleave source file and create queries file (CDR3 sequences) and peptides file
        cleaved_peptides_files = generate_cleaved_peptide_files(rules_string.replace(' ','_').replace('-','_'),father_path,input_file,rules, clone_identifier = clone_identifier, continuation = continue_CDR3_element)
        
        #split peptides file to heavy/light chains files
        associated_peps_dict = DataProcessing.create_isobaric_associated_peps_dict(cleaved_peptides_files['outpath'],cleaved_peptides_files['peptides_file'])
        splited_pep_file = DataProcessing.split_file_by_chain(cleaved_peptides_files['outpath'],cleaved_peptides_files['peptides_file'],associated_peps_dict)
        del(associated_peps_dict)
        
        isobaric_peps_clones_dict = DataProcessing.create_isobaric_peps_clones_dict(cleaved_peptides_files['outpath'],cleaved_peptides_files['peptides_file'])
        os.remove(cleaved_peptides_files['outpath']+cleaved_peptides_files['peptides_file'])
        
        #first - crate a dict from query_file of all clones exists and number sources they belong to
        clones_dict = DataProcessing.create_clones_dict(father_path,input_file,'heavy',clone_identifier = clone_identifier)
        
        #create informative/non informative peptides file from heavy chains
        inf_heavy = DataProcessing.create_coverage_file(rules_string.replace(' ','_').replace('-','_'), rules, splited_pep_file['outpath'], splited_pep_file['heavy_chains_file'], isobaric_peps_clones_dict, clones_dict, clone_identifier = clone_identifier)
        
        #create db general stats
        Charts.plot_db_sizes(inf_heavy['outpath'], [cleaved_peptides_files['total_sequences'],cleaved_peptides_files['total_peptides'],cleaved_peptides_files['unique_peptides']], 
                      ['number of sources','total peptides generated','unique peptides'], 'Peptides DB size with ' + rules_string, 'number of sequences [*10^6]', 'db_creation_for_'+rules_string.replace(' ','_'))
        #create db chain types peptides stats
        Charts.plot_chains_stats(inf_heavy['outpath'], 
                                [splited_pep_file['IGH_cnt'],splited_pep_file['IGK_cnt'],splited_pep_file['IGL_cnt'],splited_pep_file['HK_cnt'],splited_pep_file['HL_cnt'],splited_pep_file['KL_cnt'],splited_pep_file['HKL_cnt']], 
                                [cleaved_peptides_files['IGH_seqs'],cleaved_peptides_files['IGK_seqs'],cleaved_peptides_files['IGL_seqs'],0,0,0,0],cleaved_peptides_files['total_sequences'],
                                ['IGH','IGK','IGL','IGH\IGK','IGH\IGL','IGK\IGL','IGH\IGK\IGL'], 
                                'Peptides Chains Statistics with ' + rules_string, 'Percentage [%]', 'peptides_chains_statistics_with_'+rules_string.replace(' ','_'))  
        #create peptides classes stats
        Charts.plot_percentages(inf_heavy['outpath'], [inf_heavy['clone_inf_pep'],inf_heavy['coverage_inf_pep'],inf_heavy['clone_isobaric_inf_pep'],inf_heavy['coverage_isobaric_inf_pep'],inf_heavy['uninf_pep']], 
                      ['clone-inf','coverage-inf','clone isobaric-inf','coverage isobaric-inf','uninformative'], 'Heavy Chains Peptides Classification\nProtease: ' + rules_string, 'Percentage [%]', 'heavy_chains_peptides_classes_for_'+rules_string.replace(' ','_'))
        
        if filter_peps:
            #create informative/non informative peptides file from heavy chains for threshold peptides only
            inf_heavy_threshold_misses = DataProcessing.create_coverage_file(rules_string.replace(' ','_').replace('-','_'), rules, splited_pep_file['outpath'], splited_pep_file['heavy_chains_file'],  isobaric_peps_clones_dict, clones_dict, clone_identifier = clone_identifier, mw_threshold = mw_threshold, miss_cleavages = miss_cleavages_filtration,print_to_file=create_filtered_database)
            #create peptides classes stats
            Charts.plot_percentages(inf_heavy_threshold_misses['outpath'], [inf_heavy_threshold_misses['clone_inf_pep'],inf_heavy_threshold_misses['coverage_inf_pep'],inf_heavy_threshold_misses['clone_isobaric_inf_pep'],
                                    inf_heavy_threshold_misses['coverage_isobaric_inf_pep'],inf_heavy_threshold_misses['uninf_pep']],['clone-inf','coverage-inf','clone isobaric-inf','coverage isobaric-inf','uninformative'],
                                    '\nHeavy Chains Peptides Classification\nProtease: ' + rules_string + filter_string,  
                                    'Percentage [%]',
                                    'heavy_chains_peptides_classes_for'+mw_string.replace('_threshold','') + misses_string +'_'+ rules_string.replace(' ','_'))
            
        os.remove(splited_pep_file['outpath']+splited_pep_file['heavy_chains_file'])    
        
        
        del(clones_dict)
        mw_list = molecular_weight_list(inf_heavy['outpath'],inf_heavy['file_name'],groups = ['clone_informative','coverage_informative'])
        #plot molecular weight distribution (w\o threshold)
        mw_dist_name = 'Informative Peptides from Heavy Chains: Molcular Weight Distribution\nProtease: ' + rules_string
        if len(mw_list):
            Charts.plot_mw_hist(inf_heavy['outpath'],mw_list,'heavy_mw_'+ rules_string.replace(' ','_').replace('-','_'),
                                mw_dist_name, 'Molecular Weight [Da]', 'Frequency',threshold=mw_threshold)
        else:
            print('no informative peptides for heavy chains using ' + rules_string)
            print('Molecular weight distribution was not created')
        del(mw_list)
        
        
        #first - crate a dict from query_file of all clones exists and number sources they belong to
        clones_dict = DataProcessing.create_clones_dict(father_path,input_file,'light',clone_identifier = clone_identifier)
        
        #create informative/non informative peptides file from light chains
        inf_light = DataProcessing.create_coverage_file(rules_string.replace(' ','_').replace('-','_'), rules, splited_pep_file['outpath'], splited_pep_file['light_chains_file'], isobaric_peps_clones_dict, clones_dict, clone_identifier = clone_identifier)
        
        #create peptides classes stats
        Charts.plot_percentages(inf_light['outpath'], [inf_light['clone_inf_pep'],inf_light['coverage_inf_pep'],inf_light['clone_isobaric_inf_pep'],inf_light['coverage_isobaric_inf_pep'],inf_light['uninf_pep']], 
                      ['clone-inf','coverage-inf','clone isobaric-inf','coverage isobaric-inf','uninformative'], 'Light Chains Peptides Classification\nProtease: ' + rules_string, 'Percentage [%]', 'light_chains_peptides_classes_for_'+rules_string.replace(' ','_'))
        if filter_peps:
            #create informative/non informative peptides file from light chains for threshold peptides only
            inf_light_threshold_misses = DataProcessing.create_coverage_file(rules_string.replace(' ','_').replace('-','_'), rules, splited_pep_file['outpath'], splited_pep_file['light_chains_file'], isobaric_peps_clones_dict, clones_dict, clone_identifier = clone_identifier, mw_threshold = mw_threshold, miss_cleavages = miss_cleavages_filtration,print_to_file=create_filtered_database)                  
            #create peptides classes stats
            Charts.plot_percentages(inf_light_threshold_misses['outpath'], [inf_light_threshold_misses['clone_inf_pep'],inf_light_threshold_misses['coverage_inf_pep'],inf_light_threshold_misses['clone_isobaric_inf_pep'],
                                    inf_light_threshold_misses['coverage_isobaric_inf_pep'],inf_light_threshold_misses['uninf_pep']],['clone-inf','coverage-inf','clone isobaric-inf','coverage isobaric-inf','uninformative'],
                                    '\nLight Chains Peptides Classification\nProtease: ' + rules_string + filter_string,  
                                    'Percentage [%]',
                                    'light_chains_peptides_classes_for'+mw_string.replace('_threshold','') + misses_string +'_'+ rules_string.replace(' ','_'))
        
        os.remove(splited_pep_file['outpath']+splited_pep_file['light_chains_file'])  
    
        del(clones_dict)
        mw_list = molecular_weight_list(inf_light['outpath'],inf_light['file_name'],groups = ['clone_informative','coverage_informative'])
        #plot molecular weight distribution (w\o threshold)
        mw_dist_name = 'Informative Peptides from Light Chains: Molcular Weight Distribution\nProtease: ' + rules_string
        if len(mw_list):
            Charts.plot_mw_hist(inf_light['outpath'],mw_list,'light_mw_'+ rules_string.replace(' ','_').replace('-','_'),
                                mw_dist_name, 'Molecular Weight [Da]', 'Frequency',threshold=mw_threshold) 
        else:
            print('no informative peptides for light chains using ' + rules_string)
            print('Molecular weight distribution was not created')
        del(mw_list)
    
        del(isobaric_peps_clones_dict)           
        
        heavy_list.append(inf_heavy)
        inf_heavy_dict = {'parameter':list(inf_heavy.keys()),'value':list(inf_heavy.values())}
        Charts.print_dict_to_xlsx(inf_heavy['outpath'],'final_heavy_chains_peptides_stats',inf_heavy_dict)
        light_list.append(inf_light)
        inf_light_dict = {'parameter':list(inf_light.keys()),'value':list(inf_light.values())}
        Charts.print_dict_to_xlsx(inf_light['outpath'],'final_light_chains_peptides_stats',inf_light_dict)
        del(inf_heavy_dict,inf_light_dict)
        
        if filter_peps:
            heavy_threshold_misses_list.append(inf_heavy_threshold_misses)
            inf_heavy_threshold_misses_dict = {'parameter':list(inf_heavy_threshold_misses.keys()),'value':list(inf_heavy_threshold_misses.values())}
            Charts.print_dict_to_xlsx(inf_heavy_threshold_misses['outpath'],'final_filtered_heavy_chains_peptides_stats',inf_heavy_threshold_misses_dict)
            light_threshold_misses_list.append(inf_light_threshold_misses)
            inf_light_threshold_misses_dict = {'parameter':list(inf_light_threshold_misses.keys()),'value':list(inf_light_threshold_misses.values())}
            Charts.print_dict_to_xlsx(inf_light_threshold_misses['outpath'],'final_filtered_light_chains_peptides_stats',inf_light_threshold_misses_dict)
            del(inf_heavy_threshold_misses_dict,inf_light_threshold_misses_dict)
        
    
    #define coverage output folder
    if not os.path.isdir(father_path + "/cleaved_peptides_from_" + input_file[:-6] + "/clones_coverage"):
        os.makedirs(father_path + "/cleaved_peptides_from_" + input_file[:-6] + "/clones_coverage") 
    clone_coverage_path = father_path + "/cleaved_peptides_from_" + input_file[:-6] + "/clones_coverage/"
    
    #plot coverage rates
    clone_inf = []
    coverage_inf = []
    for x in heavy_list:
        if x['total_clones']:    
            clone_inf.append(round(100*x['clones_covered_by_clone_inf_pep']/x['total_clones'],2))
            coverage_inf.append(round(100*x['clones_covered_by_coverage_inf_pep']/x['total_clones'],2))
        else:
            clone_inf.append(0)
            coverage_inf.append(0) 
    Charts.plot_coverage_rate(clone_coverage_path, clone_inf, coverage_inf, labels_for_cov_chart,
                       'Clones Coverage Potential - Heavy Chains Clones', 'Protease Used','Coverage Percentage (%)',
                       strftime("%Y%m%d%H%M%S", localtime()) + '_heavy_clones_coverage')
    
    clone_inf = []
    coverage_inf = []
    for x in light_list:
        if x['total_clones']:
            clone_inf.append(round(100*x['clones_covered_by_clone_inf_pep']/x['total_clones'],2))
            coverage_inf.append(round(100*x['clones_covered_by_coverage_inf_pep']/x['total_clones'],2))
        else:
            clone_inf.append(0)
            coverage_inf.append(0)
    Charts.plot_coverage_rate(clone_coverage_path, clone_inf, coverage_inf, labels_for_cov_chart,
                       'Clones Coverage Potential - Light Chains Clones', 'Protease Used','Coverage Percentage (%)',
                       strftime("%Y%m%d%H%M%S", localtime()) + '_light_clones_coverage')
    
    #plot coverage rates of filtered data-base
    if filter_peps:
        
        clone_inf = []
        coverage_inf = []
        for x in heavy_threshold_misses_list:
            if x['total_clones']:    
                clone_inf.append(round(100*x['clones_covered_by_clone_inf_pep']/x['total_clones'],2))
                coverage_inf.append(round(100*x['clones_covered_by_coverage_inf_pep']/x['total_clones'],2))
            else:
                clone_inf.append(0)
                coverage_inf.append(0)                
        Charts.plot_coverage_rate(clone_coverage_path, clone_inf, coverage_inf, labels_for_cov_chart,
                                  'Clones Coverage Potential - Heavy Chains Clones\n(by peptides within MW threshold: ' +str(mw_threshold[0])+'-'+str(mw_threshold[1])+misses_string.replace('_mc_',', miss cleavages: ').replace('_',', ')+')\n', 'Protease Used','Coverage Percentage (%)',
                                  strftime("%Y%m%d%H%M%S", localtime()) + '_heavy_clones_coverage_filtered' + mw_string + misses_string)
        clone_inf = []
        coverage_inf = []
        for x in light_threshold_misses_list:
            if x['total_clones']:    
                clone_inf.append(round(100*x['clones_covered_by_clone_inf_pep']/x['total_clones'],2))
                coverage_inf.append(round(100*x['clones_covered_by_coverage_inf_pep']/x['total_clones'],2))
            else:
                clone_inf.append(0)
                coverage_inf.append(0) 
        Charts.plot_coverage_rate(clone_coverage_path, clone_inf, coverage_inf, labels_for_cov_chart,
                                  'Clones Coverage Potential - Light Chains Clones\n(by peptides within MW threshold: ' +str(mw_threshold[0])+'-'+str(mw_threshold[1])+misses_string.replace('_mc_',', miss cleavages: ').replace('_',', ')+')\n', 'Protease Used','Coverage Percentage (%)',
                                  strftime("%Y%m%d%H%M%S", localtime()) + '_light_clones_coverage_filtered' + mw_string + misses_string)
    
    
    return heavy_list,heavy_threshold_misses_list,light_list,light_threshold_misses_list


if __name__ == '__main__':
    
    new_father_path = "%r"%father_path
    new_father_path = father_path.replace('\\','/')
    if new_father_path[-1] != '/':
        new_father_path = new_father_path + '/'
    cleave_Igs(new_father_path, input_file, rules_list, mw_threshold = mw_threshold_filtration, miss_cleavages_filtration = miss_cleavages_filtration,create_filtered_database=create_filtered_database)
    
    print('\n')
    print('proteomics simulation completed at ' + strftime("%Y/%m/%d, %H:%M:%S", localtime()) + ', for ' + str(len(rules_list)) + ' proteases (or proteases combinations)')

    
    