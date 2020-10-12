import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from CleaveFastaPeptidesFunctions import get_clone, get_chain, cleavage_sites


"""
create_isobaric_associated_peps
goes through peptides_from_fasta output file (containing peptides)
create a dictionary of all isobaric associated peptids (I\L in peptieds, and therefore could not be fully recognized)
output - a dictionary in that format
{pep:{igh_associated_peps:[pep1_id,pep2_id],igk_associated_peps:[],igl_associated_peps:[],ambigouos_chain_associated_peps:[]}}
"""
#apd = create_isobaric_associated_peps(files_dir + 'cleaved_peptides/','trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta')
def create_isobaric_associated_peps_dict(path, file_name):
     
    temp_associated_peps_dict = {}
    
    print('\n')
    print('creating isobaric-associations dictionary from file:')
    print(file_name)
    for record in SeqIO.parse(open(path + file_name, "r"), "fasta"):
        
        chains_list = []
        igh_pep = 0
        igk_pep = 0
        igl_pep = 0
        ambigouos_chain_pep = 0
        
        #consider only peptides that could cause an isobaric problem
        if any(substring in str(record.seq) for substring in ['I','L']):
        
            if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
                des = record.description[len(record.id):]
                data = eval(des)
            else:
                data = eval(record.description)
                
            isobaric_pep = str(record.seq).replace('I','X').replace('L','X')
       
        
            #create chains set for peptide        
            for source in data:
                chains_list.append(source['chain'])
            chains_set = set(chains_list)
            
            #chain type is unique
            if len(chains_set) == 1:
                if chains_list[0] == 'IGH':
                    igh_pep += 1
                if chains_list[0] == 'IGK':
                    igk_pep += 1
                if chains_list[0] == 'IGL':
                    igl_pep += 1
            else:
                ambigouos_chain_pep +=1
                
            
            #create/update associated_peps_dict
            if isobaric_pep in temp_associated_peps_dict:
                if igh_pep:
                    temp_associated_peps_dict[isobaric_pep]['igh_associated_peps'].append(str(record.id))
                if igk_pep:
                    temp_associated_peps_dict[isobaric_pep]['igk_associated_peps'].append(str(record.id))
                if igl_pep:
                    temp_associated_peps_dict[isobaric_pep]['igl_associated_peps'].append(str(record.id))
                if ambigouos_chain_pep:
                    temp_associated_peps_dict[isobaric_pep]['ambigouos_chain_associated_peps'].append(str(record.id))              
            else:
                if igh_pep:
                    temp_associated_peps_dict.update({isobaric_pep:{'igh_associated_peps':[str(record.id)],
                                                               'igk_associated_peps':[],
                                                               'igl_associated_peps':[],
                                                               'ambigouos_chain_associated_peps':[]}})
                if igk_pep:
                    temp_associated_peps_dict.update({isobaric_pep:{'igh_associated_peps':[],
                                                               'igk_associated_peps':[str(record.id)],
                                                               'igl_associated_peps':[],
                                                               'ambigouos_chain_associated_peps':[]}})  
                if igl_pep:
                    temp_associated_peps_dict.update({isobaric_pep:{'igh_associated_peps':[],
                                                               'igk_associated_peps':[],
                                                               'igl_associated_peps':[str(record.id)],
                                                               'ambigouos_chain_associated_peps':[]}})  
                if ambigouos_chain_pep:
                    temp_associated_peps_dict.update({isobaric_pep:{'igh_associated_peps':[],
                                                               'igk_associated_peps':[],
                                                               'igl_associated_peps':[],
                                                               'ambigouos_chain_associated_peps':[str(record.id)]}})
    
    #create updated_associated_peps_dict with peps that actually have associated peptides but themselves
    associated_peps_dict = {}
    for isobaric_pep in temp_associated_peps_dict:
        k = (len(temp_associated_peps_dict[isobaric_pep]['igh_associated_peps'])+
             len(temp_associated_peps_dict[isobaric_pep]['igk_associated_peps'])+
             len(temp_associated_peps_dict[isobaric_pep]['igl_associated_peps'])+
             len(temp_associated_peps_dict[isobaric_pep]['ambigouos_chain_associated_peps']))
        if k>1:
            associated_peps_dict.update({isobaric_pep:temp_associated_peps_dict[isobaric_pep]})
    

    print('isobaric associations created')
    print(str(len(associated_peps_dict)) + ' isobaric associations in file (association = 2 or more unique peptides that could not be distinguished due to an isobaric problem)') 
    
    del(temp_associated_peps_dict)
    
    return associated_peps_dict



"""
create_isobaric_peps_clone_dict
check for each peptide if could cause an isobaric problem (I\L in pep)
if so, add isobaric:clone to dictionary
output - dict containning only isobaric peptides that that have more than one real peptide and can relate to more than one clone

"""
#ipcd = create_isobaric_peps_clones_dict(files_dir + 'cleaved_peptides/','trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta')
#ipcd = create_isobaric_peps_clones_dict(files_dir + 'cleaved_peptides/','trypsin_1_miss_peptides_from_simulation.fasta')
def create_isobaric_peps_clones_dict(path, file_name):
    
    temp_isobaric_peps_clone_dict = {}
    isobaric_peps_clone_dict = {}
    
    print('\n')
    print('creating isobaric-peptide--clones dictionary from file:')
    print(file_name)
    for record in SeqIO.parse(open(path + file_name, "r"), "fasta"):
        clones_list = []
        
        if any(substring in str(record.seq) for substring in ['I','L']):
            #get data from record's description
            if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
                des = record.description[len(record.id):]
                data = eval(des)
            else:
                data = eval(record.description)
            
            #isobaric peptide is:
            isobaric_pep = str(record.seq).replace('I','X').replace('L','X')
            
            #create clones_list of all clones peptide relates to
            for source in data:
                clones_list.append(source['clone'])   
            clones_list = list(set(clones_list))  
            
            #if isobaric_pep already exists in temp_isobaric_peps_clone_dict, add the clones from the current peptide clones list and add one the associated peptide counter
            if isobaric_pep in temp_isobaric_peps_clone_dict:
                temp_isobaric_peps_clone_dict[isobaric_pep]['clones'] += clones_list
                temp_isobaric_peps_clone_dict[isobaric_pep]['clones'] = list(set(temp_isobaric_peps_clone_dict[isobaric_pep]['clones']))
                temp_isobaric_peps_clone_dict[isobaric_pep]['peptides_counter'] += 1
            else:
                temp_isobaric_peps_clone_dict.update({isobaric_pep:{'clones':clones_list,'peptides_counter':1}})
                
    #create a new isobaric_peps_clonen_dict containing only isobaric_peptides that relate to more than one clone
    for isobaric_pep in temp_isobaric_peps_clone_dict:
        if temp_isobaric_peps_clone_dict[isobaric_pep]['peptides_counter'] > 1:
            if len(temp_isobaric_peps_clone_dict[isobaric_pep]['clones']) > 1:
                isobaric_peps_clone_dict.update({isobaric_pep:temp_isobaric_peps_clone_dict[isobaric_pep]})
    
    del(temp_isobaric_peps_clone_dict)
    
    print('isobaric-peptide vs clones dictionary created')
    print(str(len(isobaric_peps_clone_dict)) + ' isobaric-peptides that relate to more than one clone were found')
    
    return isobaric_peps_clone_dict
    



"""
create_clones_dict
based on query sequences file, create a dict of all clones and 
for each clone - the value is list of sources that relate to it
"""
#clones_dict = create_clones_dict(files_dir + 'cleaved_peptides/chain_sorted_files/','light_chains_from_trypsin_1_miss_CDR3_from_simulation.fasta')
def create_clones_dict(path, sources_file, chains, clone_identifier = 'CDR3: '):
    
    
    flag_chain = False
    light_chains = ['IGK','IGL']
    heavy_chains = ['IGH']
    
    clones_dict = {}
    
    print('\n')
    print('creating clones dictionary from file: ' + sources_file)
    for record in SeqIO.parse(open(path + sources_file, "r"), "fasta"):
        flag_chain = False

        if chains == 'light':
            if get_chain(str(record.description)) in light_chains:
                flag_chain = True
        if chains == 'heavy':
            if get_chain(str(record.description)) in heavy_chains:
                flag_chain = True

        
        if flag_chain:
            clone = get_clone(str(record.description),identifier=clone_identifier)
        
            if clone in clones_dict:   
                clones_dict[clone] += 1
            else:
                clones_dict.update({clone:1}) 
    
    print(str(len(clones_dict)) + ' ' + chains + ' chain clones in file')
    
#    return clones_dict, all_source_list
    return clones_dict


"""
check_informative
based on sources_data, determine if peptide is informative
"""
def check_informative(peptide_id, sources_data,heavy_chain_query_inf_treshold=3):
    
    sources = 0
    clone_informative = 0
    coverage_informative = 0
    clones_list = []
          
    #create clones set - all clones that peptide relates to
    for source in sources_data:    
        sources +=1
        clones_list.append(source['clone']) 
    clones_set = set(clones_list)
            
    #pep relates to one clone
    if len(clones_set)==1 and 'unknown' not in clones_set:
                
        for source in sources_data:
            for overlap in source['overlaps']:
                if overlap['is_overlap']==1:
                    if source['chain'] == 'IGH': #for heavy chains, 3 first AA in CDR3 sequence are generic
                        if overlap['end'] > (heavy_chain_query_inf_treshold-1):
                            clone_informative += 1 #count in how many sources peptide overlap the CDR3 element - shold be all
                        else:
                            coverage_informative += 1
                    else:
                        clone_informative += 1
                else:
                    coverage_informative += 1 #count in how many sources peptide doesnt overlap CDR3 element 

                
        #clone informative - relates to only one clone sequence and overlaps it
        if clone_informative and str(clones_set) != "{'unknown'}":
            inf = 'clone_informative'
                    
        #coverage informative - relates to one CDR3 sequence and doesnt overlap it
        elif coverage_informative and str(clones_set) != "{'unknown'}":
            inf = 'coverage_informative'
                  
        #not informative at all - doesnt relate to any known CDR3
        if str(clones_set) == "{'unknown'}":
            inf = 0
                   
        #for QA - shouldnt enter condition at all   
        if clone_informative and coverage_informative:
            print('peptide ' + peptide_id + ' could be clone_informative and coverage_informative')
            
    #not informative at all - peptide relates to multiple CDR3 sequences          
    else:
        inf = 0
      

    return inf, sources




"""
CreateMultiCoverageFiles
operate on multiple files with CreateCoverageFile
files must be suplied in a list of dictionaries
each dictionaey has:
    query_file - file of query elements (created by CleaveFastaPeptides)
    pep_file - file of peptides (created by CleaveFastaPeptides)
    miss_cleavages - misscleavages allowed in CleaveFastaPeptides when files were created
    chain - the chains of all peptides (assuming CleaveFastaPeptides outputs went through SpltFileByChain)
"""
#cov_list = CreateMultiCoverageFiles(files_dir+"cleaved_peptides/", files_list1)
#sim_cov_list = CreateMultiCoverageFiles(files_dir+"cleaved_peptides/", sim_files_list)
def create_multi_coverage_files(path, files_list, heavy_chain_query_inf_treshold = 3):
    
    coverage_list = []
    
    for files_pack in files_list:
        
        print('Creating Peptides Library - Informative/Uninformative Peptides for file:')
        print(files_pack['pep_file'])
        c = create_coverage_file(path, files_pack['pep_file'], files_pack['query_file'], heavy_chain_query_inf_treshold = heavy_chain_query_inf_treshold)
        c.update({'chain':files_pack['chain'],'miss_cleavages':files_pack['miss_cleavages'],'file':files_pack['pep_file'],'query_file':files_pack['query_file']})
        coverage_list.append(c)
        
    return coverage_list 
        
        

"""
CreateCoverageFile
for a source file containing all peptides after in-silico cleavage by some protease(s)
create a file with those same peptides adding the folowwong information:
if peptide belong to many sources from the same clone (determined by number identical CDR3 element in all sources) - Is informative
if Informative:
    if overlap the CDR3 element - CDR3-informative
    if not:
        count the number of sources it belongs to in that clone

heavy_chain_query_inf_treshold (optional) - default = 3 for CDR3 query sequence of heavy chains where overlap is defined only from the 3ed position of CDR3 element as supplied in source file (with two upstream AA)
also in the case of heavy chains CDR3 - first 3 letters are generic and therfore heavy_chain_query_inf_treshold
defines the overlap for peptides and CDR3 sequence that would not mark the peptide as CDR3_informative
even if it relates to sources with that specific CDR3 element only
"""
#lt0 = CreateCoverageFile(files_dir+"cleaved_peptides/",'light_chains_from_trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta','light_chains_from_trypsin_0_miss_CDR3_from_all_subset_unique_sequences_non_single.fasta')
#ht0 = CreateCoverageFile(files_dir+"cleaved_peptides/",'heavy_chains_from_trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta','heavy_chains_from_trypsin_0_miss_CDR3_from_all_subset_unique_sequences_non_single.fasta')
#sl0 = CreateCoverageFile(files_dir+"cleaved_peptides/",'light_chains_from_trypsin_0_miss_aspen_0_miss_peptides_from_simulation.fasta','light_chains_from_trypsin_0_miss_aspen_0_miss_CDR3_from_simulation.fasta')
#h0 = CreateCoverageFile(files_dir+"cleaved_peptides/chain_sorted_files/",'heavy_chains_from_trypsin_0_miss_aspen_0_miss_peptides_from_simulation.fasta','heavy_chains_from_trypsin_0_miss_aspen_0_miss_CDR3_from_simulation.fasta')
#sh0 = CreateCoverageFile(files_dir+"cleaved_peptides/chain_sorted_files/",'heavy_chains_from_trypsin_0_miss_aspen_0_miss_peptides_from_simulation.fasta','heavy_chains_from_trypsin_0_miss_aspen_0_miss_CDR3_from_simulation.fasta')
def create_coverage_file(rule_name, rule, pep_path, pep_file, isobaric_peps_clone_dict, clones_dict, clone_identifier = 'CDR3: ', heavy_chain_query_inf_treshold = 3, mw_threshold = [], miss_cleavages = [], print_to_file = True):
    
#    #define output folder
#    if not os.path.isdir(pep_path + "final_peptides_files/" + rule_name + '/'):
#        os.makedirs(pep_path + "final_peptides_files/" + rule_name + '/') 
#    outpath = pep_path + "final_peptides_files/" + rule_name + '/'
      
    outpath = pep_path
    
    n = 0
    n_clone_inf_pep = 0
    n_coverage_inf_pep = 0
    n_uninf_pep = 0
    n_clone_isobaric_inf_pep = 0
    n_coverage_isobaric_inf_pep = 0
    clones_covered_by_clone_inf_pep_list = []
    clones_covered_by_coverage_inf_pep_list = []
    clones_covered_by_clone_isobaric_inf_pep_list = []
    clones_covered_by_coverage_isobaric_inf_pep_list = []
           
    #new file name:     
    if len(mw_threshold):
        new_file_str = str(mw_threshold[0]) + '_' + str(mw_threshold[1]) + '_mw_'
    else:
        new_file_str = ''
        
    if len(miss_cleavages):
        for miss in miss_cleavages:
            new_file_str = new_file_str + str(miss) + '_'
        new_file_str = new_file_str + 'mc_only_'
        
    new_file_name = new_file_str + 'final_'+ pep_file
#    new_file_name = new_file_str

    print('\n')
    print('checking for informative peptides in file: ' + pep_file)
    #create inf_lib file - library of all peptides - each is classified as clone\coverage informative or uninformative (0)  
    with open(outpath + new_file_name, "w") as handle:
        
        #check for each record description in pep_file if relate to sequences with same CDR3
        for record in SeqIO.parse(open(pep_path + pep_file, "r"), "fasta"):
            
            number_of_cs = len(cleavage_sites(str(record.seq),rule))
            
            #first - check if peptide is in mw_threshold
            if len(mw_threshold):
                if mw_threshold[0] <= ProteinAnalysis(str(record.seq)).molecular_weight() <=  mw_threshold[1]:
                    pep_mw_in_threshold = True
                else:
                    pep_mw_in_threshold = False
            else:
                pep_mw_in_threshold = True
             
            #second - check if peptides contain misscleavages needed
            if len(miss_cleavages):
                if number_of_cs in miss_cleavages:
                    pep_miss_cleavages_is_needed = True
                else:
                    pep_miss_cleavages_is_needed = False
            else:
                pep_miss_cleavages_is_needed = True
                
            #for peptides that are in mw_threshold and have desired number of miss_cleavages check if they are clone\coverage informative and 
            if pep_mw_in_threshold and pep_miss_cleavages_is_needed:
                n+=1
    
                if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
                    des = record.description[len(record.id):]
                else:
                    des = record.description
                data = eval(des)

                #check if peptide is informative
                inf, num_of_sources_from_clone = check_informative(str(record.id),data['sources_data'],heavy_chain_query_inf_treshold)
                #isobaric peptide is:
                isobaric_pep = str(record.seq).replace('I','X').replace('L','X')
                
                #clone informative - relates to one CDR3 sequence and overlap it 
                if inf == 'clone_informative':
                    if isobaric_pep in isobaric_peps_clone_dict: #peptide is only isobaric informative. clone is sorted to a clones_covered_by_clone_isobaric_inf_pep_list
                        inf = 'clone_isobaric_informative'
                        clone_coverage = round(num_of_sources_from_clone/clones_dict[data['sources_data'][0]['clone']],2)
                        n_clone_isobaric_inf_pep += 1
                        clones_covered_by_clone_isobaric_inf_pep_list.append(data['sources_data'][0]['clone'])
                    else:
                        clone_coverage = round(num_of_sources_from_clone/clones_dict[data['sources_data'][0]['clone']],2)
                        n_clone_inf_pep += 1
                        clones_covered_by_clone_inf_pep_list.append(data['sources_data'][0]['clone'])
                    
                #coverage informative - relates to one CDR3 sequence and doesnt overlap it 
                elif inf == 'coverage_informative':
                    if isobaric_pep in isobaric_peps_clone_dict: #peptide is only isobaric informative. clone is sorted to a clones_covered_by_coverage_isobaric_inf_pep_list
                        inf = 'coverage_isobaric_informative'
                        clone_coverage = round(num_of_sources_from_clone/clones_dict[data['sources_data'][0]['clone']],2)
                        n_coverage_isobaric_inf_pep += 1
                        clones_covered_by_coverage_isobaric_inf_pep_list.append(data['sources_data'][0]['clone'])
                    else:
                        clone_coverage = round(num_of_sources_from_clone/clones_dict[data['sources_data'][0]['clone']],2)
                        n_coverage_inf_pep += 1
                        clones_covered_by_coverage_inf_pep_list.append(data['sources_data'][0]['clone'])

                #not informative at all - peptide relates to multiple CDR3 sequences          
                else:
                    clone_coverage = 0
                    n_uninf_pep += 1
            
                #add information to data from description
                new_data = {'informative':inf,'clone_coverage':clone_coverage,'cleavage_sites':number_of_cs}
                new_data.update(data)
            
                #write peptide record to new file
                if print_to_file:
                    rec = SeqRecord(Seq(str(record.seq),generic_protein), id = str(record.id), description = str(new_data))
                    SeqIO.write(rec, handle, "fasta")
                    
                    
        handle.close()
        
        if not print_to_file:
            os.remove(outpath + new_file_name)
    
    total_clones = len(clones_dict)
#    del(clones_dict)
    
    print('\n')
    print('informaive peptides check completed')
    print('file creaeted: ' + new_file_name)
    if len(mw_threshold):
        print('molecular weight threshold for peptides in file: ' + str(mw_threshold[0]) + '-' + str(mw_threshold[1]) + ' Da')
    if len(miss_cleavages):
        print('number of miss cleavages for peptides in file: ' + str(miss_cleavages))
    
    #create sets of all sources, sources covered by qurey inf pep, sources covered only by coverage inf pep,  sources uncovered by inf pep
    clone_inf_pep_covered_clones_set = set(clones_covered_by_clone_inf_pep_list)
    coverage_inf_pep_covered_clones_set = set(clones_covered_by_coverage_inf_pep_list)
    clones_covered_by_clone_isobaric_inf_pep_set = set(clones_covered_by_clone_isobaric_inf_pep_list)
    clones_covered_by_coverage_isobaric_inf_pep_set = set(clones_covered_by_coverage_isobaric_inf_pep_list)
            
    #calculate number of clones from each coverage type group     
    n_clones_covered_by_clone_inf_pep = len(clone_inf_pep_covered_clones_set)    
    n_clones_covered_by_coverage_inf_pep_only = len(coverage_inf_pep_covered_clones_set) - len(coverage_inf_pep_covered_clones_set & clone_inf_pep_covered_clones_set)    
    unification_clone_coverage_clones = set.union(clone_inf_pep_covered_clones_set,coverage_inf_pep_covered_clones_set)
    n_clones_covered_by_clone_isobaric_inf_pep = len([x for x in clones_covered_by_clone_isobaric_inf_pep_set if x not in unification_clone_coverage_clones])
    unification_clone_coverage_clone_isobaric_clones = set.union(clone_inf_pep_covered_clones_set,coverage_inf_pep_covered_clones_set,clones_covered_by_clone_isobaric_inf_pep_set)
    n_clones_covered_by_coverage_isobaric_inf_pep = len([x for x in clones_covered_by_coverage_isobaric_inf_pep_set if x not in unification_clone_coverage_clone_isobaric_clones])
    n_uncovered_clones = total_clones - sum([n_clones_covered_by_clone_inf_pep,n_clones_covered_by_coverage_inf_pep_only,n_clones_covered_by_clone_isobaric_inf_pep,n_clones_covered_by_coverage_isobaric_inf_pep])

    #calculate rates
    if total_clones:
        clones_covered_by_clone_inf_pep_rate = round(100*n_clones_covered_by_clone_inf_pep/total_clones,2)
        clones_covered_by_coverage_inf_pep_rate = round(100*n_clones_covered_by_coverage_inf_pep_only/total_clones,2)
        clones_covered_by_clone_isobaric_inf_pep_rate = round(100*n_clones_covered_by_clone_isobaric_inf_pep/total_clones,2)
        clones_covered_by_coverage_isobaric_inf_pep_rate = round(100*n_clones_covered_by_coverage_isobaric_inf_pep/total_clones,2)
        uncovered_clones_rate = round(100*n_uncovered_clones/total_clones,2)
    else:
        clones_covered_by_clone_inf_pep_rate = 0
        clones_covered_by_coverage_inf_pep_rate = 0
        clones_covered_by_clone_isobaric_inf_pep_rate = 0
        clones_covered_by_coverage_isobaric_inf_pep_rate = 0
        uncovered_clones_rate = 0
        

    if n == 0:
        clone_inf_pep_rate = 'NaN'
        coverage_inf_pep_rate = 'NaN'
        clone_isobaric_inf_pep_rate = 'NanN'
        coverage_isobaric_inf_pep_rate = 'NaN'
        uninf_pep_rate = 'NaN'
    else:
        clone_inf_pep_rate = round(100*n_clone_inf_pep/n,2)
        coverage_inf_pep_rate = round(100*n_coverage_inf_pep/n,2)
        clone_isobaric_inf_pep_rate = round(100*n_clone_isobaric_inf_pep/n,2)
        coverage_isobaric_inf_pep_rate = round(100*n_coverage_isobaric_inf_pep/n,2)
        uninf_pep_rate = round(100*n_uninf_pep/n,2)
     
            
    print('total clones: ' + str(total_clones))
    print('clones covered by clone informative peptides: ' + str(clones_covered_by_clone_inf_pep_rate) + '%')
    print('clones covered by coverage informative peptides: ' + str(clones_covered_by_coverage_inf_pep_rate) + '%')
    print('clones covered by clone isobaric-informative peptides: ' + str(clones_covered_by_clone_isobaric_inf_pep_rate) + '%')
    print('clones covered by coverage isobaric-informative peptides: ' + str(clones_covered_by_coverage_isobaric_inf_pep_rate) + '%')
    print('uncovered clones: ' + str(uncovered_clones_rate) + '%')

    print('total peptieds: '+str(n))
    print('clone informative peptides: ' + str(clone_inf_pep_rate) + '%')
    print('coverage informative peptides: ' + str(coverage_inf_pep_rate) + '%')
    print('clone isobaric-informative peptides: ' + str(clone_isobaric_inf_pep_rate) + '%')
    print('coverage isobaric-informative peptides: ' + str(coverage_isobaric_inf_pep_rate) + '%')
    print('uninformative peptides: ' + str(uninf_pep_rate) + '%')
    
    #for QA - should print 100%
    if n != 0:
        print('success rate in classifying peptides: ' + str(round(clone_inf_pep_rate + coverage_inf_pep_rate + clone_isobaric_inf_pep_rate + coverage_isobaric_inf_pep_rate + uninf_pep_rate,0)) + '%')
    else:
        print('non of the peptides met the molecular-weight or miss-cleavages conditions')
    
    return({'file_name':new_file_name,
            'outpath':outpath,
            'total_peptides':n,
            'total_clones':total_clones,
            'clones_covered_by_clone_inf_pep':n_clones_covered_by_clone_inf_pep,
            'clones_covered_by_coverage_inf_pep':n_clones_covered_by_coverage_inf_pep_only,
            'clones_covered_by_clone_isobaric_inf_pep':n_clones_covered_by_clone_isobaric_inf_pep,
            'clones_covered_by_coverage_isobaric_inf_pep':n_clones_covered_by_coverage_isobaric_inf_pep,
            'uncovered_clones':n_uncovered_clones,
            'clone_inf_pep':n_clone_inf_pep,
            'coverage_inf_pep':n_coverage_inf_pep,
            'clone_isobaric_inf_pep':n_clone_isobaric_inf_pep,
            'coverage_isobaric_inf_pep':n_coverage_isobaric_inf_pep,
            'uninf_pep':n_uninf_pep,
            })
    
    


"""
SpltFileByChain
for files containing Lib of peptides as created by CreateLib
Seperate sequences by chains to different files
peptieds that relate to IGH sequences
peptieds that relate to IGK sequences
peptieds that relate to IGL sequences
peptides relate to IGK and IGL sequences
peptides relate to IGH and (IGK or IGL) sequences
"""

#i = SpltFileByChain(files_dir+"cleaved_peptides/", ['trypsin_0_miss_peptides_from_all_subset_unique_sequences_non_single.fasta','trypsin_1_miss_peptides_from_all_subset_unique_sequences_non_single.fasta','trypsin_2_miss_peptides_from_all_subset_unique_sequences_non_single.fasta'])
#sim_split = SpltFileByChain(files_dir+"cleaved_peptides/", ['trypsin_0_miss_peptides_from_simulation.fasta','trypsin_1_miss_peptides_from_simulation.fasta','trypsin_2_miss_peptides_from_simulation.fasta'])
#sim_split = SpltFileByChain(files_dir+"cleaved_peptides/",'trypsin_0_miss_aspen_0_miss_peptides_from_simulation.fasta')
#sim_split = SpltFileByChain(files_dir+"cleaved_peptides/",'trypsin_0_miss_aspen_0_miss_CDR3_from_simulation.fasta')
def split_file_by_chain(path, file, associated_peps_dict, inf_checked = False):
    
#    #define output folder inside input file folder
#    if not os.path.isdir(path + "chain_sorted_files/"):
#        os.makedirs(path + "chain_sorted_files/") 
#    outpath = path + "chain_sorted_files/"
    
    outpath = path
    
    print('\n')
    print('splitting_file:')
    print(file)
          
    n = 0
    IGH_cnt = 0
    IGK_cnt = 0
    IGL_cnt = 0
    HK_cnt = 0
    HL_cnt = 0
    KL_cnt = 0
    HKL_cnt = 0
    different_chains = 0
    sources_list = []
        
    #create files for partition of records based on chain type
    heavy_file = open(outpath+"heavy_chains_"+file, "w")
    light_file = open(outpath+"light_chains_"+file, "w")
    ambiguous_chain_file = open(outpath+"ambiguous_chains_"+file, "w")
        
    for record in SeqIO.parse(open(path + file, "r"), "fasta"):
            
        n+=1
        chains_list = []
            
        if record.description.startswith(record.id): #for some reason the description also contains the id and we need to eliminate it
            des = record.description[len(record.id):]
        else:
            des = record.description
        
        data = eval(des)
            
        #create a list of all chain types that peptide relate to (sources chains)
        if inf_checked:
            sources_data = data['sources_data']
        else:
            sources_data = data
        
        for source in sources_data:
            chains_list.append(source['chain'])
            sources_list.append(source['source_id'])
                    
        chains_set = set(chains_list)
        
        isobaric_pep = str(record.seq).replace('I','X').replace('L','X')
        
        #check for assos
        if isobaric_pep in associated_peps_dict:
            associated_peps = {}
            for chain in associated_peps_dict[isobaric_pep]:
                associated_peps_of_chain = [v for v in associated_peps_dict[isobaric_pep][chain] if v != str(record.id)]
                associated_peps.update({chain:associated_peps_of_chain})
        else:
            associated_peps = 'none'        
        
        
        new_description = {'isobaric_associations':associated_peps,'sources_data':sources_data}
        
        rec = SeqRecord(Seq(str(record.seq),generic_protein), id = str(record.id), description = str(new_description))
            
        #chain type is unique
        if len(chains_set) == 1:
                
            if chains_list[0] == 'IGH':
                IGH_cnt += 1
                SeqIO.write(rec, heavy_file, "fasta")
            elif chains_list[0] == 'IGK':
                IGK_cnt += 1
                SeqIO.write(rec, light_file, "fasta")
            elif chains_list[0] == 'IGL':
                IGL_cnt += 1
                SeqIO.write(rec, light_file, "fasta")
            else:
                SeqIO.write(rec, ambiguous_chain_file, "fasta")
                
        #chain type is not unique
        else:
            if 'IGH' in chains_list and 'IGK' in chains_list and 'IGL' not in chains_list:
                HK_cnt +=1
            elif 'IGH' in chains_list and 'IGL' in chains_list and 'IGK' not in chains_list:
                HL_cnt +=1
            elif 'IGL' in chains_list and 'IGK' in chains_list and 'IGH' not in chains_list:
                KL_cnt +=1
            elif 'IGL' in chains_list and 'IGK' in chains_list and 'IGH' in chains_list:
                HKL_cnt +=1

            SeqIO.write(rec, ambiguous_chain_file, "fasta")
                
    heavy_file.close()
    light_file.close()
    ambiguous_chain_file.close()
    sources = len(set(sources_list))
        
    print(str(sources) + ' total sources in source file')
    print(str(n) + ' total peptides in file')
    print(str(IGH_cnt) + ' IGH peptides (' + str(100*round(IGH_cnt/n,4)) + '%)')
    print(str(IGK_cnt) + ' IGK peptides (' + str(100*round(IGK_cnt/n,4)) + '%)')
    print(str(IGL_cnt) + ' IGL peptides(' + str(100*round(IGL_cnt/n,4)) + '%)')
    print(str(HK_cnt) + ' IGH + IGK peptides (' + str(100*round(HK_cnt/n,3)) + '%)')
    print(str(HL_cnt) + ' IGH + IGL peptides (' + str(100*round(HL_cnt/n,3)) + '%)')
    print(str(KL_cnt) + ' IGL + IGK peptides (' + str(100*round(KL_cnt/n,3)) + '%)')
    print(str(HKL_cnt) + ' IGH + IGL + IGK peptides (' + str(100*round(HKL_cnt/n,3)) + '%)')
    print(str(n - IGH_cnt - IGK_cnt - IGL_cnt - HK_cnt - HL_cnt - KL_cnt - HKL_cnt) + ' are neither of these groups')
        
    split_dict_list= {'unsplitted_file':file,
                      'outpath':outpath,
                      'heavy_chains_file':"heavy_chains_"+file,
                      'light_chains_file':"light_chains_"+file,
                      'ambiguous_file':"ambiguous_chains_"+file,
                      'peptides':n,
                      'sources':sources,
                      'IGH_cnt':IGH_cnt,
                      'IGK_cnt':IGK_cnt,
                      'IGL_cnt':IGL_cnt,
                      'HK_cnt':HK_cnt,
                      'HL_cnt':HL_cnt,
                      'KL_cnt':KL_cnt,
                      'HKL_cnt':HKL_cnt}         
    
    return split_dict_list
