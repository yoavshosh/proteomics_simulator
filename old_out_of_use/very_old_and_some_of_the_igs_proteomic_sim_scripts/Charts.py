import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import statistics
from pandas import DataFrame, ExcelWriter
from matplotlib.ticker import FuncFormatter
import os

"""
print_dict_to_xls
print dictionary of lists to excel. key as field name and its list is printes to the column bellow it.
arryas must be same length
"""
#print_dict_to_xlsx('C:/Users/user/Downloads/','exceltest',dicttest)
def print_dict_to_xlsx(path,file_name,dictionary,sheet = 'sheet1',writer = ''):
        
    #create output path
    if not os.path.isdir(path + '/excel_files/'):
        os.makedirs(path + '/excel_files/')
    output_path = path + '/excel_files/'
    
    if writer == '':
    #print data frame to excel
        DataFrame(dictionary).to_excel(output_path + file_name + '.xlsx', sheet_name=sheet, index=False)
    else:
        DataFrame(dictionary).to_excel(writer, sheet_name=sheet, index=False)
    

"""
to_percent
format y axis as percent
"""
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'

"""
to_millions
format y axis as M
"""    
def to_millions(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(y/10)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\M$'
    else:
        return s + 'M'



"""
plot_reads_per_seq
plot a graph of all clones reads per seq (given as a list) 
and a reference line representing average reads per clone (given as a float)
""" 

#plot_reads_per_seq(father_path,dtest[1],3,'Reads_per_Clone_from_Bla','bfs')
def plot_reads_per_seq(path,reads_per_seq,population_reads_per_seq,dist_name,dist_file_name):
    
    #creating x axos in which every int is a clone
    x = list(range(0,len(reads_per_seq)))
    
    plt.scatter(x,reads_per_seq, marker='.', color = 'black')
    
    #determinnig clone closest to average reads_per_clone
#    closest_clone = min(range(len(x)), key=lambda i: abs(reads_per_seq[i]-population_reads_per_seq))
#    if reads_per_seq[closest_clone] < population_reads_per_seq:
#        closest_clone = closest_clone-1
    
    #setting average reads per sequence line and number of clones with higher than average RPS:
    plt.axhline(population_reads_per_seq, color='blue', linestyle='--')
    plt.annotate('average RPS = ' + str(round(population_reads_per_seq,2)), xy=(statistics.median(x), population_reads_per_seq+0.7),xytext=(statistics.median(x), population_reads_per_seq+0.5))

    plt.xlim(0,max(x))
    plt.ylim(0,max(reads_per_seq)+1)
    plt.title(dist_name)
    plt.ylabel('Reads per Sequence')
    plt.xlabel('Clone Index')
    plt.savefig(path + dist_file_name + '.png', bbox_inches='tight', dpi=500)
    plt.close('all')
    
    
#    dictionary = {'index':x,'reads_per_seq':reads_per_seq}
#    print_dict_to_xlsx(path,dist_file_name,dictionary)


"""
plot_expansion_polarization
plot a scatter graph of all clones as a point in two dim space
xaxis = clone expantion
yaxis = clone polarization
"""


def plot_expansion_polarization_rps_residuals(path,clone_exp_pol_dict,dist_name,dist_file_name):
    
    clones_exp_pol_list = []
    for clone in clone_exp_pol_dict:
        clones_exp_pol_list.append((clone_exp_pol_dict[clone]['expansion'],clone_exp_pol_dict[clone]['polarization']- clone_exp_pol_dict[clone]['expansion']))

# =============================================================================
#     print(str(len(clones_exp_pol_list)) + 'clones')
#     print(str(sum(x[1]>0 for x in clones_exp_pol_list)) + ' above RPS line')
#     print(str(sum(x[1]==0 for x in clones_exp_pol_list)) + ' on RPS line')
#     print(str(sum(x[1]<0 for x in clones_exp_pol_list)) + ' under RPS line')
# =============================================================================
    
    arrs = list(zip(*clones_exp_pol_list))
    plt.scatter(arrs[0],arrs[1], marker='.', color = 'black')
    plt.xlim(0,max(arrs[0])+0.1)
    plt.ylim(min(arrs[1])-0.1,max(arrs[1])+0.1)
    plt.title(dist_name)
    plt.ylabel('Polarization residual from average RPS line')
    plt.xlabel('Expansion [% from total sequences]')
    plt.axhline(0, color='blue', linestyle='--')
    plt.savefig(path + dist_file_name + '.png', bbox_inches='tight', dpi=500)
    plt.close('all')
    



#plot_expansion_polarization(files_dir,dtest,'dist_name')
def plot_expansion_polarization(path,clone_exp_pol_dict,dist_name,dist_file_name):
    
    clones_exp_pol_list = []
    for clone in clone_exp_pol_dict:
        clones_exp_pol_list.append((clone_exp_pol_dict[clone]['expansion'],clone_exp_pol_dict[clone]['polarization']))

    arrs = list(zip(*clones_exp_pol_list))
    
    plt.scatter(arrs[0],arrs[1], marker='.', color = 'black')
    plt.xlim(0,max(arrs[0])+0.1)
    plt.ylim(0,max(arrs[1])+0.1)
    plt.title(dist_name)
    plt.ylabel('Polarization [% from total sequences reads]')
    plt.xlabel('Expansion [% from total sequences]')
    ax = plt.axes()
    ax.plot([0, 1], [0, 1], linestyle = '--')
    plt.savefig(path + dist_file_name + '.png', bbox_inches='tight', dpi=500)
        
#    plt.show()
    plt.close('all')

#    dictionary = {'expansion':list(arrs[0]),'polarization':list(arrs[1])}
#    print_dict_to_xlsx(path,dist_file_name,dictionary)



"""
plot_mw_hist
given an array of molecular weights - plot hist of molecular weight frequency
also - the function takes list of two values as range of realistic mol-weights for proteomics and add the threshold lines to the chart 
"""
#a,b = plot_mw_hist(father_path, mw_list, 'dist_name_mw', 'title', 'xlabel', 'ylabel')
#plot_mw_hist('C:/Users/user/Google_Drive/NGS_project/Proteomics_project/', mw_list, 'dist_name_mw', 'title', 'xlabel', 'ylabel')
def plot_mw_hist(path, mw_list, dist_name, title, xlabel, ylabel, bins = 100, threshold = []):

    
    frequencies = []
    bins_borders = list(np.linspace(0,max(mw_list),bins+1))
        
    for i in range(0, len(bins_borders)-1):
        frequencies.append(len([x for x in mw_list if bins_borders[i] < x <= bins_borders[i+1]])/len(mw_list))
    
    x = np.delete(bins_borders,0)

    width = -(x[1]-x[0])
        
    plt.figure()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.bar(x, frequencies, width, color='blue', align = 'edge')

    #set xticks based on the distance between threshold boundaries
    if len(threshold)==2:
        xticks_list = [i for i in threshold]
        add_tick = False
        if threshold[1]-threshold[0]>1500:
            middle = (threshold[1] + threshold[0])/2
            maxdif = abs(threshold[1]-threshold[0])
            for i in x:
                if abs(i-middle)<maxdif:
                    add_tick = True
                    tick_to_add = i
                    maxdif = abs(tick_to_add-middle)
            if add_tick:
                xticks_list.append(tick_to_add)
                add_tick = False
                theotick = threshold[1] + tick_to_add
                mindif = 100
                for i in [e for e in x if abs(e-theotick)<100]:
                    if abs(i-theotick)<mindif:
                        mindif = abs(i-theotick)
                        add_tick = True
                        tick_to_add = i
                if add_tick:
                    xticks_list.append(tick_to_add)
        
#        if all([abs(x[30]-i)>1000 for i in threshold]):
#            xticks_list = xticks_list + [x[30]]
#        if all([abs(x[99]-i)>1000 for i in threshold]):
#            xticks_list = xticks_list + [x[99]]        
    else:
        xticks_list = [x[0],x[20],x[40],x[60],x[80],x[99]]    
        
    plt.xticks(xticks_list)
    plt.annotate('bin width = ' + str(round(-width,2)),xy=(0.7, 0.95),xycoords='axes fraction')
    plt.tight_layout()
    
    
    if len(threshold):
        for j in threshold:
            plt.axvline(j, color='black', linestyle='--')
    
    plt.savefig(path + dist_name, bbox_inches='tight', dpi=500)
    plt.close('all')
#    print(sum(frequencies))

    print('\n')
    print('molecular weight distribution created:' + dist_name)

    writer = ExcelWriter(path + '/excel_files/' + dist_name.replace(' ','_').replace('-','_') + '.xlsx')
    Molecular_weight_range = [str(i*x[0]) + '-' + str(x[i]) for i in range(len(x))]
    dictionary = {'Molecular_weight_range':Molecular_weight_range,'frequencies':frequencies}
     
    print_dict_to_xlsx(path,dist_name.replace(' ','_').replace('-','_'),dictionary,writer = writer,sheet = 'MW_histogram_100_bins')
    print_dict_to_xlsx(path,dist_name.replace(' ','_').replace('-','_'),{'Molecular Weight':mw_list},writer = writer,sheet = 'all_peptides_list')
    del(writer)
    
    return x, frequencies
    
"""
plot_coverage_rate
for multiple Fasta file Generated by CreateCoverageFile
save a stacked bar charts for clone coverage

input - two lists of percentage one
        one list of lables
        title, xlabel, ylabel
""" 
#plot_coverage_rate('C:/Users/user/Google_Drive/', list1, list2, labels_list, 'Clones Coverage Potential - Heavy Chains Clones\n(by peptides within MW threshold: 350-6000 Da)', 'Protease Used', 'Coverage Percentage (%)', 'clones_coverage_mw_threshold')
def plot_coverage_rate(outpath, list1, list2, label_list, title, xlabel, ylabel, dist_name):
        
    mpl_fig = plt.figure(figsize=(10,6))
    ax = mpl_fig.add_subplot(111)
    
    N = len(list1)
    ind = np.arange(N)
    width = 0.35
    
    p1 = plt.bar(ind, list1, width, color=(0.2,0.3,0.6))
    p2 = plt.bar(ind, list2, width, color=(0.6,0.8,0.5), bottom = list1)
    
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylim(1,100)
    
    plt.xticks(ind, label_list)
    plt.ylim([0,110])
    
    for i, v in enumerate([sum(x) for x in zip(*[list1,list2])]):
#        print(str(i) +' ' +  str(v))
        bar_val = str(v)[:4] + '%'
        ax.text(i-0.1, v+1, bar_val, color='black', fontweight='bold')

    lgd = plt.legend((p1[0], p2[0]), ('clones covered by\nclone Informative peptides', 'clones covered by\ncoverage Informative peptides'),
                     bbox_to_anchor=(1.35, 1), loc='upper right')
#    ax.legend(loc='lower right', bbox_to_anchor=(0, 1))
#    plt.show()
    mpl_fig.savefig(outpath + dist_name + '.png',bbox_extra_artists=(lgd,), bbox_inches='tight')
#    plt.show()
    
    plt.close('all')
    print('\n')
    print('clones coverage rates chart created: ' + dist_name)
    
    dictionary = {'cleavage rule':label_list,'clones covered by clone-inf peptides':list1,'clones covered by coverage-inf peptides':list2}
    print_dict_to_xlsx(outpath,dist_name.replace(' ','_').replace('-','_'),dictionary)



"""
plot_db_sizes
plot a bar chrats of three bars
number of unique sources
number of peptides created from sources
number of unique peptides
"""
#plot_db_sizes('C:/Users/user/Google_Drive/',[m,l,k],['unique sources','total peptides\n generated','unique peptides'],'Peptides DB from Trypsin Simulation (1 miss cleavages allowes)\nGeneral Statistics','Number of Sequences [Millions]','db_general_stats')
def plot_db_sizes(outpath, numbers_list, label_list, title, ylabel, dist_name):
       
    numbers_list_millions = [x/1000000 for x in numbers_list]
    
    mpl_fig = plt.figure()
    ax = mpl_fig.add_subplot(111)
    
    N = len(numbers_list_millions)
    ind = np.arange(N)
    width = 0.35
    
    plt.bar(ind, numbers_list_millions, width, color=(0.2,0.3,0.6))
    
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    plt.xticks(ind, label_list)
    plt.ylim([0,max(numbers_list_millions)*1.12])
    
    for i,v in enumerate(numbers_list_millions):
        ax.text(i+width-0.54,v + max(numbers_list_millions)*0.03, str(int(v*1000000)), color='black', fontweight='bold')
    
#    plt.show()

#    ax.legend(loc='lower right', bbox_to_anchor=(0, 1))
#    plt.show()
    mpl_fig.savefig(outpath + dist_name + '.png', bbox_inches='tight')
    print('\n')
    print('data-base general statistic created: ' + dist_name)
    
    plt.close('all')
    
    dictionary = {'Peptides':label_list,'number of peptides':numbers_list}
    print_dict_to_xlsx(outpath,dist_name.replace(' ','_').replace('-','_'),dictionary)



#plot_percentages('C:/Users/user/Google_Drive/', [3,6,2,7,8], ['a','b','c','d','e'], 'titel', 'y', 'bla')
def plot_percentages(outpath, numbers_list, label_list, title, ylabel, dist_name):
       
    tot_peptides = sum(numbers_list)
    
    try:
        percentage_list = [100*x/tot_peptides for x in numbers_list]
    except ZeroDivisionError:
        percentage_list = [0 for x in numbers_list]
    
    
    mpl_fig = plt.figure()
    ax = mpl_fig.add_subplot(111)
    
    N = len(percentage_list)
    ind = np.arange(N)
    width = 0.35
    
    plt.bar(ind, percentage_list, width, color=(0.2,0.3,0.6))
    
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    plt.xticks(ind, label_list)
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylim([0,105])
    
    for i,v in enumerate(percentage_list):
        ax.text(i+width-0.6,v+1.5, str(round(v,2))+'%', color='black', fontweight='bold')


    mpl_fig.savefig(outpath + dist_name + '.png', bbox_inches='tight')
    print('\n')
    print(dist_name + ' was crerated')
#    plt.show()
    plt.close('all')

    dictionary = {'Peptides':label_list,'frequency [%]':percentage_list}
    print_dict_to_xlsx(outpath,dist_name.replace(' ','_').replace('-','_'),dictionary)
    
    
    
#plot_chain_stats('C:/Users/user/Google_Drive/', [25.18,72.94,1.82,0.05,0,0,01,0,01], [723941,2947157,16485,0,0,0,0], 3687583, ['IGH','IGK','IGL','IGH\IGK','IGH\IGL','IGK\IGL','IGH\IGK\IGL'], 'Peptides Chains Statistics with trypsin-1mc', 'Percentage [%]', 'bla')
def plot_chains_stats(outpath, peptides_list, chains_list, total_sequences, label_list, title, ylabel, dist_name):
       
    tot_peptides = sum(peptides_list)
    
    try:
        percentage_list1 = [100*x/tot_peptides for x in peptides_list]
    except ZeroDivisionError:
        percentage_list1 = [0 for x in numbers_list]
    
    try:
        percentage_list2 = [100*x/total_sequences for x in chains_list]
    except ZeroDivisionError:
        percentage_list2 = [0 for x in chains_list]    
    
    
    mpl_fig = plt.figure()
    ax = mpl_fig.add_subplot(111)

    N = len(percentage_list1)
    ind = np.arange(N)
    width = 0.35
    
    line1 = plt.bar(ind, percentage_list1, width, color=(0.1,0.3,1))
    line2 = plt.plot(ind,percentage_list2, 'rv')
#    mpl_fig.ledgend(loc='upper right')
#    seqs.legend([line1, line2], ['peptides', 'sequences'], loc='upper right')
    
    plt.legend((line2[0],line1[0]), ('sequences','peptides'), loc = 'upper right')
    
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    plt.xticks(ind, label_list)
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylim([0,105])
    
    for i,v in enumerate(percentage_list1):
        ax.text(i+width-0.6,v+1.5, str(round(v,2))+'%', color='black', fontweight='bold')


    mpl_fig.savefig(outpath + dist_name + '.png', bbox_inches='tight')
    print('\n')
    print(dist_name + ' was crerated')
#    plt.show()
    plt.close('all')

    dictionary = {'Chain':label_list,'sequences frequency [%]':percentage_list2,'peptides frequency [%]':percentage_list1}
    print_dict_to_xlsx(outpath,dist_name.replace(' ','_').replace('-','_'),dictionary)

    
# =============================================================================
#     
# #plot_db_sizes2('C:/Users/user/Google_Drive/',[700436,7289034,1019825],['unique sources','total peptides\n generated','unique peptides'],'Peptides DB size with lys-c-1mc','Number of Sequences [Millions]','db_general_stats1')
# def plot_db_sizes2(outpath, numbers_list, label_list, title, ylabel, dist_name):
#        
#     numbers_list_millions = [x/1000000 for x in numbers_list]
#     
#     mpl_fig = plt.figure()
#     ax = mpl_fig.add_subplot(111)
#     
#     N = len(numbers_list_millions)
#     ind = np.arange(N)
#     width = 0.35
#     
#     plt.bar(ind, numbers_list_millions, width, color=(0.2,0.3,0.6))
#     
#     ax.set_ylabel(ylabel)
#     ax.set_title(title)
#     
#     plt.xticks(ind, label_list)
#     plt.ylim([0,7.289034*1.12])
#     
#     for i,v in enumerate(numbers_list_millions):
#         ax.text(i+width-0.54,v + max(numbers_list_millions)*0.03, str(int(v*1000000)), color='black', fontweight='bold')
#         
#     mpl_fig.savefig(outpath + dist_name + '.png', bbox_inches='tight')
# =============================================================================
