import re
import os
import sys
import pandas as pd
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statsmodels.stats.multitest as p_adjust
from collections import Counter
from pylab import text
from scipy import stats
from collections import deque
from functools import reduce
from matplotlib import colors
from matplotlib.colors import LogNorm
from heapq import nsmallest
from Bio import SeqIO
import xlsxwriter

all_mm = ['AG','AC','AT','CA','CG','CT','GA','GC','GT','TA','TG','TC']


def get_sites_list(df):
    
    sites = []
    for sample in list(set(df['sample_name'])):
        samples_df = df[df['sample_name'] == sample]
        row_gene_sites = []
        for i, row in samples_df.iterrows():
            for mm_type in eval(row['genomic_keys']):
                row_gene_sites = row_gene_sites + [row['gene_name']+';'+s for s in mm_type]
        sites = sites + list(set(row_gene_sites))
    return sites


def plot_editings_sites_repetitions(fig_path,names, sites, groups = [1,2,3,4,5,6,7,8,9,10]):
    col = ['b','g','salmon','y','m','tan','r','silver','yellow','crimson','teal','k','orange','brown','limegreen','navy','gold','c','olive']
    
    bars = []
    sites_counts = [list(Counter(s).values()) for s in sites]
    for s in sites_counts:
        b=[]
        for g in sorted(groups[:-1]):
            b.append(len([c for c in s if c==g]))
        b.append(len([c for c in s if c>=groups[-1]]))
        bars.append(b)
    
    # set width of bar
    barWidth = 0.2
    
    groups_labels = []
    for i in groups:
        if i < groups[-1]:
            groups_labels.append(str(i))
    groups_labels.append(str(groups[-1])+'<=')
    

    # Set position of bar on X axis
    rs = [np.arange(1,len(bars[0])+1)]
    r_temp = rs[0]
    for i in names[1:]:
        r_temp = [x + barWidth for x in r_temp]
        rs.append(np.array(r_temp))

     
    # Make the plot
    for i in range(len(bars)):
        plt.bar(list(rs[i]), bars[i], color=col[i], width=barWidth, edgecolor='white', label=names[i])
      
    # Add xticks on the middle of the group bars
    plt.xlabel('Repetitions', fontweight='bold')
    plt.ylabel('Sites', fontweight='bold')
    plt.xticks([g+0.2 for g in groups], groups_labels) 
    plt.yticks(np.arange(0,20,5))
    plt.title('Discovered Editing Sites')
    
    # Create legend & Show graphic
    plt.legend()
    plt.savefig(fig_path)
#    plt.show()
    plt.close()


def scatter_discovered_peptides(df, name, path):
    color_markers = [('y', '+'), ('r', '*'), ('b', '>'), ('g', '<'), ('m', '^'), ('tan', 'D'), ('silver', '.'), ('yellow', 'o'), ('crimson', '<'), ('teal', '>'), ('k', 'D'), ('orange', '^'), ('brown', 'D'), ('limegreen', 'D'), ('navy', '^'), ('gold', 'o'), ('c', '.'), ('olive', '.')]
    
    tissues = Counter(list(df['tissue']))
    
    for i,t in enumerate(tissues.items()):
        print( t[0]+' '+str(t[1]))
        tissues_dfs = df[df['tissue']==t[0]]
        plt.scatter(list(tissues_dfs['total_peptides']),list(tissues_dfs['edited_peptides']), label = t[0]+' '+str(t[1]), color = color_markers[i][0], marker = color_markers[i][1])
    plt.legend(loc='center left', bbox_to_anchor=(0.8, 0.5))
    plt.title('Peptides for Samples - ' + name + ' proteom')
    plt.xlabel('total peptides')
    plt.savefig(path)
    plt.close()
    
    
if __name__ == '__main__':
    
#    summ1 = 'E:/RNA_editing_Large_files/MQ/results_from_ag_finished.txt'
    summ2 = 'E:/RNA_editing_Large_files/MQ/results_from_all_non_ag_finished.txt'
    summ3 = 'E:/RNA_editing_Large_files/MQ/results_from_random_ag_finished.txt'
    summ_files = [summ1,summ2,summ3]
    
    peps1 = 'E:/RNA_editing_Large_files/MQ/peptides_lists_from_ag_finished.txt'
    peps2 = 'E:/RNA_editing_Large_files/MQ/peptides_lists_from_all_non_ag_finished.txt'
    peps3 = 'E:/RNA_editing_Large_files/MQ/peptides_lists_from_random_ag_finished.txt'
    peps_files = [peps1,peps2,peps3]
    
    names = ['AG','NonAG','randomAG']
    
    
#    summ1 = 'E:/RNA_editing_Large_files/MQ/results_from_ag_finished.txt'
#    summ2 = 'E:/RNA_editing_Large_files/MQ/results_from_all_non_ag_finished.txt'
#    summ_files = [summ1,summ2]
#    
#    peps1 = 'E:/RNA_editing_Large_files/MQ/peptides_lists_from_ag_finished.txt'
#    peps2 = 'E:/RNA_editing_Large_files/MQ/peptides_lists_from_all_non_ag_finished.txt'
#    peps_files = [peps1,peps2]
#    
#    names = ['AG','NonAG']
    
    
    fig_path = '/'.join(peps_files[0].split('/')[:-1]) + '/Discovered_Editing_Sites.jpg'
    
    

    
    peps_dfs_dict = {}
    for i,file in enumerate(peps_files):
        peps_dfs_dict.update({names[i]:pd.read_csv(file,sep = '\t', header = 0)})
    
    summaries_dfs_dict = {}
    for i,file in enumerate(summ_files):
        summaries_dfs_dict.update({names[i]:pd.read_csv(file,sep = '\t', header = 0)})
        
    samples_in_all = list(set.intersection(*map(set,[list(summaries_dfs_dict[n]['sample_name']) for n in names])))
    for k,v in summaries_dfs_dict.items():
        summaries_dfs_dict[k] = v[v['sample_name'].isin(samples_in_all)]
    for k,v in peps_dfs_dict.items():
        peps_dfs_dict[k] = v[v['sample_name'].isin(samples_in_all)]
    
    
    print('\nAnalysis for samples:')
    for s in samples_in_all:
        print(s)
    print('\n')
    
    #list of sites (as a non-unique list) discovered in each file
    sites = [get_sites_list(peps_dfs_dict[n]) for n in names]
    tissues = list(set(summaries_dfs_dict[names[0]]['tissue']))
    
    plot_editings_sites_repetitions(fig_path, names, sites)
    
    for i, n in enumerate(names):
        path = '/'.join(summ_files[i].split('/')[:-1]) + '/peptides_for_samples_' + n + '.jpg' 
        c = Counter(sites[i]).items()
        print('\n\n'+n+' Sites (' + str(len(c)) + '):')
        for key, val in Counter(sites[i]).items():
#            print(key)
            print(key+ ': ' + str(val))
        print('\n')
        scatter_discovered_peptides(summaries_dfs_dict[n], n, path)
    
   
    
    
    

