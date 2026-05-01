#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 16:30:29 2023

@author: Sebastian
"""



import os 
import scanpy as sc
# import scFates as scf
import warnings
import pandas as pd
# import palantir
import seaborn as sns
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import matplotlib as mpl

import sys
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout  
mpl.rcParams['text.color'] = 'black' 
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['xtick.color'] = 'black'
mpl.rcParams['ytick.color'] = 'black'


path = "/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/"
os.chdir(path)

def names_changes_list(name_list): 
    new_list = []
    names = pd.read_csv("data/gene_names_for_scrna.csv", sep = ",", index_col = 0)
    for n in name_list: 
        temp1 = names[names["gene_ids"]==n]
        if len(temp1) > 0:
            nn = temp1.iloc[0,1]
            new_list.append(nn)
            print(nn)
        else: 
            new_list.append(n)
    return new_list

########################
### GRN per cluster 
#####################

### AUXIN - RELATED GRN in EP 
tf_evidence = pd.read_csv("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/red_oficial_pares_cruce.csv", sep = "\t")
tf_evidence = tf_evidence[tf_evidence['brady_CHIPSEQconnecTF_DAPSEQ_Gaudinier_TARGETconnectTF_DHalvarez'].str.contains('DAPSEQ')] ## only considering DAPSEQ information 

genes = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx")
EP_genes = genes[genes['group'] == 'EP']

adata = sc.read_h5ad(path + '/adata_filtered_norm_clusters.h5ad')

 
tf_network_merged =  tf_evidence[tf_evidence['TF'].isin(EP_genes.names) & tf_evidence['TARGET'].isin(EP_genes.names)]


tf_network_filter_by_TF = tf_network_merged.groupby(['TF']).size().reset_index(name='count')
tf_network_filter_by_TF = tf_network_filter_by_TF.sort_values(by = 'count', ascending = False)
tf_network_filter_by_TF['TF_name'] = names_changes_list(tf_network_filter_by_TF['TF'])
tf_network_filter_by_TF['percentage'] = 100 * tf_network_filter_by_TF['count'] / len(EP_genes)
# tf_network_filter_by_TF_marker = tf_network_filter_by_TF[tf_network_filter_by_TF['TF'].isin()]
f, ax = plt.subplots(figsize=(2, 3), dpi = 500)
ax.grid(False)
sns.barplot(x="percentage", y="TF_name", data=tf_network_filter_by_TF[0:10], edgecolor = 'black', linewidth = 1.5, 
            palette='inferno')
plt.xlabel('% target genes \nin cluster')
plt.ylabel('Top TFs')
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/TF_' + 'EP' +'_.pdf', bbox_inches='tight')
plt.show()

tf_network_merged['TF_name'] = names_changes_list(tf_network_merged['TF'])
tf_network_merged['TARGET_name'] = names_changes_list(tf_network_merged['TARGET'])
tf_network_merged.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/'  +  'EP_GRN.xlsx')
 










