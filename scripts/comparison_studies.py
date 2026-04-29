#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:05:22 2025

@author: Sebastian
"""


import pandas as pd 
import matplotlib.cm
import matplotlib.pyplot as plt  
import seaborn as sns
import scanpy as sc
import numpy as np
import os 
import matplotlib.colors as mcolors

from pySankey.sankey import sankey
import matplotlib as mpl


mpl.rcParams['font.family'] = ['Arial']
sns.set_style("whitegrid", {'axes.grid' : False})

os.chdir('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/')

def names_changes_list(name_list): 
    new_list = []
    names = pd.read_csv("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/gene_names_for_scrna.csv", sep = ",", index_col = 0)
    for n in name_list: 
        temp1 = names[names["gene_ids"]==n]
        if len(temp1) > 0:
            nn = temp1.iloc[0,1]
            new_list.append(nn)
            print(nn)
        else: 
            new_list.append(n)
    return new_list
   


# %% Import marker genes from available single-cell datasets
sn = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)


# %% Comparing against Guo et al., 2025. Early flower dataset 
guo2025 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Guo(2024)early_flower_marker_genes_TableS1_4l.xlsx')
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

rename_dict = {  # cluster should be in the orther from 0 to final cluster
    '0': ' 0_cortex',
    '1': '1_tapetum',
    '2':  '2_tapetum', 
    '3': '3_tapetum' , 
    '4': '4_guard cell', 
    '5': '5_phloem parenchyma',
    '6':  '6_pollen' , 
    '7': '7_procambium',   
    '8':  '8_cortex',
    '9': '9_epidermis' , 
    '10': '10_companion cells' , 
    '11':  '11_dividing cells',
    '12': '12_microsporocyte',
    '13': '13_epidermis', 
    '14':  '14_phloem parenchyma' , 
    '15': '15_xylem', 
    '16': '16_tapetum',
    '17': '17_tapetum', 
    '18': '18_tapetum', 
    '19': '19_tapetum', 
    '20': '20_tapetum', 
    '21': '21_tapetum', 
    '22' : '22_pollen'
}
guo2025['cluster'] = guo2025['cluster'].astype('string')
guo2025['cluster'] = guo2025['cluster'].map(rename_dict)

## preparing dataset
sankey_moreno_guo = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = guo2025[guo2025['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['guo_cluster'] = domains
        sankey_moreno_guo = pd.concat([sankey_moreno_guo, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['guo_cluster'] = 'unmatched'
        sankey_moreno_guo = pd.concat([sankey_moreno_guo, temp1_sn])


## Ploting sankey 
# fig, ax = plt.subplots(figsize=(3, 3))  # Adjust the size as needed
ax = sankey(
    sankey_moreno_guo['group'], sankey_moreno_guo['guo_cluster'], 
    fontsize=20, palette_colors= 'tab20'
)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/clv_snrnaseq/analysis/sankey.pdf',  bbox_inches='tight')
plt.show()

## without unmatched 
sankey_moreno_guo = sankey_moreno_guo[sankey_moreno_guo['guo_cluster'] != 'unmatched']

## Ploting sankey 
ax = sankey(
    sankey_moreno_guo['group'], sankey_moreno_guo['guo_cluster'], 
    fontsize=20, palette_colors= 'tab20', strip_threshold= 30, 
)
# plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/clv_snrnaseq/analysis/sankey.pdf',  bbox_inches='tight')
plt.show()


# %%  Lee et al., 2025. Flower
lee2025 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Lee(2025)flower_marker_genes_Table4.xlsx')
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)


rename_dict = {  # cluster should be in the orther from 0 to final cluster
    '0': 'Unknown',
    '1': 'Unknown',
    '2':  'Epidermal', 
    '3': 'Tapetum' , 
    '4': 'developing ovule', 
    '5': 'Unknown',
    '6':  'Anther' , 
    '7': 'Anther',   
    '8':  'Tapetum',
    '9': 'Male miocyte' , 
    '10': 'Gynoecium and developing ovule' , 
    '11':  'Pollen',
    '12': 'Vascular',

}
lee2025['cluster'] = lee2025['cluster'].astype('string')
lee2025['cluster'] = lee2025['cluster'].map(rename_dict)



## preparing dataset
sankey_moreno_lee = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = lee2025[lee2025['Gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['lee_cluster'] = domains
        sankey_moreno_lee = pd.concat([sankey_moreno_lee, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['lee_cluster'] = 'unmatched'
        sankey_moreno_lee = pd.concat([sankey_moreno_lee, temp1_sn])

sankey_moreno_lee.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Lee(2025)cluster_per_gene.xlsx')
## without unmatched 
sankey_moreno_lee = sankey_moreno_lee[sankey_moreno_lee['lee_cluster'] != 'unmatched']

## Ploting sankey 
ax = sankey(
    sankey_moreno_lee['group'], sankey_moreno_lee['lee_cluster'], 
    fontsize=20, palette_colors= 'tab20', strip_threshold= 100, 
)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/lee_flower_sankey.pdf',  bbox_inches='tight')
plt.show()

 # %%  Neumann (2022) 
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

neumann2022_names = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/neumann(2022)cluster_names.xlsx')
neumann2022_markers = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/neumann(2022)cluster_markers.xlsx') 

name_map = dict(zip(sorted(neumann2022_markers['cluster'].unique()), neumann2022_names.annotation.tolist()))
neumann2022_markers['cluster'] = neumann2022_markers['cluster'].replace(name_map)

## preparing dataset
sankey_moreno_neumann = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = neumann2022_markers[neumann2022_markers['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['neumann_cluster'] = domains
        sankey_moreno_neumann = pd.concat([sankey_moreno_neumann, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['neumann_cluster'] = 'unmatched'
        sankey_moreno_neumann = pd.concat([sankey_moreno_neumann, temp1_sn])

sankey_moreno_neumann.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Neumann(2022)_cluster_per_gene.xlsx')

## without unmatched 
sankey_moreno_neumann = sankey_moreno_neumann[sankey_moreno_neumann['neumann_cluster'] != 'unmatched']

## Ploting sankey 
ax = sankey(
    sankey_moreno_neumann['group'], sankey_moreno_neumann['neumann_cluster'], 
    fontsize=20, palette_colors= 'tab20c', strip_threshold= 50, 
)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/neumann_sankey.pdf',  bbox_inches='tight')
plt.show()


# %% Yadav  2014  - dataset 

yadav2014 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/yadav(2014)_DEG.xlsx')
cols = [
    'CLV3_LAS_KAN1',
    'HMG_HDG4_WUS',
    'CLV3_FIL_KAN1_LAS',
    'HMG_HDG4_AtML1',
    'AtHB8_S17Soot_HDG4'
]

yadav2014_filtered = yadav2014[   ## ONLY DEG IN YADAV
    yadav2014[cols].notna().any(axis=1)
]
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)
# markers = markers[markers['group'].isin(['undifferentiated cells'])]



yadav_markers = pd.DataFrame()  ## obtain where the genes are differentially expressed in yadav (2009, 2014) datasets 
for gene in yadav2014_filtered.AGI.unique(): 
    temp1 = yadav2014_filtered[yadav2014_filtered['AGI'] == gene]
    temp1 = temp1[['AtHB8', 'CLV3', 'FIL', 'HDG4', 'HMG',
           'KAN1', 'AtML1', 'LAS', 'S17Shoot', 'WUS']]
    sorted_cols = temp1.iloc[0].sort_values(ascending=False)
    temp2 = sorted_cols.index[0]
    d1 = {'gene': gene, 'domain': temp2}
    df1 = pd.DataFrame(data=d1, index = [0])
    yadav_markers = pd.concat([yadav_markers, df1])
        
yadav_markers['name'] =  names_changes_list(yadav_markers.gene)


## preparing dataset
sankey_y_sm = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = yadav_markers[yadav_markers['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['domain'].tolist()[0]
        temp1_sn['yadav'] = domains
        sankey_y_sm = pd.concat([sankey_y_sm, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['yadav'] = 'unmatched'
        sankey_y_sm = pd.concat([sankey_y_sm, temp1_sn])


sankey_y_sm.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/yadav(2014)_cluster_per_gene.xlsx')

# sankey_y_sm_fil = sankey_y_sm[sankey_y_sm['yadav'] != 'unmatched']
sankey_y_sm_fil = sankey_y_sm
# sankey_y_sm_fil = sankey_y_sm_fil[sankey_y_sm_fil['yadav'].isin(['WUS', 'CLV3', 'HDG4', 'unmatched'])]
## Ploting sankey 
ax = sankey(
    sankey_y_sm_fil['group'], sankey_y_sm_fil['yadav'], palette_colors = 'tab20', 
    fontsize=20, strip_threshold= 25, 
)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/yadav_sarkey.pdf',  bbox_inches='tight')
plt.show()

# %%  Tiang (2019) comparison Vegetative apex 

tiang2019 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Tian(2019)TableII.xlsx')
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)
# markers = markers[markers['group'].isin(['undifferentiated cells', 'procambium', 'IM epidermis', 'stem epidermis'])]



## preparing dataset
sankey_tiang_sm = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = tiang2019[tiang2019['GID'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['Specifically expressed domain'].tolist()[0]
        temp1_sn['tiang'] = domains
        sankey_tiang_sm = pd.concat([sankey_tiang_sm, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['tiang'] = 'unmatched'
        sankey_tiang_sm = pd.concat([sankey_tiang_sm, temp1_sn])


sankey_tiang_sm.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Tiang(2019)cluster_per_gene.xlsx')


sankey_tiang_sm = sankey_tiang_sm[sankey_tiang_sm['tiang'] != 'unmatched']


## Ploting sankey 
# fig, ax = plt.subplots(dpi = 500)  # Adjust the size as needed
ax = sankey(
    sankey_tiang_sm['group'], sankey_tiang_sm['tiang'], palette_colors = 'tab20', 
    fontsize=20, strip_threshold= 10,  figure_name  = 'Tian'
)
# plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/tian(2019).pdf',  bbox_inches='tight')
plt.show()

# 

# %%  ### Zhang (2021) comparison Vegetative apex (STM+ cells)
zhang2021 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/zhang(2021)STM+.xlsx')
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

rename_dict = {  # cluster should be in the orther from 0 to final cluster
    '0': 'undifferentiated cells',
    '1': 'transport process',
    '2':  'phloem', 
    '3': 'response to stress' , 
    '4': 'BD', 
    '5': 'Xylem',
    '6':  'Epidermis' , 
    '7': 'Vasculature',   
    '8':  'Abaxial domain',

}
zhang2021['cluster'] = zhang2021['cluster'].astype('string')
zhang2021['cluster'] = zhang2021['cluster'].map(rename_dict)

sankey_z_sm = pd.DataFrame()
## Combining marker genes 
for g in sn.names.unique(): 
    temp1_sn = sn[sn['names'] == g]
    temp1_query = zhang2021[zhang2021['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['zhang'] = domains
        sankey_z_sm = pd.concat([sankey_z_sm, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['zhang'] = 'unmatched'
        sankey_z_sm = pd.concat([sankey_z_sm, temp1_sn])
        
        
## preparing dataset
sankey_moreno_zhang = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = zhang2021[zhang2021['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['zhang2021_cluster'] = domains
        sankey_moreno_zhang = pd.concat([sankey_moreno_zhang, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['zhang2021_cluster'] = 'unmatched'
        sankey_moreno_zhang = pd.concat([sankey_moreno_zhang, temp1_sn])

sankey_moreno_zhang.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Zhang(2021)cluster_per_gene.xlsx')

sankey_moreno_zhang = sankey_moreno_zhang[sankey_moreno_zhang['zhang2021_cluster'] != 'unmatched']


## Ploting sankey 
ax = sankey(
    sankey_moreno_zhang['group'], sankey_moreno_zhang['zhang2021_cluster'], palette_colors = 'tab20c', 
    fontsize=20, strip_threshold= 50, 
)
# plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/zhang_sarkey.pdf',  bbox_inches='tight')
plt.show()


# %%  ### Zhang (2021) comparison Vegetative apex 
zhang2021 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/zhang(2021)_all_clusters.xlsx')
markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

rename_dict = {  # cluster should be in the orther from 0 to final cluster
    '0': 'mesophyll cell (MC)_1',
    '1': 'shoot meristematic cells (SMC)_1',
    '2':  'epidermal cell (EC)_1', 
    '3': 'mesophyll cell (MC)_2' , 
    '4': 'mesophyll cell (MC)_3', 
    '5': 'proliferating cell (PC) 1',
    '6':  'vascular cell (VC) 1' , 
    '7': 'vascular cell (VC) 2',   
    '8':  'undefined cell 1',
    '9': 'proliferating cell (PC) 2',
    '10': 'vascular cell (VC) 3',
    '11':  'guard cell (GC)', 
    '12': 'epidermal cell (EC)_2' , 
    '13': 'shoot meristematic cells (SMC)_2', 
    '14': 'epidermal cell (EC)_3',
    '15':  'companion cell (CC)' , 
    '16': 'mesophyll cell (MC)_4',   
    '17':  'proliferating cell (PC) 3',
    '18': 'epidermal cell (EC)_4', 
    '19': 'proliferating cell (PC) 4',
    '20':  'shoot endodermis (SEn)' , 
    '21': 'undefined cell 1',   
    '22':  'vascular cell (VC) 4',

}
zhang2021['cluster'] = zhang2021['cluster'].astype('string')
zhang2021['cluster'] = zhang2021['cluster'].map(rename_dict)

sankey_z_sm = pd.DataFrame()
## Combining marker genes 
for g in sn.names.unique(): 
    temp1_sn = sn[sn['names'] == g]
    temp1_query = zhang2021[zhang2021['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['zhang'] = domains
        sankey_z_sm = pd.concat([sankey_z_sm, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['zhang'] = 'unmatched'
        sankey_z_sm = pd.concat([sankey_z_sm, temp1_sn])
        
        
## preparing dataset
sankey_moreno_zhang = pd.DataFrame()
## Combining marker genes 
for g in markers.names.unique(): 
    temp1_sn = markers[markers['names'] == g]
    temp1_query = zhang2021[zhang2021['gene'] == g]
    if len(temp1_query) > 0: 
        domains = temp1_query['cluster'].tolist()[0]
        temp1_sn['zhang2021_cluster'] = domains
        sankey_moreno_zhang = pd.concat([sankey_moreno_zhang, temp1_sn])
    if len(temp1_query) == 0: 
        temp1_sn['zhang2021_cluster'] = 'unmatched'
        sankey_moreno_zhang = pd.concat([sankey_moreno_zhang, temp1_sn])

sankey_moreno_zhang.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Zhang(2021)cluster_per_gene_all.xlsx')

sankey_moreno_zhang = sankey_moreno_zhang[sankey_moreno_zhang['zhang2021_cluster'] != 'unmatched']

## Ploting sankey 
# fig, ax = plt.subplots(dpi = 300)  # Adjust the size as needed
ax = sankey(
    sankey_moreno_zhang['group'], sankey_moreno_zhang['zhang2021_cluster'], palette_colors = 'tab20c', 
    fontsize=20, strip_threshold= 50, 
)
# plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/zhang_sarkey.pdf',  bbox_inches='tight')
plt.show()
# %% Generate multidimensional dataset comparison
## Previous DEG genes
zhang2021_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Zhang(2021)cluster_per_gene.xlsx', index_col= 0)
neumann2022_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Neumann(2022)_cluster_per_gene.xlsx', index_col= 0)
yadav2014_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/yadav(2014)_cluster_per_gene.xlsx', index_col= 0)
tiang2019_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Tiang(2019)cluster_per_gene.xlsx', index_col= 0)
# lee2025_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Lee(2025)cluster_per_gene.xlsx', index_col= 0)

## DEG genes from this study
moreno_DEG_cluster = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

# moreno_DEG_cluster['lee2025_flower'] = lee2025_DEG_cluster['lee_cluster']
moreno_DEG_cluster['neumann2022'] = neumann2022_DEG_cluster['neumann_cluster']
moreno_DEG_cluster['zhang2021'] = zhang2021_DEG_cluster['zhang2021_cluster']
moreno_DEG_cluster['tiang2019'] = tiang2019_DEG_cluster['tiang']
moreno_DEG_cluster['yadav2014'] = yadav2014_DEG_cluster['yadav']


moreno_DEG_cluster.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final_comparison.xlsx")


# %% Generate multidimensional dataset comparison
## Previous DEG genes
zhang2021_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Zhang(2021)cluster_per_gene.xlsx', index_col= 0)
neumann2022_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Neumann(2022)_cluster_per_gene.xlsx', index_col= 0)
yadav2014_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/yadav(2014)_cluster_per_gene.xlsx', index_col= 0)
tiang2019_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Tiang(2019)cluster_per_gene.xlsx', index_col= 0)
lee2025_DEG_cluster = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Lee(2025)cluster_per_gene.xlsx', index_col= 0)

## DEG genes from this study
moreno_DEG_cluster = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

moreno_DEG_cluster['lee2025_flower'] = lee2025_DEG_cluster['lee_cluster']
moreno_DEG_cluster['neumann2022'] = neumann2022_DEG_cluster['neumann_cluster']
moreno_DEG_cluster['zhang2021'] = zhang2021_DEG_cluster['zhang2021_cluster']
moreno_DEG_cluster['tiang2019'] = tiang2019_DEG_cluster['tiang']
moreno_DEG_cluster['yadav2014'] = yadav2014_DEG_cluster['yadav']


moreno_DEG_cluster.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final_comparison.xlsx")


# %% Multidimensional dataset to obsere overlap between samples (Supplementary Figure 5)
###########

comparison = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final_comparison.xlsx')
comparison = comparison[['group', 'names',  'yadav2014','neumann2022', 'zhang2021' ]]
comparison.yadav2014.unique()
comparison.neumann2022.unique()
comparison.zhang2021.unique()
comparison.group.unique()

comparison['yadav_new'] = comparison['yadav2014'].replace({ 'AtML1': 'epidermis', 
                                                           'AtHB8': 'vasculature',
                                                           'S17Shoot' : 'vasculature',
                                                           'WUS' : 'undifferentiated cells',
                                                           'CLV3' : 'undifferentiated cells',
                                                            'HDG4' : 'undifferentiated cells',
                                                           'HMG' : 'epidermis', 
                                                           'LAS' : 'BD' , 
                                                           'FIL' : 'EP'})

comparison['neumann_new'] = comparison['neumann2022'].replace({ 'xylem parenchyma': 'vasculature',
                                                           '(pro-)cambium' : 'vasculature',
                                                           'phloem' : 'vasculature',
                                                           'epidermis/dividing' : 'epidermis',
                                                           'floral meristem' : 'floral identity',
                                                           'S-phase cells': 'cell cycle',
                                                           'dividing cells': 'cell cycle', 
                                                           'mesophyll' : 'photosynthetic cells',
                                                           })

comparison['zhang_new'] = comparison['zhang2021'].replace({ 'response to stress': 'stress-response',
                                                       'Epidermis' : 'epidermis',
                                                       'Xylem' : 'vasculature',
                                                       'phloem' : 'vasculature',
                                                       'Vasculature' : 'vasculature'
                                                          })

comparison['moreno_new'] = comparison['group'].replace({ 'stem epidermis': 'epidermis',
                                                        'xylem parenchyma' : 'vasculature' ,
                                                           'procambium' : 'vasculature' ,
                                                               'IM epidermis' : 'epidermis',    
                                                               'S-phase': 'cell cycle',
                                                               'G2/M-phase' : 'cell cycle'
                                                          })


comparison.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_comparisons_miltidim.xlsx')


comparison = comparison.replace("unmatched", np.nan)

def consensus_with_datasets(row):
    row = row.dropna()
    if row.empty:
        return pd.Series([0, []])
    
    # Count the frequency of each label
    counts = row.value_counts()
    top_label = counts.idxmax()
    
    # Datasets that agree on top label
    agreeing_datasets = row[row == top_label].index.tolist()
    
    # Score = how many agree / total annotated
    score = len(agreeing_datasets) / len(row)
    
    return pd.Series([score, agreeing_datasets])

# Apply to your dataframe
comparison[['agreement_score', 'agreeing_datasets']] = comparison[['yadav_new', 'neumann_new', 'zhang_new', 'moreno_new']].apply(consensus_with_datasets, axis=1)

### Yadav (2014) comparison
yadav2014 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/yadav(2014)_DEG.xlsx')
yadav2014_fil =  yadav2014[['AGI', 'CLV3_LAS_KAN1',
                        'HMG_HDG4_WUS', 'CLV3_FIL_KAN1_LAS', 'HMG_HDG4_AtML1',
                        'AtHB8_S17Soot_HDG4']]
yadav2014_fil2 = yadav2014_fil.dropna(subset=yadav2014_fil.columns[1:], how='all') ## Marker genes in Yadav


neumann2022_names = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/neumann(2022)cluster_names.xlsx')
neumann2022_markers = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/neumann(2022)cluster_markers.xlsx') 


zhang2021 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/zhang(2021)STM+.xlsx')
    

marker_counts = {
    'yadav_new': yadav2014_fil2.shape[0],
    'neumann_new': neumann2022_markers.shape[0],
    'zhang_new': zhang2021.shape[0],
    'moreno_new': comparison.shape[0]   # <-- whatever your object is
}

all_marker_genes = list(set(yadav2014_fil2.AGI.tolist() +   neumann2022_markers.gene.tolist()  +  zhang2021.gene.tolist() +  comparison.names.tolist()))


datasets = ['yadav_new', 'neumann_new', 'zhang_new', 'moreno_new']

# Initialize square matrix
pair_matrix = pd.DataFrame(0, index=datasets, columns=datasets)

# Iterate through the agreeing_datasets column
for lst in comparison['agreeing_datasets']:
    # Only consider lists with more than one dataset (i.e., actual agreement)
    if len(lst) > 1:
        # Generate all unique pairs (combinations) from the list
        for i in range(len(lst)):
            for j in range(i+1, len(lst)):
                d1 = lst[i]
                d2 = lst[j]
                # Increment only the off-diagonal cells
                pair_matrix.loc[d1, d2] += 1
                pair_matrix.loc[d2, d1] += 1

for ds in datasets:
    pair_matrix.loc[ds, ds] = marker_counts[ds]


pair_matrix = pair_matrix.rename(index={
    'yadav_new': 'FACS-sorted',
    'neumann_new': 'floral meristem',
    'zhang_new': 'STM+ veg.apex',
    'moreno_new': 'This study'
})

pair_matrix = pair_matrix.rename(columns={
    'yadav_new': 'FACS-sorted',
    'neumann_new': 'floral meristem',
    'zhang_new': 'STM+ veg.apex',
    'moreno_new': 'This study'
})


# Plot heatmap number of shared genes in similar cluster 
plt.figure(figsize=(6,5), dpi = 300)
ax = sns.heatmap(pair_matrix, annot=True, fmt='d',
            vmax=1500, 
            annot_kws={'size': 15},   
            cmap='viridis', 
            cbar_kws={'label': "Shared genes in similar clusters"})
# plt.title("Shared genes across datasets in similar clusters")
ax.tick_params(axis='both', labelsize=17)

plt.subplots_adjust(right=0.85)  # increase space for colorbar
plt.show()



# Dataset sizes (diagonal)
dataset_sizes = pair_matrix.values.diagonal()
sizes = dict(zip(pair_matrix.index, dataset_sizes))


percent_relative = pd.DataFrame(0.0, index=pair_matrix.index.tolist(), columns=pair_matrix.columns.tolist())


for d1 in pair_matrix.index.tolist():
    for d2 in pair_matrix.index.tolist():
        inter = pair_matrix.loc[d1, d2]
        percent_relative.loc[d1, d2] = inter / sizes[d2] * 100  # relative to B

print(percent_relative)

plt.figure(figsize=(6,5), dpi=300)
ax = sns.heatmap(
    percent_relative, 
    annot=True, 
    fmt=".2f", 
    cmap="viridis", 
    vmax=70, 
    annot_kws={'size': 15},   
    cbar_kws={'label': "% shared marker genes in similar cluster"}
)

ax.tick_params(axis='both', labelsize=17)
# Increase colorbar title size
cbar = ax.collections[0].colorbar

cbar.set_label("% shared marker genes \nin similar cluster", size=20)

plt.subplots_adjust(right=0.85)  # increase space for colorbar
plt.show()


# %%  Lee (2025) stem 
# sn = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)
# lee2025_stem = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Lee(2025)stem_marker_genes_Table4.xlsx')


# rename_dict = {  # cluster should be in the orther from 0 to final cluster
#     '0': 'Unknown',
#     '1': 'Unknown',
#     '2':  'Unknown', 
#     '3': 'Xylem' , 
#     '4': 'Vascular', 
#     '5': 'Guard',
#     '6':  'Unknown' , 
#     '7': 'Unknown',   
#     '8':  'Xylem',
#     '9': 'Phloem' , 
#     '10': 'Epidermal' , 
#     '11':  'Cambium',
#     '12': 'Epidermal',
#     '13' : 'Xylem',
#     '14' : 'Phloem',
#     '15' : 'Trichome'

# }

# lee2025_stem['cluster'] = lee2025_stem['cluster'].astype('string')
# lee2025_stem['cluster'] = lee2025_stem['cluster'].map(rename_dict)

# ## preparing dataset
# sankey_lee_stem_sm = pd.DataFrame()
# ## Combining marker genes 
# for g in sn.names.unique(): 
#     temp1_sn = sn[sn['names'] == g]
#     temp1_query = lee2025_stem[lee2025_stem['Gene'] == g]
#     if len(temp1_query) > 0: 
#         domains = temp1_query['cluster'].tolist()[0]
#         temp1_sn['lee'] = domains
#         sankey_lee_stem_sm = pd.concat([sankey_lee_stem_sm, temp1_sn])
#     if len(temp1_query) == 0: 
#         temp1_sn['lee'] = 'unmatched'
#         sankey_lee_stem_sm = pd.concat([sankey_lee_stem_sm, temp1_sn])


# sankey_lee_stem_sm = sankey_lee_stem_sm[sankey_lee_stem_sm['lee'] != 'unmatched']


# ## Ploting sankey 
# fig, ax = plt.subplots(figsize=(3, 3))  # Adjust the size as needed
# ax = sankey(
#     sankey_lee_stem_sm['group'], sankey_lee_stem_sm['lee'],
#     fontsize=20, strip_threshold= 10, 
# )
# plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/lee_stem_sankey.pdf',  bbox_inches='tight')
# plt.show()

# # %%  Guo (2025) stem 
# sn = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)
# guo2025_stem = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Guo(2025)stem_marker_genes_TableS1.xlsx')


# rename_dict = {  # cluster should be in the orther from 0 to final cluster
#     '0': 'pith',
#     '1': 'cortex',
#     '2':  'epidermis', 
#     '3': 'cortex' , 
#     '4': 'procambium', 
#     '5': 'cortex',
#     '6':  'phloem parenchyma' , 
#     '7': 'xylem/interfascicular fiber',   
#     '8':  'companion cell',
#     '9': 'unknown' , 
#     '10': 'starch sheath' , 
#     '11':  'vascular cell',
#     '12': 'cortex',
#     '13' : 'guard cell',
#     '14' : 'xylem/interfascicular fiber',

# }
# guo2025_stem['cluster'] = guo2025_stem['cluster'].astype('string')
# guo2025_stem['cluster'] = guo2025_stem['cluster'].map(rename_dict)

# ## preparing dataset
# sankey_guo_stem_sm = pd.DataFrame()
# ## Combining marker genes 
# for g in sn.names.unique(): 
#     temp1_sn = sn[sn['names'] == g]
#     temp1_query = guo2025_stem[guo2025_stem['gene'] == g]
#     if len(temp1_query) > 0: 
#         domains = temp1_query['cluster'].tolist()[0]
#         temp1_sn['guo'] = domains
#         sankey_guo_stem_sm = pd.concat([sankey_guo_stem_sm, temp1_sn])
#     if len(temp1_query) == 0: 
#         temp1_sn['guo'] = 'unmatched'
#         sankey_guo_stem_sm = pd.concat([sankey_guo_stem_sm, temp1_sn])


# sankey_guo_stem_sm = sankey_guo_stem_sm[sankey_guo_stem_sm['guo'] != 'unmatched']


# ## Ploting sankey 
# fig, ax = plt.subplots(figsize=(3, 3))  # Adjust the size as needed
# ax = sankey(
#     sankey_guo_stem_sm['group'], sankey_guo_stem_sm['guo'],
#     fontsize=12, strip_threshold= 10, 
# )
# plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/lee_stem_sankey.pdf',  bbox_inches='tight')
# plt.show()


# # %% Comparing against Guo et al., 2025. Middle flower dataset 
# guo2025 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Guo(2024)middle_flower_marker_genes_TableS1_4l.xlsx')
# markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

# rename_dict = {  # cluster should be in the orther from 0 to final cluster
#     '0': ' 0_cortex',
#     '1': '1_cortex',
#     '2':  '2_pollen', 
#     '3': '3_sperm' , 
#     '4': '4_pollen (vegetative nucleus)', 
#     '5': '5_pollen',
#     '6':  '6_guard cell' , 
#     '7': '7_dividing cell',   
#     '8':  '8_companion cell',
#     '9': '9_epidermis' , 
#     '10': '10_dividing cell' , 
#     '11':  '11_phloem parenchyma',
#     '12': '12_phloem parenchyma',
#     '13': '13_unknown', 
#     '14':  '14_stigma' , 
#     '15': '15_pollen', 
#     '16': '16_nectary',
#     '17': '17_procambium-xylem', 

# }
# guo2025['cluster'] = guo2025['cluster'].astype('string')
# guo2025['cluster'] = guo2025['cluster'].map(rename_dict)

# ## preparing dataset
# sankey_moreno_guo = pd.DataFrame()
# ## Combining marker genes 
# for g in markers.names.unique(): 
#     temp1_sn = markers[markers['names'] == g]
#     temp1_query = guo2025[guo2025['gene'] == g]
#     if len(temp1_query) > 0: 
#         domains = temp1_query['cluster'].tolist()[0]
#         temp1_sn['guo_cluster'] = domains
#         sankey_moreno_guo = pd.concat([sankey_moreno_guo, temp1_sn])
#     if len(temp1_query) == 0: 
#         temp1_sn['guo_cluster'] = 'unmatched'
#         sankey_moreno_guo = pd.concat([sankey_moreno_guo, temp1_sn])


# # ## Ploting sankey 
# # # fig, ax = plt.subplots(figsize=(3, 3))  # Adjust the size as needed
# # ax = sankey(
# #     sankey_moreno_guo['group'], sankey_moreno_guo['guo_cluster'], 
# #     fontsize=20, palette_colors= 'tab20'
# # )
# # plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/clv_snrnaseq/analysis/sankey.pdf',  bbox_inches='tight')
# # plt.show()

# ## without unmatched 
# sankey_moreno_guo = sankey_moreno_guo[sankey_moreno_guo['guo_cluster'] != 'unmatched']

# ## Ploting sankey 
# ax = sankey(
#     sankey_moreno_guo['group'], sankey_moreno_guo['guo_cluster'], 
#     fontsize=20, palette_colors= 'tab20', strip_threshold= 50, 
# )
# # plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/clv_snrnaseq/analysis/sankey.pdf',  bbox_inches='tight')
# plt.show()

# # %% Comparing against Guo et al., 2025. Late flower dataset 
# guo2025 = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/Guo(2024)late_flower_marker_genes_TableS1_4l.xlsx')
# markers = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx", index_col= 0)

# rename_dict = {  # cluster should be in the orther from 0 to final cluster
#     '0': ' 0_cortex',
#     '1': '1_pollen',
#     '2':  '2_pollen', 
#     '3': '3_sperm' , 
#     '4': '4_integument', 
#     '5': '5_pollen(vegetative nucleus)',
#     '6':  '6_cortex' , 
#     '7': '7_companion cell',   
#     '8':  '8_guard cell',
#     '9': '9_dividing cell' , 
#     '10': '10_epidermis' , 
#     '11':  '11_pollen (vegetative nucleus',
#     '12': '12_procambium-xylem',
#     '13': '13_cortex', 
#     '14':  '14_dividing cell' , 
#     '15': '15_abscission zone', 

# }
# guo2025['cluster'] = guo2025['cluster'].astype('string')
# guo2025['cluster'] = guo2025['cluster'].map(rename_dict)

# ## preparing dataset
# sankey_moreno_guo = pd.DataFrame()
# ## Combining marker genes 
# for g in markers.names.unique(): 
#     temp1_sn = markers[markers['names'] == g]
#     temp1_query = guo2025[guo2025['gene'] == g]
#     if len(temp1_query) > 0: 
#         domains = temp1_query['cluster'].tolist()[0]
#         temp1_sn['guo_cluster'] = domains
#         sankey_moreno_guo = pd.concat([sankey_moreno_guo, temp1_sn])
#     if len(temp1_query) == 0: 
#         temp1_sn['guo_cluster'] = 'unmatched'
#         sankey_moreno_guo = pd.concat([sankey_moreno_guo, temp1_sn])


# # ## Ploting sankey 
# # # fig, ax = plt.subplots(figsize=(3, 3))  # Adjust the size as needed
# # ax = sankey(
# #     sankey_moreno_guo['group'], sankey_moreno_guo['guo_cluster'], 
# #     fontsize=20, palette_colors= 'tab20'
# # )
# # plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/clv_snrnaseq/analysis/sankey.pdf',  bbox_inches='tight')
# # plt.show()

# ## without unmatched 
# sankey_moreno_guo = sankey_moreno_guo[sankey_moreno_guo['guo_cluster'] != 'unmatched']

# ## Ploting sankey 
# ax = sankey(
#     sankey_moreno_guo['group'], sankey_moreno_guo['guo_cluster'], 
#     fontsize=20, palette_colors= 'Accent', strip_threshold= 50, 
# )
# # plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/clv_snrnaseq/analysis/sankey.pdf',  bbox_inches='tight')
# plt.show()
