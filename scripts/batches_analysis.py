#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 12:29:18 2025

@author: Sebastian
"""

import pandas as pd 
import matplotlib.cm
import matplotlib.pyplot as plt  
import seaborn as sns
import scanpy as sc
import numpy as np
import os 

import matplotlib as mpl
mpl.rcParams['font.family'] = ['Arial']


# %% UMAP per sample(batch) 

adata =  sc.read_h5ad(OUTPUT_DATA / 'adata_filtered_norm_clusters.h5ad') ## normalized counts single-nuclei


### Metrics per batch 

count = adata.to_df()
# count['sample'] = adata.obs['Sample']
detected_genes_per_row = pd.DataFrame((count > 0).sum(axis=1))
detected_genes_per_row['sample'] = adata.obs['Sample']


### Separate by batch 

sc.set_figure_params(figsize= (4,4), frameon=True, dpi=300)
batch1 = adata[adata.obs['Sample'] == 'sittc6']
sc.pl.umap(batch1, color=['leiden' ], title = 'batch1', 
           add_outline=False,alpha =0.8, size = 30,
       frameon=True, legend_fontsize=7, legend_fontoutline=2,
       palette='tab20')

batch2 = adata[adata.obs['Sample'] == 'sitta1']
sc.pl.umap(batch2, color=['leiden' ], title= 'batch2',
           add_outline=False,alpha =0.8, size = 40,
       frameon=True, legend_fontsize=7, legend_fontoutline=2,
       palette='tab20')

batch3 = adata[adata.obs['Sample'] == 'sitte2']
sc.pl.umap(batch3, color=['leiden' ], title = 'batch3', 
           add_outline=False,alpha =0.8, size = 40,
       frameon=True, legend_fontsize=7, legend_fontoutline=2,
       palette='tab20')

len(adata.obs[adata.obs['Sample'] == 'sittc6'])
len(adata.obs[adata.obs['Sample'] == 'sitta1'])
len(adata.obs[adata.obs['Sample'] == 'sitte2'])

# %%  Correlation between replicates 

adata = sc.read_h5ad(OUTPUT_DATA / 'adata_filtered_norm_clusters.h5ad') ## normalized counts single-nuclei

adata1 = adata[adata.obs['Sample'] == 'sittc6']
adata2 = adata[adata.obs['Sample'] == 'sitta1']
adata3 = adata[adata.obs['Sample'] == 'sitte2']

count1 = adata1.to_df().sum(axis = 0)
count2 = adata2.to_df().sum(axis = 0)
count3 = adata3.to_df().sum(axis = 0)

count1= np.log2(count1 +1 )
count2= np.log2(count2 +1 )
count3= np.log2(count3 +1 )

fig, ax = plt.subplots(figsize=(3, 3), dpi = 500)
ax.grid(False)
sns.scatterplot(
    x = count1, 
    y= count2, 
    alpha = 0.1,
    size = 0.2, 
    color = 'purple'
    # color=".6", line_kws=dict(color="purple"),  ci = 99,
    # fit_reg = True
)
ax.get_legend().set_visible(False)
plt.xlabel('Rep. 1')
plt.ylabel('Rep. 2')
plt.show()

fig, ax = plt.subplots(figsize=(3, 3), dpi = 500)
ax.grid(False)
sns.scatterplot(
    x = count1, 
    y= count3, 
    alpha = 0.1,
    size = 0.2, 
    color = 'purple'
    # color=".6", line_kws=dict(color="purple"),  ci = 99,
    # fit_reg = True
)
ax.get_legend().set_visible(False)
plt.xlabel('Rep. 1')
plt.ylabel('Rep. 3')
plt.show()

fig, ax = plt.subplots(figsize=(3, 3), dpi = 500)
ax.grid(False)
sns.scatterplot(
    x = count2, 
    y= count3, 
    alpha = 0.1,
    size = 0.2 ,
    color = 'purple'
    # color=".6", line_kws=dict(color="purple"),  ci = 99,
    # fit_reg = True
)
ax.get_legend().set_visible(False)
plt.xlabel('Rep. 2')
plt.ylabel('Rep. 3')
plt.show()

from scipy import stats
print(stats.spearmanr(count1,count2))
print(stats.spearmanr(count1,count3))
print(stats.spearmanr(count2,count3))

# %% TOTAL COUNT AND NUMBER OF GENES PER CLUSTER
fig, ax = plt.subplots(figsize=(6, 4), dpi = 500)
obs = pd.DataFrame(adata.obs)
ax = sns.violinplot(data = obs, 
              x='leiden',
              y='total_counts',  inner = None, 
              linecolor = 'black',   palette= batch1.uns['leiden_colors'],  linewidth= 1,  width=1, 
               )
ax.grid(False)
# Calculate means
group_means = obs.groupby('leiden')['total_counts'].mean()
ax.tick_params(axis='x', rotation=90)
plt.ylabel('Total counts', fontsize =15)
plt.ylim(-200,15000) 
plt.savefig(OUTPUT_FIGURES / 'total_counts.pdf', bbox_inches='tight')
plt.show()

# %%  VIOLIN PLOT PER BATCH 

sc.set_figure_params(frameon=False, dpi=500)
plt.rcParams['axes.grid'] = False
ax = sc.pl.violin(adata, keys='n_genes_by_counts', groupby='Sample', stripplot=False, jitter=False, 
             order = ('sittc6', 'sitta1', 'sitte2'),palette='magma', ylabel='Genes per cell' )

sc.set_figure_params(frameon=False, dpi=500)
plt.rcParams['axes.grid'] = False
ax = sc.pl.violin(adata, keys='total_counts', groupby='Sample', stripplot=False, jitter=False, 
             order = ('sittc6', 'sitta1', 'sitte2'),palette='magma', ylabel='Read count' )

# %%  Ploting % cell per cluster
temp1 = pd.DataFrame(adata.obs['leiden'].value_counts())
temp1['cluster'] = temp1.index
temp1.columns = ['count', 'cluster']
temp1['percentage'] = temp1['count']/len(adata)*100
temp1['cluster'] = temp1.index
plt.figure(figsize=(3, 3))
sns.set(style="ticks")
ax = sns.barplot(x='cluster', y="count",
                 data=temp1, linewidth=1, edgecolor='black',
                 palette=batch1.uns['leiden_colors'])
# ax.bar_label(ax.containers[0]) # only 1 container needed unless using `hue`
plt.setp(ax.artists, edgecolor='k', facecolor='w')
ax.tick_params(axis='x', rotation=90)
ax.set(xlabel='Cell type', ylabel='Number of cells')
plt.setp(ax.lines, color='k')
sns.despine(offset=3, trim=False)
plt.savefig(OUTPUT_FIGURES / 'percentage_clusters.pdf',  bbox_inches='tight')
plt.show()

