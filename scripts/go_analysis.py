#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 11:31:18 2023

@author: Sebastian
"""
import glob
import pandas as pd 
import os 
import matplotlib.pyplot as plt  
import seaborn as sns
import numpy as np
from scipy import stats
import matplotlib as mpl
import anndata as ad
from glob import glob

# %% Folder paths
REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DATA = REPO_ROOT / "input_data"
OUTPUT_FIGURES = REPO_ROOT / "figures"
OUTPUT_DATA = REPO_ROOT / "output_data"

# %% Functions
def names_changes_list(name_list): 
    new_list = []
    names = pd.read_csv(INPUT_DATA / "gene_names_for_scrna.csv", sep = ",", index_col = 0)
    for n in name_list: 
        temp1 = names[names["gene_ids"]==n]
        if len(temp1) > 0:
            nn = temp1.iloc[0,1]
            new_list.append(nn)
            print(nn)
        else: 
            new_list.append(n)
    return new_list

#%% Gene ontology for all clusters 
import pandas as pd
from gprofiler import GProfiler
import numpy as np
trajectory = 'EP'
cc_var = pd.read_excel(
    OUTPUT_DATA / 'marker_genes_reg_final.xlsx',
    index_col=0)

# Instantiate the object
gp = GProfiler(return_dataframe=True)
# window_size = int(len(cc_var) * 0.25)
# num_windows = (len(cc_var) + window_size - 1) // window_size

all_results = []  # To store results of all windows
for cluster in cc_var.group.unique().tolist():
    # window_df = cc_var.iloc[start:end]
    temp1 = cc_var[cc_var['group'] == cluster]
    
    # Extract gene IDs (TAIR IDs)
    gene_list = temp1['names'].tolist()[0:100]
    
    # Run GO enrichment for Arabidopsis
    results = gp.profile(
        organism='athaliana',
        query=gene_list,
        significance_threshold_method='fdr', 
        sources=['GO:BP'], 
        user_threshold=0.05, 
        no_evidences=False
    )
    results['cluster'] = cluster
    results['fold_enrichment'] = ( ## important metric
                                (results['intersection_size'] / results['query_size']) /
                                (results['term_size'] / results['effective_domain_size'])
                                    )
    all_results.append(results)

final_results = pd.concat(all_results, ignore_index=True)
final_results['-log10_pval'] = -np.log10(final_results['p_value'])
## Incresing threhold? 
final_results = final_results[(final_results['p_value'] < 0.01) & (final_results['intersection_size'] > 5)
                              & (final_results['term_size'] <1000)]

## Removing not interesting columns 
final_results.drop(columns=['source', 'significant', 'evidences', 'parents','description'], inplace=True)
final_results.to_excel(OUTPUT_DATA / 'GO_all_clusters.xlsx')

# %%  Plotting all GO in one scatterplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
top = 7
## Only select top ontologies per cluster 
top_per_cluster = (
    final_results
    .groupby('cluster', group_keys=False)
    .apply(lambda x: x.nlargest(top, 'fold_enrichment'))
)

top_per_cluster

# Scaling factor for marker sizes
s_incr = 30
s_values = top_per_cluster['intersection_size'] * s_incr
s_values = pd.to_numeric(s_values)
# # Optional: clip extremely large marker sizes for readability
# s_values = np.clip(s_values, 20, 500)
 
# Create figure
fig, ax = plt.subplots(figsize=(4.5, 20), dpi=300)
ax.grid(axis='x', linestyle='--', alpha=0.6)
ax.grid(False, axis='y')

# Scatter plot
sc = ax.scatter(
    x=top_per_cluster['cluster'],
    y=top_per_cluster['name'],
    c=top_per_cluster['-log10_pval'],
    cmap='magma',
    s=s_values,
    edgecolor='k',
    linewidth=1,
    alpha=0.7,
    zorder=3
)


# ---- FIX 1: X-axis clipping ----
# x_vals = np.sort(final_results['window_index'].unique())
# ax.set_xlim(x_vals.min() - 1, x_vals.max() + 0.5)

# ax.set_xticks(x_vals)
# ax.set_xticklabels(x_vals, fontsize=12)

# ---- FIX 2: Y-axis excess spacing ----
y_vals = top_per_cluster['name'].unique()
ax.set_ylim(-0.5, len(y_vals) - 0.5)

ax.tick_params(axis='y', labelsize=15)
# Dynamic size legend based on actual marker sizes
legend_sizes = [s_values.min(), s_values.mean(), s_values.max()]
for size in legend_sizes:
    ax.scatter([], [], c='darkgrey', alpha=1, s=size, label=str(int(size / s_incr)))
# Move legend outside to the right
ax.legend(
    title='Number of genes',
    scatterpoints=1,
    frameon=True,
    fontsize=12,
    bbox_to_anchor=(1.1, 1),  # x,y coordinates outside axes
    loc='upper left'
)
# Colorbar
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('-log10(FDR)', fontsize=30)
cbar.ax.tick_params(labelsize=30)  # increase tick numbers


# Formatting x-axis
# ax.set_xticks(sorted(final_results['window_index'].unique()))
# ax.set_xticklabels(sorted(final_results['window_index'].unique()),  fontsize=12)
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15, rotation = 90)



# ax.set_xlabel('clusters', fontsize=14, rotation = 90)1
ax.set_ylabel('GO term', fontsize=14)
ax.set_title('GO enrichment across windows', fontsize=16)

# Save and show
plt.tight_layout()
plt.savefig(
    OUTPUT_FIGURES / 'GO_all_clusters.pdf',
    bbox_inches='tight'
)
plt.show()


# %%  ## Replacing AGI numbers by gene name

final_results_names = pd.DataFrame()
for row in list(range(len(final_results))):
    temp1 = final_results.iloc[row, :]
    temp1['names'] = names_changes_list(temp1['intersections'])
    temp1 = temp1.to_frame().T
    final_results_names = pd.concat([final_results_names, temp1], axis=0)

final_results_names.to_excel(OUTPUT_DATA / 'GO_all_clusters_genes_names.xlsx')

# %% # PLOT PER CLUSTER 

from matplotlib.lines import Line2D

go_df = pd.DataFrame()
for cluster in final_results['cluster'].unique().tolist(): 
    temp1 = final_results[final_results['cluster'] == cluster]
    temp1 = temp1[0:10]
    s_incr = 50
    s_values = temp1['intersection_size'] * s_incr
    s_values = pd.to_numeric(s_values)
    fig, ax = plt.subplots(figsize=(3, 5), dpi = 300)   
    ax.grid(False)       
    
    sc = ax.scatter(data=temp1,  x='-log10_pval',
                    y='name',
                    c='fold_enrichment',
                    cmap='inferno',
                    s=s_values,
                    edgecolor='k',
                    linewidth=1,
                    alpha=0.8,
                    zorder=3)
    # Colorbar
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Fold Enrichment', fontsize=15)
    ax.tick_params(axis='both', labelsize=15, width=2)  # thicker tick marks
    plt.title(str(cluster), fontsize = 15)
    plt.xlabel('-log10_pval', fontsize=15)
    # plt.savefig(str(cluster) + '__GO.pdf', bbox_inches='tight')
    plt.show()
