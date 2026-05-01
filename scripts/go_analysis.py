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
# from imgmisc import listdir

os.chdir("/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/figures/GO_clusters/")
from glob import glob

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


#%% Gene ontology for all clusters 
import pandas as pd
from gprofiler import GProfiler
import numpy as np
trajectory = 'EP'
cc_var = pd.read_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx',
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


final_results.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_all_clusters.xlsx')

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
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/GO_all_clusters.pdf',
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

final_results_names.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_all_clusters_genes_names.xlsx')

# %% # PLOT PER CLUSTER 

from matplotlib.lines import Line2D


go_df = pd.DataFrame()
for cluster in final_results['cluster'].unique().tolist(): 
    temp1 = final_results[final_results['cluster'] == cluster]
    # temp1['cluster'] = f.split('_')[0]
    # temp1['-log(FDR)'] = np.log(temp1['Enrichment FDR'])*-1
    # temp1 = temp1[temp1['Enrichment FDR'] <= 0.05]
    # temp1 = temp1.sort_values('Fold Enrichment', ascending = False)
    # temp1 = temp1.sort_values('-log(FDR)', ascending = False)
    # temp1['Pathway'] = temp1['Pathway']
    temp1 = temp1[0:10]
    # temp1['label'] = temp1['Pathway'] + '   ' + temp1['cluster']
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


     # ---------- FIXED SIZE LEGEND ----------
    # sizes_for_legend = sorted(temp1['intersection_size'].unique())

    # legend_elements = [
    #     Line2D(
    #         [0], [0],
    #         marker='o',
    #         linestyle='',
    #         label=str(int(size)),
    #         markerfacecolor='darkgrey',
    #         markeredgecolor='k',
    #         markersize=np.sqrt(size * s_incr)  # area → diameter
    #     )
    #     for size in sizes_for_legend
    # ]

    # legend1 = ax.legend(
    #     handles=legend_elements,
    #     title='nGenes',
    #     loc='center left',
    #     bbox_to_anchor=(1.1, 1),  # outside, right side
    #     frameon=False
    # )

    # ax.add_artist(legend1)
    # Colorbar
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Fold Enrichment', fontsize=15)
    ax.tick_params(axis='both', labelsize=15, width=2)  # thicker tick marks
    plt.title(str(cluster), fontsize = 15)
    plt.xlabel('-log10_pval', fontsize=15)
    # plt.savefig(str(cluster) + '__GO.pdf', bbox_inches='tight')
    plt.show()

# %% 

#     go_df = pd.concat([go_df,temp1[0:5]])  ### number of top GO to plot later 
    
# # go_df['Pathway'] = go_df['Pathway'].astype(str)
# # go_df = go_df.dropna(subset=['Pathway'])

  


# ### Plotting all GO in one scatterplot 
# # Create figure and axis
# fig, ax = plt.subplots(figsize=(5, 20), dpi = 1000)   
# # Add grid behind scatterplot
# ax.grid(axis='x', zorder=0)
# # Scatter plot
# s_incr = 25     
# sc = ax.scatter(data=go_df, x='cluster', y='Pathway', 
#                 c='-log(FDR)', 
#                 cmap='inferno', 
#                 s=go_df['nGenes'] * s_incr, 
#                 edgecolor='k', linewidth=0.5, alpha=0.95, zorder=3)

# # Add size legend for number of genes
# for size in range(int(go_df['nGenes'].min() * s_incr), int(go_df['nGenes'].max() * s_incr), int(s_incr)):
#     ax.scatter([], [], c='darkgrey', alpha=1, s=size, label=str(int(size / s_incr)))

# # Colorbar
# cbar = plt.colorbar(sc, ax=ax)
# cbar.set_label('-log(FDR)', fontsize=12)

# # Formatting
# ax.set_xticks(go_df['cluster'].unique())  # Ensure x-ticks are correct
# ax.set_xticklabels(go_df['cluster'].unique(), rotation=90, fontsize=15)
# ax.tick_params(axis='y', labelsize=15)

# # Save and show
# plt.savefig(INPUT_DIR + 'ALL-CLUSTERS__GO.pdf', bbox_inches='tight')
# plt.show()
    
    
    
# # ##### GO for TARGETS OF DOF5.8 and MP 


# # INPUT_DIR = "/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/figures/GO_clusters/"
# # # from snuclei_run_def import scrublet_processing, qc_visualization, qc_processing , scar_ambient_denoising , scvi_doublet_removal, names_changes_list
# # # from snuclei_run_def import separate_floral_meristem, change_names_adata ,  hvg_pca, names_change_dataframe, umap, regressing_out_cell_cycle, population_per_domain , removing_outlier_clusters
# # # from snuclei_run_def import coexpression_cluster, coexpression_arboreto, coexpression_arboreto_by_cluster, umap_3d
# # # plt.rcParams["font.family"] = "Arial"
# # sc.set_figure_params(dpi=100)

# # ## PLOT PER CLUSTER 
# # files = os.listdir(INPUT_DIR)
# # files = [ x for x in files if not "pdf"  in x ] 
# # files = [ x for x in files if not "xlsx"  in x ] 
# # files = [ x for x in files if  "target_genes_GO.csv"  in x ] 



# # go_df = pd.DataFrame()

# # for f in files: 
# #     temp1 = pd.read_csv(f)
# #     temp1['cluster'] = f.split('_')[0]
# #     temp1['-log(FDR)'] = np.log(temp1['Enrichment FDR'])*-1
# #     temp1 = temp1[temp1['Enrichment FDR'] <= 0.05]
# #     temp1 = temp1.sort_values('Fold Enrichment', ascending = False)
# #     temp1['Pathway'] = temp1['Pathway'].str.extract(r'(GO:\d+)\s+(.*)')[1]
       
# #     temp1['label'] = temp1['Pathway'] + '   ' + temp1['cluster']
# #     go_df = pd.concat([go_df,temp1[0:25]])  ### number of top GO to plot later 
       
# #     temp1 = temp1.sort_values('-log(FDR)', ascending = True)
# #     # if len(temp1) > 0: 
# #     # tissue_color = cm.get_cmap('viridis')   
# #     plt.figure(figsize=(3,3))   
# #     plt. grid(False)
# #     s_incr = 20     
# #     sc = plt.scatter(data = temp1, x = '-log(FDR)', y ='Pathway' , c = 'Fold Enrichment',
# #                 cmap = 'cividis', s = temp1['nGenes']*s_incr,
# #                 edgecolor='k', linewidth=1, alpha=1)
# #     # plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)
# #     for size in list(range(temp1.nGenes.min()*s_incr, temp1.nGenes.max()*s_incr,temp1.nGenes.min()*s_incr*2 ))[0::]:
# #         plt.scatter([], [], c='grey', alpha=1, s=size,
# #                        label=str(size/s_incr))
# #     plt.legend(scatterpoints=1, frameon=False,
# #                   labelspacing=0.5, title='Number of genes per GO', loc='upper center', bbox_to_anchor=(1, 1.3), ncol=10)
# #     # for c in temp1.cluster.unique(): 
# #     #     temp3 = temp1[temp1["cluster"] == c]
# #     #     # print(temp3, c)
# #     #     line = max(temp3.index)
# #     #     # plt.hlines(line+0.5,xmin=10, xmax=temp1['Fold Enrichment'].max()+10, color="black", linestyle="--" ) 
# #     plt.colorbar(sc, label = 'Fold Enrichment' )
# #     # print(temp1['Fold Enrichment'].max()+10, f)
# #     # plt.xlim(0 , temp1['-log(FDR)'].max()+10)
# #     plt.title(str(f))
# #     plt.xlabel('-log(FDR)')
# #     plt.savefig( INPUT_DIR + str(f) + '_GO.pdf', bbox_inches='tight')
# #     plt.show()        
# #     temp3 = pd.DataFrame()
# #     for i in list(range(len(temp1))): 
# #         temp2 = temp1.iloc[i,::]
# #         genes = temp2['Genes'].split(' ')
# #         names = names_changes_list(genes)
# #         temp2['AGI_name'] = names
# #         temp3 = pd.concat([temp3,temp2], axis=1)
# #     temp3.to_excel(INPUT_DIR  +str(f)+ ".xlsx")
# #     # go_df = pd.concat[(go_df, temp1)]
    


# # ### Plotting all GO in one scatterplot 
# # plt.figure(figsize=(1,5))   
# # plt. grid(False)
# # s_incr = 10     
# # ax = plt.scatter(data = go_df, x = 'cluster', y ='Pathway' , c = '-log(FDR)', 
# #             cmap = 'RdBu_r', s = go_df['nGenes']*s_incr,
# #             edgecolor='k', linewidth=0.5, alpha=1)
# # for size in list(range(int(go_df['nGenes'].min() * s_incr), int(go_df['nGenes'].max() * s_incr), int(s_incr))):
# #     plt.scatter([], [], c='darkgrey', alpha=1, s=size, label=str(int(size / s_incr)))

# # plt.legend(scatterpoints=1, frameon=False,
# #               labelspacing=0.5, title='Number of genes', loc='upper center', bbox_to_anchor=(1, 1.3), ncol=10)
# # plt.colorbar(ax, label = '-log(FDR)' )
# # plt.xticks(rotation=90, ha='right', fontsize = 7)
# # plt.yticks( fontsize = 7)
# # plt.margins(0.5,0.1)
# # plt.savefig( INPUT_DIR +  'EP_MP_and_DOF5.8_targets_GO.pdf', bbox_inches='tight')
# # plt.show()        
   











