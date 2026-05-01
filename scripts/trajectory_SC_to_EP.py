#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:30:32 2026

@author: Sebastian
"""

import os, sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"

import scFates as scf

import scanpy as sc
import warnings
import numpy as np
import pandas as pd
import palantir
import seaborn as sns
warnings.filterwarnings("ignore")
import sys
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout   


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

def scfates_trajectories_alignment(adata, clusters, root_gene, nodes, sigma, eigen):
    # Select clusters
    if len(clusters) == 1:
        adata = adata[adata.obs['leiden'] == clusters[0]]
    if len(clusters) > 1:
        adata = adata[adata.obs['leiden'].isin(clusters)]
    # Learn curve using ElPiGraph algorithm
    # sc.pl.pca(adata,color="leiden",cmap="tab20")
    pca_projections = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
    # Run Palantir to obtain multiscale diffusion space
    palantir.utils.run_diffusion_maps(adata,n_components=20)
    ms_data = palantir.utils.determine_multiscale_space(adata, n_eigs=eigen)

    # MAGIC imputation 
    imputed_X = palantir.utils.run_magic_imputation(adata)
    adata.layers["magic"] = imputed_X

    # generate neighbor draph in multiscale diffusion space
    adata.obsm["X_palantir"] = ms_data.values
    sc.pp.neighbors(adata, n_neighbors=100,
                    use_rep="X_palantir", method='umap')
    # draw ForceAtlas2 embedding using 2 first PCs as initial positions
    adata.obsm["X_pca2d"] = adata.obsm["X_pca"][:, 0:2]
    sc.tl.draw_graph(adata, init_pos='X_pca2d')
    title = str(nodes) + ' number of nodes ' + str(sigma) + \
        ' sigma value_eigen: ' + str(eigen)
    sc.set_figure_params(dpi_save=600)
    sc.pl.draw_graph(adata, color=['leiden'],  add_outline=True,  legend_fontsize=10, legend_fontoutline=2,
                     title=title,
                     palette=sns.color_palette('colorblind'), save=title + '.pdf')
    sc.pl.draw_graph(adata, color=[root_gene],  add_outline=True,  legend_fontsize=10, legend_fontoutline=2,
                     cmap='magma', save="_FA_diffusion_trajectory_root_gene.pdf")
    scf.tl.tree(adata, method="ppt", Nodes=nodes, use_rep="palantir", plot=False,
                device="cpu", seed=1,  ppt_lambda=100, ppt_nsteps=200,
                ppt_sigma=sigma)
    scf.pl.graph(adata, title=str(nodes) + ' number of nodes ' + str(sigma) + ' sigma value', color_cells='leiden',
                 palette=sns.color_palette('colorblind'), save="_NODES_FA_diffusion_trajectory_root_gene.pdf")
    return adata

def scfates_trajectories_dendogram(adata, root_segment):
    # Selecting root from graph  and pseudotime
    scf.tl.root(adata, root_segment)
    scf.tl.pseudotime(adata, n_jobs=10, n_map=10, seed=42)
    scf.pl.trajectory(adata, color_seg="t", basis="draw_graph_fa",
                      frameon=False, s=100, scale_path=.6, save="_pseudotime_paths.pdf")
    # scf.pl.trajectory(adata,color_seg="seg",basis="draw_graph_fa",frameon=False,s=100,scale_path=.6,save="_seg_paths.png",
    #                   palette=sns.color_palette('colorblind'))
    sc.pl.draw_graph(adata, color=["milestones", 'leiden'], palette=sns.color_palette(
        'colorblind'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2)
    scf.tl.dendrogram(adata)
    scf.pl.dendrogram(adata, color="seg", palette=sns.color_palette('colorblind'), legend_fontoutline=True, legend_loc="on data",
                      save="_seg.pdf")
    scf.pl.dendrogram(adata, color="milestones", palette=sns.color_palette('colorblind'), legend_fontoutline=True, legend_loc="on data",
                      save="_milestones.pdf")

    return adata


# %%  scfates pipeline to obtain DEG for internal layers  
adata = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_filtered_norm_clusters.h5ad')  ### WITH CELL CYCLE CLUSTERS 
# adata = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_filtered_norm_clusters_without_cell_cycle.h5ad')  ### WITH CELL CYCLE CLUSTERS =
adata.obs['new_cell_identity'] = adata.obs['leiden']  ## Create this new obs column for processing 

for eigen in [4]:
        for nodes in [200]:
            for sigma in [0.05]:
                adata1 = scfates_trajectories_alignment(adata, clusters = [
                                                                      'undifferentiated cells',
                                                                      'EP',
                                                                     ],
                                                        eigen = eigen,
                                                        root_gene = 'AT2G17950', 
                                                        nodes = nodes, 
                                                        sigma = sigma,
                                                        )

sc.set_figure_params(figsize=(3,3),dpi=300,frameon=True)
sc.pl.draw_graph(adata1,color=["leiden"], palette=['rebeccapurple', 'pink'], add_outline=True,  
                 legend_fontsize=10, legend_fontoutline=2, size = 250, 
                     save="_seg_leiden.pdf")  

sc.set_figure_params(figsize=(3,3),dpi_save=300,frameon=False)
sc.pl.draw_graph(adata1, color = ['AT2G17950','AT2G27250' ,'AT1G80100'],  
                 layer = 'magic',
                 add_outline=True, legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = ['WUS', 'CLV3', 'AHP6'])



# %% Running milestone and pseudotime defininig root tip using CLV3 and WUS expression levels
adata2 = scfates_trajectories_dendogram(adata1, 
                                        root_segment = 18)  

adata2.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_SC-to-EP_layers.h5ad')    ## save with new cluster names 


# %%  chage milestones names 
adata2 = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_SC-to-EP_layers.h5ad')

milestone = pd.DataFrame(adata2.obs.groupby('milestones'))[0].tolist()
new_names_milestone = ['EP', 'SC']

scf.tl.rename_milestones(adata2,new_names_milestone)
sc.pl.draw_graph(adata2,color='milestones',color_map = 'tab20', add_outline=True, size =200, legend_fontsize=10, legend_fontoutline=2)             

from scipy.stats import rankdata
sc.set_figure_params(figsize=(3,4),dpi_save=600,frameon=False)
adata_Ic=scf.tl.subset_tree(adata2,root_milestone='SC',milestones=['EP'],copy=True)
sc.pl.draw_graph(adata_Ic,color=["leiden"],frameon=True, palette=sns.color_palette('Dark2'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2)
sc.pl.draw_graph(adata_Ic,color=["milestones"], palette=sns.color_palette('colorblind'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2) 
pt = adata_Ic.obs['t'].values
# # Step 2: Rank and rescale to [0, 1]
pt_smooth = rankdata(pt, method='average') / len(pt)
# # Step 3: Store in AnnData
adata_Ic.obs['t'] = pt_smooth
sc.pl.draw_graph(adata_Ic,color=["t"], palette=sns.color_palette('viridis'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2) 
scf.tl.test_association(adata_Ic,n_jobs=1,A_cut=.3)
adata_Ic.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_SC-to-EP_layers'  + 'EP' + '_association.h5ad') 
adata_Ic = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_SC-to-EP_layers'  + 'EP' + '_association.h5ad')
adata_Ic.var['signi'] =  (adata_Ic.var['p_val'] < 0.01) & (adata_Ic.var['A'] > 0)  ### Use this filter instead the default one 
scf.tl.fit(adata_Ic,n_jobs=1)
adata_Ic.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_SC-to-EP_layers'  + 'EP' + '_association_fitted.h5ad') 
           

# %% Importing anndata with DEG per trajectory only
trans = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_SC-to-EP_layers'  + 'EP' + '_association_fitted.h5ad')
     
from scipy.stats import pearsonr
trans.var["corr"]=list(map(lambda g: pearsonr(trans.obs.t,trans[:,g].layers["fitted"].flatten())[0],trans.var_names))
trans.var["up"]=trans.var["corr"]>0
# cortex.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/stem_cortex_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( trans[:, trans.var_names].layers["fitted"], index=trans.obs_names, columns=trans.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: trans.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = trans.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/SCtoEP_DEG_ordered_peak_expression.xlsx")
                
                    
# %% 
from pathlib import Path
import os

project_dir = Path("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024")
os.chdir(project_dir)

trans.var.index = names_changes_list(trans.var['gene_ids'])

#%% 
sc.set_figure_params(figsize=(2,10),dpi_save=600,frameon=False)
scf.pl.trends(trans,
                highlight_features = [
                                        'CKX3',
                                       'WUS', 
                                        'GH3.3',
                                        'CLV3', 
                                        'YAB1', 
                                        'IAA30',                        
                                        'AHP6'
                                            ],
                layer="magic" , 
                # n_features = len(phloem), 
              style="italic",
              add_outline=True,
              basis="dendro",
              # feature_cmap = 'inferno', 
              show_segs = True, 
              fontsize = 10,
              figsize = (3,5),
                ordering="max",
                ord_thre=0.99,
              save="SC-to-EP_DEF_selected.pdf")           
         

# %% single-gene trends
sc.set_figure_params(figsize=(2,10),dpi=300,frameon=False)
scf.pl.single_trend(trans,'WUS',ylab= 'WUS expression',layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (10,3))
scf.pl.single_trend(trans,'CLV3',ylab= 'CLV3 expression',layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (10,3))
scf.pl.single_trend(trans,'AHP6',ylab= 'AHP expression',layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (10,3))

sc.pl.draw_graph(trans,color=["t"],  add_outline=True,  legend_fontsize=10, legend_fontoutline=2) 


#%% Gene ontology analysis over cell cycle trajectory (dataframe)
import pandas as pd
from gprofiler import GProfiler
import numpy as np

cc_var = pd.read_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/SCtoEP_DEG_ordered_peak_expression.xlsx')

# Instantiate the object
gp = GProfiler(return_dataframe=True)

window_size = int(len(cc_var) * 0.10)
num_windows = (len(cc_var) + window_size - 1) // window_size

all_results = []  # To store results of all windows
for i in range(num_windows):
    start = i * window_size
    end = start + window_size
    window_df = cc_var.iloc[start:end]
    
    # Extract gene IDs (TAIR IDs)
    gene_list = window_df['gene_ids'].tolist()
    
    # Run GO enrichment for Arabidopsis
    results = gp.profile(
        organism='athaliana',
        query=gene_list,
        significance_threshold_method='fdr', 
        sources=['GO:BP'], 
        user_threshold=0.05, 
        no_evidences=False
    )
       
    # Add a column to keep track of the window
    results['window_index'] = i + 1
    all_results.append(results)

final_results = pd.concat(all_results, ignore_index=True)
final_results['-log10_pval'] = -np.log10(final_results['p_value'])
## Incresing threhold? 
final_results = final_results[(final_results['p_value'] < 0.01) & (final_results['intersection_size'] > 2) & 
                               (final_results['term_size'] < 1000)]

## Removing not interesting columns 
final_results.drop(columns=['source', 'significant', 'evidences', 'parents','description'], inplace=True)


final_results.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_SCtoEP_DEG_ordered_peak_expression_genes.xlsx')

# %%  Plotting all GO in one scatterplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Ensure intersection_size is numeric and drop NaNs
# final_results['intersection_size'] = pd.to_numeric(final_results['intersection_size'], errors='coerce')
# final_results = final_results.dropna(subset=['intersection_size', 'window_index', 'name', '-log10_pval'])

# Scaling factor for marker sizes
s_incr = 40
s_values = final_results['intersection_size'] * s_incr
s_values = pd.to_numeric(s_values)
# # Optional: clip extremely large marker sizes for readability
# s_values = np.clip(s_values, 20, 500)
 
# Create figure
fig, ax = plt.subplots(figsize=(2, 12), dpi=300)
ax.grid(axis='x', zorder=0)

# Scatter plot
sc = ax.scatter(
    x=final_results['window_index'],
    y=final_results['name'],
    c=final_results['-log10_pval'],
    cmap='inferno',
    s=s_values,
    edgecolor='k',
    linewidth=1,
    alpha=0.8,
    zorder=3
)


# ---- FIX 1: X-axis clipping ----
x_vals = np.sort(final_results['window_index'].unique())
ax.set_xlim(x_vals.min() - 1, x_vals.max() + 0.5)

ax.set_xticks(x_vals)
ax.set_xticklabels(x_vals, fontsize=12)

# ---- FIX 2: Y-axis excess spacing ----
y_vals = final_results['name'].unique()
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
ax.set_xticks(sorted(final_results['window_index'].unique()))
ax.set_xticklabels(sorted(final_results['window_index'].unique()),  fontsize=12)
ax.tick_params(axis='y', labelsize=15)

ax.set_xlabel(str(window_size) + ' DEG GO Pseudtime from SC to EP', fontsize=14)
ax.set_ylabel('GO term', fontsize=14)
ax.set_title('GO enrichment across windows', fontsize=16)

# Save and show
plt.tight_layout()
plt.savefig(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/GO_SCtoEP_DEG_ordered_peak_expression_genes.pdf',
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

final_results_names.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_SCtoEP_DEG_ordered_peak_expression_genes_names.xlsx')



# %% hypoxia-response from SC to EP
sc.set_figure_params(dpi=300,frameon=False)
scf.pl.single_trend(trans,'TIL',ylab= 'TIL expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))
scf.pl.single_trend(trans,'SUS1',ylab= 'SUS1 expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))
scf.pl.single_trend(trans,'AT3G10040',ylab= 'HRA1 expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))
scf.pl.single_trend(trans,'AT4G27450',ylab= 'HUP54 expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))
scf.pl.single_trend(trans,'GDPD2',ylab= 'GDPD2 expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))

scf.pl.single_trend(trans,'EIF3G1',ylab= 'EIF3G1 expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))
scf.pl.single_trend(trans,'KAN1',ylab= 'KAN1 expression',plot_emb =  False,layer ='magic',color_exp = 'indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.02,  figsize = (3,3))





sc.set_figure_params(figsize=(3,3),dpi_save=300,frameon=False)
sc.pl.draw_graph(trans, color = ['AT3G10040','AT4G27450', 'TIL','SUS1' ,],  
                 layer = 'magic',
                 add_outline=True, legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = [ 'HRA1','HUP54' ,'TIL', 'SUS1',])





