#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 17:06:16 2023

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

for eigen in [3,5,7,9]:
        for nodes in [200]:
            for sigma in [0.05]:
                adata1 = scfates_trajectories_alignment(adata, clusters = [
                                                                      'undifferentiated cells',
                                                                      'xylem parenchyma',
                                                                      'procambium',
                                                                      'EP',
                                                                      # 'rib zone',
                                                                      'cortex',
                                                                      # 'BD'
                                                                     ],
                                                        eigen = eigen,
                                                        root_gene = 'AT2G17950', 
                                                        nodes = nodes, 
                                                        sigma = sigma,
                                                        )

sc.pl.draw_graph(adata1,color=["new_cell_identity"], palette=sns.color_palette('tab20'), add_outline=True,  
                 legend_fontsize=10, legend_fontoutline=2, size = 50, 
                     save="_seg_leiden.pdf")  

sc.pl.draw_graph(adata1,color=["leiden"], palette=sns.color_palette('Accent'), add_outline=True,  
                 legend_fontsize=10, legend_fontoutline=2, size = 100, 
                     save="_seg_leiden.pdf")  



# %% Running milestone and pseudotime defininig root tip using CLV3 and WUS expression levels
adata2 = scfates_trajectories_dendogram(adata1, 
                                        root_segment = 145) 

# %%  chage milestones names 
milestone = pd.DataFrame(adata2.obs.groupby('milestones'))[0].tolist()
new_names_milestone = ['xylem', 'EP', 'xylem parenchyma', '136', 'SC', 'cambium', '18', '184', 'phloem', '53', '63', 'cortex']
scf.tl.rename_milestones(adata2,new_names_milestone)
sc.pl.draw_graph(adata2,color='milestones',palette= 'tab20', add_outline=True, size =200, legend_fontsize=10, legend_fontoutline=2,
                  )

# %% Plotting specific genes for main figure 6  

sc.set_figure_params(figsize=(4,4),frameon=False, dpi=200)
sc.pl.draw_graph(adata2, color = ['AT2G17950'], layer = 'magic',  add_outline=True, legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'WUS')
sc.pl.draw_graph(adata2, color = ['AT2G27250'], layer = 'magic', add_outline=True, legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'CLV3')
sc.pl.draw_graph(adata2, color = ['AT4G32880'],layer = 'magic',  add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size = 300, 
                cmap = 'RdGy_r', title = 'ATHB8')
sc.pl.draw_graph(adata2, color = ['AT1G80100'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'AHP6')
sc.pl.draw_graph(adata2, color = ['AT3G14570'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'CALS8')
sc.pl.draw_graph(adata2, color = ['AT5G61480'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 200, 
                cmap = 'RdGy_r', title = 'PXY')
sc.pl.draw_graph(adata2, color = ['AT5G19260'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 200, 
                cmap = 'RdGy_r', title = 'FAF3')
sc.pl.draw_graph(adata2, color = ['AT5G03150'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 200, 
                cmap = 'RdGy_r', title = 'JKD')

adata2.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers.h5ad')    ## save with new cluster names 

# %%  Plotting segments figures 
adata2 = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers.h5ad')

sc.pl.draw_graph(adata2,color=["milestones"], palette=sns.color_palette('tab20'), add_outline=True,  
                 legend_fontsize=10, legend_fontoutline=2, size = 50)

### REPLACE MILESTONE FOR SOMETHING MORE MEANINGFUL 
seg = adata2.obs['milestones'].tolist()

seg = ['procambium' if cell == '136' else cell for cell in seg]
seg = ['procambium' if cell == '18' else cell for cell in seg]
seg = ['Trans-amplifying cells' if cell == '184' else cell for cell in seg]
seg = ['xylem parenchyma' if cell == '63' else cell for cell in seg]
seg = ['cortex' if cell == '53' else cell for cell in seg]

adata2.obs['milestones'] = seg

sc.set_figure_params(figsize=(7,7),frameon=False, dpi=300)
sc.pl.draw_graph(adata2,color=["milestones"], palette='PRGn', add_outline=True, 
                 size = 200, legend_fontsize=10, legend_fontoutline=2, 
                 save = 'segments_trajectory.pdf')
                 
# %% 
leiden = adata2.obs['leiden'].tolist()
leiden = ['EP 1' if cell == 'new' else cell for cell in leiden]
leiden = ['EP 2' if cell == 'EP' else cell for cell in leiden]
leiden = ['Cortex' if cell == 'cortex' else cell for cell in leiden]
leiden = ['Stem cells' if cell == 'iSC' else cell for cell in leiden]
leiden = ['Vasculature 1' if cell == 'vasculature 1' else cell for cell in leiden]
leiden = ['Vasculature 2' if cell == 'vasculature 2' else cell for cell in leiden]
adata2.obs['leiden'] = leiden

import matplotlib.colors as mcolors
 # Get the original tab20 colors
tab20_palette = sns.color_palette("tab20", 20)
# Function to darken colors
def darken_color(color, factor=0.6):
    rgb = np.array(mcolors.to_rgb(color))  # Convert to NumPy array
    return np.clip(rgb * factor, 0, 1)  # Multiply and clip values between 0 and 1
# Apply the darkening function
dark_tab20 = [darken_color(c, factor=0.9) for c in tab20_palette]
# Convert back to hex
dark_tab20_hex = [mcolors.to_hex(c) for c in dark_tab20]

sc.set_figure_params(figsize=(5,5),frameon=False, dpi=300)
sc.pl.draw_graph(adata2,color=["new_cell_identity"], palette=dark_tab20_hex, add_outline=True, 
                 size = 100, legend_fontsize=10, legend_fontoutline=2, frameon = True, 
                 save = 'leiden_trajectory.pdf')

# adata2.write_h5ad(path + '/adata_scfates_' + str(name_file) + '.h5ad')    ## save with new cluster names 

# %% Plotting DEG per cell fate 
adata2 = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers.h5ad')

from scipy.stats import rankdata

### only DEG per branch
for milestone in [ 'xylem parenchyma', 'xylem', 'cambium', 'phloem', 'cortex', 'EP']: 
    sc.set_figure_params(figsize=(3,4),dpi_save=600,frameon=False)
    adata_Ic=scf.tl.subset_tree(adata2,root_milestone='SC',milestones=[milestone],copy=True)
    sc.pl.draw_graph(adata_Ic,color=["leiden"],frameon=True, palette=sns.color_palette('Dark2'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2)
    sc.pl.draw_graph(adata_Ic,color=["milestones"], palette=sns.color_palette('colorblind'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2) 
    pt = adata_Ic.obs['t'].values
    # Step 2: Rank and rescale to [0, 1]
    pt_smooth = rankdata(pt, method='average') / len(pt)
    # Step 3: Store in AnnData
    adata_Ic.obs['t'] = pt_smooth
    # adata_Ic = scfates_trajectories_dendogram(adata_Ic,  root_segment=131)
    sc.pl.draw_graph(adata_Ic,color=["t"], palette=sns.color_palette('viridis'), add_outline=True,  legend_fontsize=10, legend_fontoutline=2) 
    scf.tl.test_association(adata_Ic,n_jobs=1,A_cut=.3)
    adata_Ic.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + str(milestone) + '_association.h5ad') 
  
    
for milestone in [ 'xylem parenchyma', 'xylem', 'cambium', 'phloem', 'cortex', 'EP']:    
    adata_Ic = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + str(milestone) + '_association.h5ad')
    adata_Ic.var['signi'] =  adata_Ic.var['p_val'] < 0.01  ### Use this filter instead the default one 
    scf.tl.fit(adata_Ic,n_jobs=1, gamma = 1)
    adata_Ic.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + str(milestone) + '_association_fitted.h5ad') 

# %% Importing anndata with DEG per trajectory only
cortex = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + 'cortex' + '_association_fitted.h5ad')
phloem = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + 'phloem' + '_association_fitted.h5ad')
xylem = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + 'xylem' + '_association_fitted.h5ad')
cambium = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + 'cambium' + '_association_fitted.h5ad')
EP = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + 'EP' + '_association_fitted.h5ad')
xylem_par = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers_'  + 'xylem parenchyma' + '_association_fitted.h5ad')


# %% pearson correlation for each trajectory for posterior ploting by max intensity 
from scipy.stats import pearsonr

# Cortex 
cortex.var["corr"]=list(map(lambda g: pearsonr(cortex.obs.t,cortex[:,g].layers["fitted"].flatten())[0],cortex.var_names))
cortex.var["up"]=cortex.var["corr"]>0
# cortex.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/stem_cortex_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( cortex[:, cortex.var_names].layers["fitted"], index=cortex.obs_names, columns=cortex.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: cortex.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = cortex.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/cortex_DEG_ordered_peak_expression.xlsx")

# Xylem 
xylem.var["corr"]=list(map(lambda g: pearsonr(xylem.obs.t,xylem[:,g].layers["fitted"].flatten())[0],xylem.var_names))
xylem.var["up"]=xylem.var["corr"]>0
# xylem.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/xylem_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( xylem[:, xylem.var_names].layers["fitted"], index=xylem.obs_names, columns=xylem.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: xylem.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = xylem.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/xylem_DEG_ordered_peak_expression.xlsx")

# Xylem parenchyma
xylem_par.var["corr"]=list(map(lambda g: pearsonr(xylem_par.obs.t,xylem_par[:,g].layers["fitted"].flatten())[0],xylem_par.var_names))
xylem_par.var["up"]=xylem_par.var["corr"]>0
# xylem_par.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/xylem_ves_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( xylem_par[:, xylem_par.var_names].layers["fitted"], index=xylem_par.obs_names, columns=xylem_par.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: xylem_par.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = xylem_par.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/xylemparenchyma_DEG_ordered_peak_expression.xlsx")

# EP 
EP.var["corr"]=list(map(lambda g: pearsonr(EP.obs.t,EP[:,g].layers["fitted"].flatten())[0],EP.var_names))
EP.var["up"]=EP.var["corr"]>0
# EP.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/EP_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( EP[:, EP.var_names].layers["fitted"], index=EP.obs_names, columns=EP.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: EP.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = EP.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/EP_DEG_ordered_peak_expression.xlsx")

# Phloem 
phloem.var["corr"]=list(map(lambda g: pearsonr(phloem.obs.t,phloem[:,g].layers["fitted"].flatten())[0],phloem.var_names))
phloem.var["up"]=phloem.var["corr"]>0
# phloem.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/phloem_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( phloem[:, phloem.var_names].layers["fitted"], index=phloem.obs_names, columns=phloem.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: phloem.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = phloem.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/phloem_DEG_ordered_peak_expression.xlsx")

# Cambium 
cambium.var["corr"]=list(map(lambda g: pearsonr(cambium.obs.t,cambium[:,g].layers["fitted"].flatten())[0],cambium.var_names))
cambium.var["up"]=cambium.var["corr"]>0
# cambium.var.sort_values(by=['corr'], ascending=True).to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/cambium_DEG_ordered_peak_expression.xlsx")
fitted = pd.DataFrame( cambium[:, cambium.var_names].layers["fitted"], index=cambium.obs_names, columns=cambium.var_names).T.copy(deep=True)
feature_order = ( fitted.apply( lambda x: cambium.obs.t[fitted.columns][(x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
df_bla = cambium.var.loc[feature_order]
df_bla['gene_name'] = names_changes_list(df_bla['gene_ids'])
df_bla.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/cambium_DEG_ordered_peak_expression.xlsx")


# %% PHLOEM trajectry plotting 
from pathlib import Path
import os

project_dir = Path("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024")
os.chdir(project_dir)
sc.set_figure_params(figsize=(2,10),dpi_save=600,frameon=False)
phloem.var.index = names_changes_list(phloem.var['gene_ids'])
scf.pl.trends(phloem,
                highlight_features = [
                                        'CLV1',
                                      'WUS', 
                                        # 'WAT1',
                                        'CALS6', 
                                        'OBP2',  
                                        # 'SWEET10',
                                        'SMXL5',
                                        'BAM3',
                                        'CALS8',
                                        'APL', 
                                        'DOF5.6',
                                        'CALS7' ,                        
                                        'DOF5.1',
                                            ],
                # annot="milestones" , 
                # n_features = len(phloem), 
              style="italic",
              add_outline=True,
              basis="dendro",
              # feature_cmap = 'inferno', 
              show_segs = True, 
              fontsize = 10,
              figsize = (3,5),
                # annot = 'milestones',
              # heatmap_space=.2,
              # offset_names=.05,
              # linewidth_seg=2,
              # basis="dendro",
                # color="t",
               # cmap="viridis",
                ordering="max",
                ord_thre=0.99,
              # pseudo_cmap="viridis",
              save="phloem_selected.pdf")

# %% CORTEX trajectry plotting 
sc.set_figure_params(figsize=(2,10),dpi_save=500,frameon=False)
cortex.var.index = names_changes_list(cortex.var['gene_ids'])
scf.pl.trends(cortex,
                highlight_features = ['WUS',  
                                        'JKD', 
                                       'AED3', 
                                        'GAMT2',
                                        'EXT3', 
                                        'C/VIF2', 
                                        'PME3'
                                       ],
                # n_features = len(stem_cortex), 
              style="italic",annot="devtime",
              add_outline=True,
              basis="dendro",
              # feature_cmap = 'magma', 
                fontsize = 10,
              figsize = (3,5),
                 ordering="max",
                 ord_thre=0.99,
              pseudo_cmap="viridis",
              save="-OC-to_stem_cortex_selected.pdf")



# %% XYLEM
sc.set_figure_params(figsize=(2,10),dpi_save=600,frameon=False)
xylem.var.index = names_changes_list(xylem.var['gene_ids'])
scf.pl.trends(xylem,
                highlight_features = [  'WUS', 
                                        'CKX3', 
                                        'CLV1',
                                        'GH3.3',
                                        'TMO5',
                                        'ATHB-8',
                                        'PFA5',
                                        # 'PXY',
                                        'ACL5',
                                        'BAM3',
                                        'IAA16', 
                                        'SAC51'

                                        
                                      ],
                # n_features = len(xylem), 
              style="italic",annot="devtime",
              add_outline=True,
              basis="dendro",
              # color = 'leiden',
              # annot= 'milestones', 
              fontsize = 10,
              figsize = (3,5),
                ordering="max",
                ord_thre=0.99,
              pseudo_cmap="viridis",
               save="-OC_to_xylem_selected.pdf")

# %%  CAMBIUM plotting trajectory heatmap 
sc.set_figure_params(figsize=(2,10),dpi_save=600,frameon=False)
cambium.var.index = names_changes_list(cambium.var['gene_ids'])
scf.pl.trends(cambium,
                highlight_features = ['WUS', 
                                        'CLV1', 
                                        'GH3.3',
                                        'PFA5',
                                        'PXY',
                                        'ANAC020',
                                        'TMO5',
                                        'DOF2.2',
                                        'WAT1', 
                                        'CDKB2-1',
                                        'CYCA3-4', 
                                        'BUB3.1',

                                      ],
                # n_features = len(cambium), 
              style="italic",annot="devtime",
              add_outline=True,
              basis="dendro",
              # color = 'leiden',
              # annot= 'milestones', 
              fontsize = 10, 
              figsize = (3,5),
                # fig_heigth=3,
              # heatmap_space=.2,
              # feature_cmap = 'magma', 
              # linewidth_seg=2,
              # basis="dendro",
                # color="t",
              # cmap="viridis",
                ordering="max",
                ord_thre=0.99,
              pseudo_cmap="viridis",
              save="-OC_to_cambium_selected.pdf")

# %%  EP heatmap trajectory
sc.set_figure_params(figsize=(2,10),dpi_save=600,frameon=False)
EP.var.index = names_changes_list(EP.var['gene_ids'])
scf.pl.trends(EP,
                highlight_features = ['WUS',  
                                        'CKX3',
                                        'YAB1',
                                        'AHK4', 
                                        'PIN1', 
                                        'AHP6', 
                                        'IAA20',
                                        'REM1',
                                        'TPPJ',
                                        'GH3.3', 
                                        'H2B', 
                                        'CRF2', 
                                        # 'DRNL',
                                        ],
                # n_features = len(EP), 
              style="italic",
              annot="devtime",
              add_outline=True,
              basis="dendro",
              # color = 'leiden',
              # feature_cmap = 'magma', 
                fontsize = 10,
              figsize = (3,5),
                    ordering="max",
                    ord_thre=0.99,
              pseudo_cmap="viridis",
               save="-OC_to_EP.pdf"
              )


# %% XYLEM PARENCHYMA trajectory 
sc.set_figure_params(figsize=(2,10),dpi_save=600,frameon=False)
xylem_par.var.index = names_changes_list(xylem_par.var['gene_ids'])
scf.pl.trends(xylem_par,
                highlight_features = ['WUS', 
                                        'TMO5', 
                                           'GH3.6',
                                            'PER34',
                                            'PFA5',
                                            'SOBIR1',
                                            'ABCC14',
                                            'LOG5',
                                            'XTH4', 
                                            'FAF3', 
                                            'GH3.3',
                                            'SAC51', 
                                            'PXY',
                                            'BGLU46'
                                        ],
              style="italic",annot="devtime",
              add_outline=True,
              basis="dendro",
  
                fontsize = 10,
              figsize = (3,5),
            
                ordering="max",
                ord_thre=0.99,
              pseudo_cmap="viridis",
              save="-OC_to_vasc2_selected.pdf")

# %% pseudotime and supplementary figures FA plot
adata2 = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_internal_layers.h5ad')

sc.pl.draw_graph(adata2,color='t',color_map = 'viridis', add_outline=True, size =200, legend_fontsize=10, legend_fontoutline=2,
                  save="_pseudotimen.pdf")

sc.set_figure_params(figsize=(4,4),frameon=False, dpi=200)
sc.pl.draw_graph(adata2, color = ['AT2G45190'], layer = 'magic',  add_outline=True, legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'YAB1')
sc.pl.draw_graph(adata2, color = ['AT5G65140'], layer = 'magic', add_outline=True, legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'TPPJ')
sc.pl.draw_graph(adata2, color = ['AT3G25710'],layer = 'magic',  add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size = 300, 
                cmap = 'RdGy_r', title = 'TMO5')
sc.pl.draw_graph(adata2, color = ['AT5G19530'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'ACL5')
sc.pl.draw_graph(adata2, color = ['AT4G35190'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 300, 
                cmap = 'RdGy_r', title = 'LOG5')
sc.pl.draw_graph(adata2, color = ['AT1G76540'], layer = 'magic', add_outline=True,  legend_fontsize=10, legend_fontoutline=2, size= 200, 
                cmap = 'RdGy_r', title = 'CDKB2;1')



#%% Gene ontology analysis from SC to EP (first bifurcation)
import pandas as pd
from gprofiler import GProfiler
import numpy as np
trajectory = 'EP'
cc_var = pd.read_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/' + str(trajectory) + '_DEG_ordered_peak_expression.xlsx')

# Instantiate the object
gp = GProfiler(return_dataframe=True)

window_size = int(len(cc_var) * 0.25)
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
final_results = final_results[(final_results['p_value'] < 0.01) & (final_results['intersection_size'] > 5)
                              & (final_results['term_size'] <1000)]

## Removing not interesting columns 
final_results.drop(columns=['source', 'significant', 'evidences', 'parents','description'], inplace=True)


final_results.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_SC_to_EP_amp_rolling_' + str(window_size) + 'genes.xlsx')

# %%  Plotting all GO in one scatterplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Ensure intersection_size is numeric and drop NaNs
# final_results['intersection_size'] = pd.to_numeric(final_results['intersection_size'], errors='coerce')
# final_results = final_results.dropna(subset=['intersection_size', 'window_index', 'name', '-log10_pval'])

# Scaling factor for marker sizes
s_incr = 50
s_values = final_results['intersection_size'] * s_incr
s_values = pd.to_numeric(s_values)
# # Optional: clip extremely large marker sizes for readability
# s_values = np.clip(s_values, 20, 500)
 
# Create figure
fig, ax = plt.subplots(figsize=(3, 15), dpi=300)
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
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/SC_to_EP_amplyfging_GO.pdf',
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

final_results_names.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_SC_to_EP_amp_rolling_' + str(window_size) + 'genes_names.xlsx')

# %% Real plotting
# sc.set_figure_params(figsize = (10,3) , dpi_save=300)
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300  # set higher DPI globally
scf.pl.single_trend(EP,'WUS',ylab= 'WUS expression',layer ='magic',color_exp = 'grey', color_line='indigo',  basis = None, fitted_linewidth = 2, wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (10,3))
scf.pl.single_trend(EP,'CLV3',ylab= 'CLV3 expression',layer ='magic',color_exp = 'grey', color_line='indigo',  basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'CLV1',ylab= 'CLV1 expression',layer ='magic',color_exp = 'grey', color_line='indigo',  basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))

scf.pl.single_trend(EP,'YAB1',ylab= 'YAB1 expression',layer ='magic',color_exp = 'grey', color_line='indigo',  basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'GH3.2',ylab= 'GH3.2 expression',layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'GH3.3',ylab= 'GH3.3 expression',layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'IAA30',ylab= 'IAA30 expression',layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'DRNL',ylab= 'DRNL expression',layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))

# %% hypoxia-response from SC to EP
scf.pl.single_trend(EP,'TIL',ylab= 'TIL expression', layer ='magic',color_exp='grey', color_line='indigo',wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (10,3))
scf.pl.single_trend(EP,'SUS1',ylab= 'SUS1 expression',layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))

scf.pl.single_trend(EP,'AT3G10040',ylab= 'HRA1 expression',plot_emb = False, layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'AT4G27450',ylab= 'HUP54 expression',plot_emb = False,layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))



scf.pl.single_trend(EP,'RPP0B',ylab= 'RPP0B expression',plot_emb = False,layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.050, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))
scf.pl.single_trend(EP,'RPP2B',ylab= 'RPP2B expression',plot_emb = False,layer ='magic',color_exp='grey', color_line='indigo', basis = 'draw_graph_fa', wspace=-.50, size_expr = 25, alpha_expr = 0.1,  figsize = (3,3))























