#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 17:06:16 2023

@author: Sebastian
"""

from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import palantir
import pandas as pd
import numpy as np
import warnings
import scanpy as sc
import scFates as scf
import os
import sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"


warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
# import scFates


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
    adata.obsm["X_pca2d"] = adata.obsm["X_pca"][:, 1::]
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
    # Normalize pseudotime for trajectory
    ###################################
    # from scipy.stats import rankdata
    # Step 1: Get pseudotime
    # pt = adata.obs['t'].values
    # # Step 2: Rank and rescale to [0, 1]
    # pt_smooth = rankdata(pt, method='average') / len(pt)
    # # Step 3: Store in AnnData
    # adata.obs['t'] = pt_smooth
    # sc.pl.draw_graph(adata2,color='pseudotime_smooth',color_map = 'viridis', add_outline=True, size =200, legend_fontsize=10, legend_fontoutline=2,
    #                  save="_pseudotimen.pdf")
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

# %% Cell-cycle dynamic

# Data with previus information about HVG and UMAP
adata = sc.read_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_filtered_norm_clusters.h5ad')

sc.set_figure_params(figsize=(4, 4), frameon=False)
adata1 = scfates_trajectories_alignment(adata, clusters=['S-phase', 'G2/M-phase'],
                                        root_gene='AT4G33270',
                                        nodes=100,
                                        sigma=0.002,
                                        eigen=3)


# %% trajectory milestones
adata2 = scfates_trajectories_dendogram(
    adata1,  root_segment=87)  # assign root segment/tip cell
# save with new cluster names
adata2.write_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle.h5ad')

# %%  Selecting main cell cycle branch
adata2 = sc.read_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle.h5ad')

adata2 = scf.tl.subset_tree(adata2, root_milestone='87', milestones=[
                            '43'], copy=True)  # Focus only on main trajectory

scf.tl.tree(adata2, method="ppt", Nodes=200, use_rep="palantir", plot=False,
            device="cpu", seed=1,  ppt_lambda=100, ppt_nsteps=200,
            ppt_sigma=0.01)
sc.set_figure_params(figsize=(4, 4), frameon=False, dpi = 500)
scf.pl.graph(adata2, color_cells='leiden',
             palette=sns.color_palette('colorblind'))

# assign root segment/tip cell to the now main cell cycle trajectory
adata2 = scfates_trajectories_dendogram(adata2,  root_segment=95)
# save with new cluster names
adata2.write_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_subtree.h5ad')

# %% DEG along cell cycle trajectory

adata2 = sc.read_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_subtree.h5ad')

## Smoothing pseudotime within trajectory
from scipy.stats import rankdata
pt = adata2.obs['t'].values
pt_smooth = rankdata(pt, method='average') / len(pt)
adata2.obs['t'] = pt_smooth

scf.tl.test_association(adata2, n_jobs=1)
adata2.write_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_ass.h5ad')

# %%  DEG along cell cycle trajectory  2
adata2 = sc.read_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_ass.h5ad')
# Use this filter instead the default one
adata2.var['signi'] = adata2.var['p_val'] < 0.01
scf.tl.fit(adata2, n_jobs=1)
adata2.write_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_ass_fitted.h5ad')

# %% Plotting some DEG along trajectory graph
adata2 = sc.read_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_ass_fitted.h5ad')
sc.set_figure_params(figsize=(4, 4), frameon=False, dpi=500)
sc.pl.draw_graph(adata2, color=['AT5G10390'],  add_outline=True, layer = 'magic',  legend_fontsize=10, legend_fontoutline=2,  title='HTR13',
                 cmap='RdGy_r', save='HTR13.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT4G33270'],  add_outline=True, layer = 'magic', legend_fontsize=10, legend_fontoutline=2,  title='CDC20-1',
                 cmap='RdGy_r', save='CDC20-1.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT5G11510'],  add_outline=True,layer = 'magic',  legend_fontsize=10, legend_fontoutline=2,  title='MYB3R4',
                 cmap='RdGy_r', save='MYB3R4.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT2G28740'],  add_outline=True,layer = 'magic',  legend_fontsize=10, legend_fontoutline=2,  title='HIS4',
                 cmap='RdGy_r', save='HIS4.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT3G25980'],  add_outline=True, layer = 'magic', legend_fontsize=10, legend_fontoutline=2,  title='MAD2',
                 cmap='RdGy_r', save='MAD2.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT1G08880'],  add_outline=True, layer = 'magic', legend_fontsize=10, legend_fontoutline=2,  title='HTA5',
                 cmap='RdGy_r', save='HTA5.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT3G20670'],  add_outline=True,layer = 'magic',  legend_fontsize=10, legend_fontoutline=2,  title='HTA13',
                 cmap='RdGy_r', save='HTA13.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT2G38810'],  add_outline=True,layer = 'magic',  legend_fontsize=10, legend_fontoutline=2,  title='HTA8',
                 cmap='RdGy_r', save='HTA8.pdf', size=300)
sc.pl.draw_graph(adata2, color=['AT1G63100'],  add_outline=True, layer = 'magic', legend_fontsize=10, legend_fontoutline=2,  title='SCL28',
                 cmap='RdGy_r', save='SCL28.pdf', size=300)


sc.set_figure_params(figsize=(4, 4), frameon=False, dpi=500)
sc.pl.draw_graph(adata2, color='t', color_map='viridis', add_outline=True, size=200, legend_fontsize=10, legend_fontoutline=2,
                 save="_pseudotimen.pdf")

sc.set_figure_params(figsize=(4, 4), frameon=False, dpi=500)
sc.pl.draw_graph(adata2, color='leiden', palette=( 'khaki','lightcoral'), add_outline=True,
                 size=200, legend_fontsize=15, legend_fontoutline=2, frameon=True,
                 save='cc_leiden.pdf')

# %% Generate datasheet with cell cycle DEG  and main heatmap figure
adata2 = sc.read_h5ad(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_scfates_cellcycle_ass_fitted.h5ad')
adata2.var["corr"] = list(map(lambda g: pearsonr(
    adata2.obs.t, adata2[:, g].layers["fitted"].flatten())[0], adata2.var_names))
adata2.var["up"] = adata2.var["corr"] > 0
fitted = pd.DataFrame(adata2[:, adata2.var_names].layers["fitted"],
                      index=adata2.obs_names, columns=adata2.var_names).T.copy(deep=True)
feature_order = (fitted.apply(lambda x: adata2.obs.t[fitted.columns][(
    x - x.min()) / (x.max() - x.min()) > 0.7].mean(), axis=1, ).sort_values().index)
cc_var = adata2.var.loc[feature_order]
cc_var['gene_name'] = names_changes_list(cc_var.gene_ids)
cc_var.to_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/cell_cycle_DEG_ordered_peak_expression.xlsx')

adata2.var.index = adata2.var['gene_ids']
adata2.var_names = names_changes_list(adata2.var_names)


# %%
scf.pl.trends(adata2,
              style="italic", annot="devtime",
              fontsize=10,
              highlight_features=['HIS4',
                                    'HTA3',
                                  'HTR2',
                                    'CYCB1-1',
                                  'CDKB2-1',
                                  'KIN7N',
                                  'CDKB2-1',
                                  'MYB3R4',
                                  'CYCA1-1',
                                  'SMR6',
                                  'FZR3',
                                  'IMPA-6',
                                  'CDC20-2',
                                  'KNOLLE',
                                  'SCL28',
                                  'CYCB1-2',
                                  'PCNA',
                                  'CSLD5',
                                  'CDKB2-2',
                                  'MED36A'
                                  ],
              # n_features= len(adata2),
              layer = 'magic', 
              figsize=(5, 5),
              add_outline=True,
              s=30,
              heatmap_space=.2,
              offset_names=.05,
              basis="draw_graph_fa",
              color="t",
              cmap="viridis",
              ord_thre=0.9,
              ordering="max",
              pseudo_cmap="viridis",
              save="-cell_cycle_amplitude.pdf")


# %% Figure 5e, gene expression profiles for some genes related to cell cycle 
scf.pl.single_trend(adata2,"HIS4",basis = 'draw_graph_fa', layer = 'magic',  wspace=-.025,save="HIS4_cell_cycle.pdf", size_expr = 50, alpha_expr = 0.02, ylab = 'HIS4 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"KNOLLE",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025,save="KNOLLE_cell_cycle.pdf", size_expr = 50, alpha_expr = 0.02, ylab = 'KNOLLE expression', figsize = (8,3))
scf.pl.single_trend(adata2,"HTA5",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025,save="HTA5_cell_cycle.pdf", size_expr = 50, alpha_expr = 0.02, ylab = 'HTA5 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"CDC20-1",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025,save="CDC20-1_cell_cycle.pdf", size_expr = 50, alpha_expr = 0.02, ylab = 'CDC20-1 expression', figsize = (8,3))


scf.pl.single_trend(adata2,"CDKB2-1",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025, size_expr = 50, alpha_expr = 0.02, ylab = 'CDKB2-1 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"HTA13",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025, size_expr = 50, alpha_expr = 0.02, ylab = 'HTA13 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"CYCB1-1",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025, size_expr = 50, alpha_expr = 0.02, ylab = 'CYCB1;1 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"MAD2",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025, size_expr = 50, alpha_expr = 0.02, ylab = 'MAD2 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"SCL28",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025, size_expr = 50, alpha_expr = 0.02, ylab = 'SCL28 expression', figsize = (8,3))
scf.pl.single_trend(adata2,"FZR3",basis = 'draw_graph_fa',layer = 'magic',wspace=-.025, size_expr = 50, alpha_expr = 0.02, ylab = 'FZR3 expression', figsize = (8,3))


#%% Gene ontology analysis over cell cycle trajectory (dataframe)
import pandas as pd
from gprofiler import GProfiler
import numpy as np

cc_var = pd.read_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/cell_cycle_DEG_ordered_peak_expression.xlsx')

# Instantiate the object
gp = GProfiler(return_dataframe=True)

window_size = int(len(cc_var) * 0.15)
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
final_results = final_results[(final_results['p_value'] < 0.01) & (final_results['intersection_size'] > 5)
                              & (final_results['term_size'] <1000)]

final_results.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_cell_cycle_rolling_' + str(window_size) + 'genes.xlsx')


# %%  Plotting all GO in one scatterplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Ensure intersection_size is numeric and drop NaNs
final_results['intersection_size'] = pd.to_numeric(final_results['intersection_size'], errors='coerce')
final_results = final_results.dropna(subset=['intersection_size', 'window_index', 'name', '-log10_pval'])

# Scaling factor for marker sizes
s_incr = 25
s_values = final_results['intersection_size'] * s_incr

# # Optional: clip extremely large marker sizes for readability
# s_values = np.clip(s_values, 20, 500)
 
# Create figure
fig, ax = plt.subplots(figsize=(10, 15), dpi=500)
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
    alpha=0.7,
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

ax.set_xlabel('Cell cycle progression (S-G2-M)', fontsize=14)
ax.set_ylabel('GO term', fontsize=14)
ax.set_title('GO enrichment across windows', fontsize=16)

# Save and show
plt.tight_layout()
plt.savefig(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/cell_cycle_GO.pdf',
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

final_results_names.to_excel('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/GO_cell_cycle_rolling_' + str(window_size) + 'genes_names.xlsx')


# %%  Chip-Seq analysis for floral homeotic genes (Figure 5H)
chip_seq = pd.read_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/chip_seq_flowers.xlsx', header=1)
chip_seq['Name'] = names_changes_list(chip_seq['Gene ID'])


N = 50   ### To use the same number of genes per cluster, using same number of top marker genes per cluster 
genes = pd.read_excel(
    "/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx")
# genes = genes.sort_values(by = ['pvals_adj'], ascending = True)
s = genes[genes['group'] == 'S-phase'][0:N]
m = genes[genes['group'] == 'G2/M-phase'][0:N]

cc_df = pd.DataFrame({
    'S-phase': s['names'].reset_index(drop=True),
    'G2/M-phase': m['names'].reset_index(drop=True)
})

tf = ['AP1',
      'PI',
      'LFY',
      'SEP3',
      'AP3',
      'SVP',
      'ETT',
      'JAG',
      'SOC1',
      'BLR',
      'FLC',
      'AP2',
      'RGA',
      'FLM',
      'AG',
      # 'SUP'
      ]

tf_targets = pd.DataFrame()
# chip_seq[chip_seq['Regulators'].str.split(',', expand = True).columns.tolist()] = chip_seq['Regulators'].str.split(',', expand = True)
for f in tf:
    temp1 = chip_seq[chip_seq['Regulators'].str.contains(f)]
    temp1 = temp1[temp1['Gene ID'].str.contains('AT')]
    # NaN_needed = 9000 - len(temp1)
    temp2 = pd.DataFrame(temp1['Gene ID'].tolist())
    temp2.columns = [f]
    tf_targets = pd.concat([tf_targets, temp2], axis=1)


intersection_df = pd.DataFrame()
gene_list = pd.DataFrame()
for col in tf_targets.columns:
    for phase in ['S-phase', 'G2/M-phase']:
        intersection_result = list(set(tf_targets[col]) & set(cc_df[phase]))
        d = {'TF': col,
             'chip_seq_target': len([x for x in tf_targets[col] if str(x) != 'nan']),
             'phase_genes': len([x for x in cc_df[phase] if str(x) != 'nan']),
             'intersection': len(intersection_result),
             'phase': str(phase),
             'genes': [','.join([str(gene) for gene in names_changes_list(intersection_result) if pd.notna(gene)])]
             }
        df = pd.DataFrame(d, index=[0])
        intersection_df = pd.concat([intersection_df, df])

intersection_df = intersection_df.reset_index(drop=True)

intersection_df['porcentage_phase'] = intersection_df['intersection'] / \
    intersection_df['phase_genes'] * 100
intersection_df['porcentage_chip_seq'] = intersection_df['intersection'] / \
    intersection_df['chip_seq_target'] * 100


intersection_df.to_excel(
    '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/intersection_chip_seq_cell_cycle.xlsx')

plt.figure(figsize=(4 , 6), dpi = 300)
sns.set(style="ticks")
plt.rcParams['lines.color'] = 'black'
ax = sns.barplot(data=intersection_df,  x="TF", y="porcentage_phase",    hue='phase',
                 palette=['lightcoral','khaki' ],  linewidth=1, width=0.8,
                 edgecolor='black',  errwidth=0.5, errcolor='b'
                 )
ax.grid(False)
plt.setp(ax.artists, edgecolor='k', facecolor='w')
plt.setp(ax.lines, color='k')
sns.despine(offset=3, trim=False)
ax.legend(bbox_to_anchor=(1, 1.05))
plt.setp(ax.artists, edgecolor='k', facecolor='w')
ax.tick_params(axis='x', rotation=90)
ax.set(xlabel='Floral TF', ylabel='% target genes per cluster')
plt.setp(ax.lines, color='k')
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/floral_genes_porcentage_s_phase' + ".pdf",  bbox_inches='tight')
plt.show()


# # %%  Cell cycle identity per cluster based on gene expression (Figure 5F) 
# genes = pd.read_excel(
#     "/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx")
# # Data with previus information about HVG and UMAP
# adata = sc.read_h5ad(
#     '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_filtered_norm_clusters.h5ad')

# df_effect_cc = pd.DataFrame()
# for cell_cycle in ['S-phase', 'G2/M-phase']:
#     phase_marker = genes[genes['group'] == cell_cycle][0:50]
#     for cluster in adata.obs['leiden'].unique().tolist():
#         if cluster not in ['S-phase', 'G2/M-phase']:
#             temp3 = adata[adata.obs['leiden'] == cluster]
#             df = temp3.to_df()
#             df2 = df.loc[:, df.columns.isin(phase_marker["names"])]
#             mean_expr_per_cell = pd.DataFrame(df2.mean(axis=1))
#             mean_expr_per_cell['cluster'] = cluster
#             mean_expr_per_cell['cell_cycle'] = cell_cycle
#             mean_expr_per_cell = mean_expr_per_cell.reset_index(drop=True)
#             df_effect_cc = pd.concat([df_effect_cc, mean_expr_per_cell])


# df_effect_cc.columns = ['gene_expression', 'cluster' ,'cell_cycle']


# # Plotting S-phase only
# sphase = df_effect_cc[df_effect_cc['cell_cycle'] == 'S-phase']

# plt.figure(figsize=(4, 3))
# sns.set(style="ticks")
# ax = sns.barplot(data=sphase,  x="cluster", y="gene_expression", # hue='cell_cycle',
#                  # hue_order=['S-phase', 'G2/M-phase'],
#                  color='lightcoral', linewidth=1.5,
#                  # order = cluster_order,
#                  edgecolor='black', capsize=0.2, errwidth=2, errcolor='black'
#                  )
# ax.grid(False)
# plt.setp(ax.artists, edgecolor='k', facecolor='w')
# plt.setp(ax.lines, color='k')
# sns.despine(offset=3, trim=False)

# plt.setp(ax.artists, edgecolor='k', facecolor='w')
# ax.tick_params(axis='x', rotation=90)
# ax.set(xlabel='Cell domains', ylabel='Mean gene expression S-phase marker genes')
# plt.setp(ax.lines, color='k')
# plt.legend(bbox_to_anchor=(1, 1.3), ncol=10)
# sns.despine(offset=3, trim=False)
# plt.savefig("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/S-phase_gene_gene_expression_clusters.pdf",  bbox_inches='tight')
# plt.show()


# # Plotting M-phase only
# mphase = df_effect_cc[df_effect_cc['cell_cycle'] == 'G2/M-phase']

# plt.figure(figsize=(4, 3))
# sns.set(style="ticks")
# ax = sns.barplot(data=mphase,  x="cluster", y="gene_expression", # hue='cell_cycle',
#                  # hue_order=['S-phase', 'G2/M-phase'],
#                  color='khaki', linewidth=1.5,
#                  # order = cluster_order,
#                  edgecolor='black', capsize=0.2, errwidth=2, errcolor='black'
#                  )
# ax.grid(False)
# plt.setp(ax.artists, edgecolor='k', facecolor='w')
# plt.setp(ax.lines, color='k')
# sns.despine(offset=3, trim=False)

# plt.setp(ax.artists, edgecolor='k', facecolor='w')
# ax.tick_params(axis='x', rotation=90)
# ax.set(xlabel='Cell domains', ylabel='Mean gene expression G2/M-phase marker genes')
# plt.setp(ax.lines, color='k')
# plt.legend(bbox_to_anchor=(1, 1.3), ncol=10)
# sns.despine(offset=3, trim=False)
# plt.savefig("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/M-phase_gene_gene_expression_clusters.pdf",  bbox_inches='tight')
# plt.show()



# ## S-phase statistic 
# from statsmodels.stats.multicomp import pairwise_tukeyhsd

# tukey_s = pairwise_tukeyhsd(
#     endog=sphase["gene_expression"],
#     groups=sphase["cluster"],
#     alpha=0.01
# )
# tukey_df_s = pd.DataFrame(
#     tukey_s.summary().data[1:],
#     columns=tukey_s.summary().data[0]
# )

# ## G2/M-phase statistic 
# from statsmodels.stats.multicomp import pairwise_tukeyhsd

# tukey_m = pairwise_tukeyhsd(
#     endog=mphase["gene_expression"],
#     groups=mphase["cluster"],
#     alpha=0.01
# )
# tukey_df_m = pd.DataFrame(
#     tukey_m.summary().data[1:],
#     columns=tukey_m.summary().data[0]
# )

# ## Save dataframes 

# tukey_df_s.to_excel(
#     '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/tukey_df_s_phase.xlsx')

# tukey_df_m.to_excel(
#     '/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/tukey_df_m_phase.xlsx')


