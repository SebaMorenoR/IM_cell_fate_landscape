
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 13:26:39 2023

@author: Sebastian
"""


import pandas as pd
import matplotlib.cm
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
import os
import scanpy.external as sce

import matplotlib as mpl
mpl.rcParams['font.family'] = ['Arial']
plt.rcParams['axes.grid'] = False

from sklearn.linear_model import LinearRegression
from scipy import stats
from statsmodels import robust
import scrublet as scr
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


os.chdir('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/')
path = "/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/"

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


# %% Importing datasets 
#  Using CELLBENDER to remove environmental reads.
data1 = sc.read_10x_h5('data/output_cellbender_filtered_resequenced_SITTC6.h5')  # submitted ##
data1.obs['Sample'] = 'sittc6'

data3 = sc.read_10x_h5('data/output_cellbender_filtered_resequenced_SITTE2.h5')
data3.obs['Sample'] = 'sitte2'

data4 = sc.read_10x_h5('data/output_cellbender_filtered_resequenced_SITTA1.h5')
data4.obs['Sample'] = 'sitta1'

data1.var_names_make_unique()
data3.var_names_make_unique()
data4.var_names_make_unique()

# %% Combine samples into one anndata
adata = data1.concatenate(data3.concatenate(data4))

# Add GENE NAMES Only if I need to remove Ribosome genes 
# new_list = names_changes_list(adata.var.index)
# adata.var['gene_ids'] = new_list

## DOUBLET DETECTION
scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
scrub.call_doublets(threshold=0.25)
scrub.plot_histogram()
adata.obs["Doublet"] = predicted_doublets
adata = adata[adata.obs.Doublet == False, :] ## Removing doublets barcodes 


adata.var["atm"] = adata.var_names.str.startswith("ATM")
adata.var["atc"] = adata.var_names.str.startswith("ATC")

sc.pp.calculate_qc_metrics(adata,qc_vars=["atm", 'atc'], percent_top = None, log1p=False, inplace=True)

ax = sc.pl.scatter(adata,x = "total_counts", y = "pct_counts_atc")
ax = sc.pl.scatter(adata,x = "total_counts", y = "pct_counts_atm")
ax = sc.pl.scatter(adata,x = "total_counts", y = "n_genes_by_counts")
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts" ], 
             jitter=  4, multi_panel=True)
sc.pl.violin(adata, ["pct_counts_atm", 'pct_counts_atc'], stripplot=False,
             ylabel='Percentage of UMI counts %') 
             # save='mito_chloro_count.pdf')

adata2 = adata

# sc.pp.filter_genes(adata2, max_counts=np.quantile(adata2.var['total_counts'], .995))

## Removing low quality cells  by counts 
gene_counts = adata2.obs['n_genes_by_counts']
median = np.median(gene_counts)
mad = robust.mad(gene_counts)
# 5x MAD range
lower = median - 5 * mad
upper = median + 5 * mad
adata2 = adata2[(gene_counts >= lower) & (gene_counts <= upper)].copy()


## Removing low quality cells by %  mitochondrial genes
adata2 = adata2[adata2.obs.pct_counts_atm < 10, :] ## Mitochondria 
sc.pp.filter_genes(adata2, min_cells= 10)

# %%  Removing mitochondrial and chloroplast genes
mito_genes = adata2.var[adata2.var['atm'] == True].index.tolist()
cloro_genes = adata2.var[adata2.var['atc'] == True].index.tolist()


adata2 = adata2[:, ~adata2.var_names.isin(mito_genes)]
adata2 = adata2[:, ~adata2.var_names.isin(cloro_genes)]

adata2.var_names_make_unique()


adata2.layers['counts'] = adata2.X.copy()
adata2.raw = adata2
adata2.write_h5ad(path + '/adata_filtered.h5ad')

sc.pp.normalize_total(adata2,
                      target_sum=1e4,
                      exclude_highly_expressed=False)  # which normalization value to use to scale reads?


sc.pp.log1p(adata2)
adata2.raw = adata2
adata2.write_h5ad(path + '/adata_norm.h5ad')


# %%  REMOVING OUTLIERS using bulk-RNA-seq 
norm_counts = pd.read_csv('/Users/Sebastian/Documentos/SLCU_lab/results/hormo_RNA_seq/norm_counts.csv', index_col=0)
norm_counts = norm_counts[['mock4h', 'mock4h.1', 'mock4h.2']]
norm_counts = np.log2(norm_counts + 1)
norm_counts['mean'] = norm_counts.mean(axis=1)

# normalized counts single-nuclei
sc_adata = sc.read_h5ad(path + '/adata_norm.h5ad')
sc_adata_df = sc_adata.to_df()

shared_genes = list(set(norm_counts.index) & set(sc_adata_df.columns))

rna_seq = norm_counts[norm_counts.index.isin(shared_genes)]
rna_seq = rna_seq['mean']
rna_seq = rna_seq.sort_index()
rna_seq = pd.DataFrame(rna_seq)
rna_seq.columns = ['rna_seq']

sc = np.log2(sc_adata_df.sum(axis=0) + 1)
sc = sc[sc.index.isin(shared_genes)]
sc = sc.sort_index()

rna_seq['sc'] = sc.tolist()

fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi = 500)
sns.scatterplot(
    data=rna_seq,
    x="sc",
    y="rna_seq",
    color="purple",
    alpha=0.3,
    ax=ax,
)

plt.ylabel('Inflorescence meristem bulk RNA-seq (log2)', size = 14)
plt.xlabel('Inflorescence meristem pseudo-bulk snRNA-seq (log2)', size = 15)
plt.show()


res = stats.spearmanr(rna_seq['sc'], rna_seq['rna_seq'])
print(res)

# Perform linear regression
X = rna_seq[['rna_seq']]  # Independent variable
Y = rna_seq['sc']  # Dependent variable
model = LinearRegression()
model.fit(X, Y)
rna_seq['sc_pred'] = model.predict(X)
# Compute residuals
rna_seq['Residuals'] = rna_seq['sc'] - rna_seq['sc_pred']
threshold = 5  # Can be adjusted to 2.5 or 3 for stricter outlier detection
rna_seq['Outlier'] = np.abs(rna_seq['Residuals']) > threshold
# sns.histplot(rna_seq['Residuals'], kde=True)

# Extract outlier data points
outliers = rna_seq[rna_seq['Outlier']]
outliers['name'] = names_changes_list(outliers.index)

# %% GO for outliers removed in this pipeline
from gprofiler import GProfiler


outliers = outliers[outliers['Residuals'] < 0]
# Instantiate the object
gp = GProfiler(return_dataframe=True)

# Extract gene IDs (TAIR IDs)
gene_list = outliers.index.tolist()

# Run GO enrichment for Arabidopsis
results = gp.profile(
    organism='athaliana',
    query=gene_list,
    significance_threshold_method='fdr', 
    sources=['GO:BP'], 
    user_threshold=0.05, 
    no_evidences=False
)

final_results = results
final_results['-log10_pval'] = -np.log10(final_results['p_value'])
## Incresing threhold? 
final_results = final_results[(final_results['p_value'] < 0.01) & (final_results['intersection_size'] > 5)
                              & (final_results['term_size'] <1000)]

## Removing not interesting columns 
final_results.drop(columns=['source', 'significant', 'evidences', 'parents','description'], inplace=True)




# %% # removing outliers outside linear regression comparing RNA-seq with snRNA-seq
adata2 = adata2[:, ~adata2.var_names.isin(outliers.index)]
adata2.var_names_make_unique()

adata2.write_h5ad(path + '/adata_filtered_norm.h5ad')

# %%  HVG
import scanpy as sc

sc.pp.highly_variable_genes(adata2,  n_top_genes=1000, subset=False,
                            layer='counts', flavor='seurat_v3', batch_key='Sample')  # default values

# %%  Calculate the percentage of total counts contributed by each gene
total_counts_sum = adata2.var['total_counts'].sum()
sorted_counts = adata2.var['total_counts'].sort_values(ascending=False)
cumulative_perc = np.cumsum(sorted_counts) / total_counts_sum * 100

# Plot it
plt.plot(range(100), cumulative_perc[:100])
# plt.axhline(y=10, color='r', linestyle='--') # Mark 10% of total library
plt.ylabel('Cumulative % of Total Reads')
plt.xlabel('Gene Rank (Top 100)')
plt.show()
# %%  Identify the top  globally expressed genes from the main object
top_genes_count = 20
top_genes_removing = adata2.var['total_counts'].sort_values(ascending=False).head(top_genes_count).index.tolist()
manual_exclude = (
    adata2.var_names.isin(top_genes_removing) 
)
adata2.var.loc[manual_exclude, 'highly_variable'] = False
pcHVG_data = adata2[:, adata2.var.highly_variable].copy()

# %%  SCALE DATA
sc.pp.scale(pcHVG_data, max_value=10)


# REGRESSING OUT SOURCE VARIATION
sc.pp.regress_out(pcHVG_data, [ "total_counts", "pct_counts_atm"])
# 5.- DATA VISUALIZATION
# COMPUTING PCA
sc.tl.pca(pcHVG_data)  # COMPUTING PCA to reduce dimension of the data
# sc.pl.pca_variance_ratio(pcHVG_data, n_pcs = 50) ## To decide how many PCA compute

# UMAP

# sc.external.pp.bbknn(pcHVG_data, batch_key='Sample', n_pcs= 50) ### integrate different batches
sce.pp.harmony_integrate(pcHVG_data, key="Sample")
pcHVG_data.obsm['X_pca'] = pcHVG_data.obsm['X_pca_harmony']

sc.pp.neighbors(pcHVG_data, n_neighbors=50)

#%% Add GENE NAMES 
new_list = names_changes_list(pcHVG_data.raw.var.index)
pcHVG_data.raw.var['gene_ids'] = new_list

#%% UMAP leiden pipeline 
for dist in [0.1]:
    for sp in [3]:
        sc.tl.umap(pcHVG_data,
                   min_dist=dist,
                   spread=sp)
        r = 1.3
        sc.set_figure_params(figsize=(7,5),frameon=True, dpi=50)
        sc.tl.leiden(pcHVG_data, resolution=r)
        sc.pl.umap(pcHVG_data, color='leiden',  # legend_loc='on data',
                   add_outline=False, alpha=0.8, size=50, 
                   #title= 'clv3-15 UMAP: n_neighbors' + str(nn) + '_PCA:' + str(pcs) +  'min_dist:' + str(dist) + '_spread:' + str(sp),
                   frameon=False, legend_fontsize=7, legend_fontoutline=2,
                   palette='tab20', 
                   )    
        iSC = [ 'WUS','CLV3', 'CLV1', 'AT1G26680','AT3G59270', 'CKX3', 'AHK4', 
               'TMO5', 'SMXL5','PFA5', 
               # 'TGG2', 'TGG1',
               # 'JKD','AED3', 
               # 'HIS4', 'HTB1', 'HTB11','HTA5', 
               # 'PI', 'AP3', 'AG', 'SEP1',
               'CUC2', 'CUC3',  'SOD7', 
               'ATML1', 'KAN1', 'MYB4', 
               'AHP6', 'GH3.3','DOF5.8',
               'FAF1', 'FAF3', 'GH3.6', 
               ] 
        sc.pl.matrixplot(pcHVG_data, var_names = iSC , groupby= 'leiden', 
                         swap_axes = True , cmap = 'Purples', 
                        standard_scale='var', 
                        gene_symbols= 'gene_ids', 
                        figsize= (8,6),
                        colorbar_title ='Standarized  \nmean expression  \nper cluster',
                        linewidth = 0.5, color = 'black')

   
# %% Plotting CLV3 and WUS marker genes 
sc.set_figure_params(figsize=(5,4),frameon=False, dpi=150)
sc.pl.umap(pcHVG_data, color=['AT2G17950', 'AT2G27250', 'AT1G26680', 'AT3G59270'], title = ['WUS', 'CLV3', 'AT1G26680', 'CKX3'],
            add_outline=False,alpha =1, size =  200,
        frameon=True, legend_fontsize=7, legend_fontoutline=5,
        cmap="magma",)

# %%  MARKER GENES 
sc.tl.rank_genes_groups(pcHVG_data, groupby='leiden', method='wilcoxon',
                          corr_method='benjamini-hochberg', tie_correct=True, pts=True)
markers = sc.get.rank_genes_groups_df(pcHVG_data, group=None)
markers["min_diff_pct"] = (markers["pct_nz_group"] - markers["pct_nz_reference"])
markers["pct_ratio"] = ( (markers["pct_nz_group"] ) /  (markers["pct_nz_reference"] ))
markers_def = markers[(markers.pvals_adj < 0.01) &
                      (markers.logfoldchanges > 1) 
                        &  (markers.pct_nz_group > 0.05)
                      ]


markers_def['names_ref'] = names_changes_list(markers_def['names'])
markers_def.to_excel(
    "/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg.xlsx")


# %% # save with without cluster names
pcHVG_data.write_h5ad(path + '/pcHVG_data_proccessed.h5ad')

#%%  RESUME FROM HERE TO SAVE ANNDATA WITH CLUSTER NAMES 
pcHVG_data = sc.read_h5ad(path + '/pcHVG_data_proccessed.h5ad')
# Re-naming clusters
rename_dict = {  # cluster should be in the orther from 0 to final cluster
    '0': 'photosynthetic cells',
    '1': 'cortex',
    '2':  'S-phase', 
    '3': 'G2/M-phase' , 
    '4': 'procambium', 
    '5': 'unknown 1',
    '6':  'rib zone' , 
    '7': 'stem epidermis',   
    '8':  'IM epidermis',
    '9': 'unknown 2' , 
    '10': 'xylem parenchyma' , 
    '11':  'EP',
    '12': 'floral identity',
    '13': 'undifferentiated cells', 
    '14':  'pollen' , 
    '15': 'BD', 
    '16': 'cell wall-related',
    '17': 'stress-response', 
}

pcHVG_data.obs['leiden'] = pcHVG_data.obs['leiden'].map(rename_dict)

####################
# Making dotplot with specific marker genes
########################
# Change names to more common ones

markers = {
    'stress-response': ['ERF109', 'WRKY40', 'LOX3'],
    'cell wall-related': ['PME5', 'PME4', 'AGP23'],
    'BD': ['CUC2', 'CUC3', 'SOD7'],  # CUC3, CUC2, WOX12
    'pollen': ['AMS', 'ABCG26', 'CYP704B1'],  # AMS, ABCG26
    'undifferentiated cells': ['WUS', 'CKX3', 'CLV3'],  # AT3G59270,  WUS,
    'floral identity': ['AP3', 'PI', 'AG'],
    'EP': ['AHP6', 'DRNL', 'DOF5.8'],  # AHP6, FIL,
    'xylem parenchyma': ['WAK1', 'GH3.6', 'FAF1'],
    'IM epidermis': ['ATML1', 'MYB4', 'HTH'],  # 'PDF1', 'ATML1'
    'stem epidermis': ['CER1', 'MLP28', 'CUT1'],
    'rib zone': [ 'PCO2','AT1G29270',  'ARD',],
    # 'hypoxia-response': [ 'MIOX2', 'GDPD2'],
    'procambium': ['PFA3', 'ATHB-8', 'PFA5'],
    'G2/M-phase': ['CDC20-1', 'CDC20-2', 'FZR3'],
    'S-phase': ['HTA6', 'HTA7', 'HTA10'],
    'cortex': ['JKD',  'DIR13','AED3',],  # JKD, pCORTEX
    'photosynthetic cells': ['LHCB1.1', 'LHCB2.2', 'GER3'],
        }
markers = dict(reversed(markers.items()))

orig_map = plt.cm.get_cmap('Purples')

# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()


sc.set_figure_params(frameon=False, dpi=100)
sc.pl.dotplot(adata=pcHVG_data, var_names=markers, groupby='leiden',
              dendrogram=False,
               standard_scale='var',
              dot_max=0.4,
              gene_symbols='gene_ids',
              figsize=(21, 4.5),
              cmap='BuPu',
              mean_only_expressed=False,        
              save='dotplot_marker_genes.pdf'
              )


# save with new cluster names
pcHVG_data.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/pcHVG_data_proccessed_names.h5ad')


#%% GENERATE DIFFERENTIALLY EXPRESSED (MARKER GENES) TABLE
pcHVG_data = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/pcHVG_data_proccessed_names.h5ad')
pcHVG_data.uns['log1p']["base"] = None

markers = sc.get.rank_genes_groups_df(pcHVG_data, group=None)
markers['group'] = markers['group'].map(rename_dict)
markers_def = markers[(markers.pvals_adj < 0.01) & (markers.logfoldchanges > 1) & (markers.pct_nz_group > 0.05)]
markers_def['names_ref'] = names_changes_list(markers_def['names'])
markers_def.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final.xlsx")

# %% Final UMAP 
# STEP 1: Get the original tab20 colors
tab20_palette = sns.color_palette("tab20", 20)
# STEP 2: Darken the colors
def darken_color(color, factor=0.9):
    rgb = np.array(mcolors.to_rgb(color))
    return np.clip(rgb * factor, 0, 1)
# Apply darkening function
dark_tab20 = [darken_color(c, factor=1) for c in tab20_palette]
# Convert to HEX
dark_tab20_hex = [mcolors.to_hex(c) for c in dark_tab20]
# STEP 6: Ensure leiden column is string (if needed)
pcHVG_data.obs['leiden'] = pcHVG_data.obs['leiden'].astype(str)
# STEP 7: Plot UMAP with custom palette
sc.set_figure_params(figsize=(5,5), frameon=False, dpi=500)
sc.pl.umap(pcHVG_data, color=['leiden'],
           add_outline=False, alpha=0.5, size=30,
           frameon=False, legend_fontsize=7, legend_fontoutline=2,
           palette=dark_tab20_hex,
           save='leiden_cluster_names.pdf')


# Ploting % cell per cluster
temp1 = pd.DataFrame(pcHVG_data.obs['leiden'].value_counts())
temp1['cluster'] = temp1.index
temp1.columns = ['count', 'cluster']
temp1['percentage'] = temp1['count']/len(pcHVG_data)*100
temp1['cluster'] = temp1.index
plt.figure(figsize=(3, 3))
sns.set(style="ticks")
ax = sns.barplot(x='cluster', y="count",
                 data=temp1, linewidth=1, edgecolor='black',
                 palette='tab20')
# ax.bar_label(ax.containers[0]) # only 1 container needed unless using `hue`
plt.setp(ax.artists, edgecolor='k', facecolor='w')
ax.tick_params(axis='x', rotation=90)
ax.set(xlabel='Cell type', ylabel='Number of cells')
plt.setp(ax.lines, color='k')
sns.despine(offset=3, trim=False)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/figures/percemntage_clusters.pdf',  bbox_inches='tight')
plt.show()

# %%  RE ASSIGN CELL CLUSTERS TO NORM DATA
adata = sc.read_h5ad('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/adata_filtered_norm.h5ad')
adata.obs['leiden'] = pcHVG_data.obs['leiden']
# adata.var['highly_variable'] = pcHVG_data.var['highly_variable']
adata.uns['umap'] = pcHVG_data.uns['umap']
adata.obsm['X_pca'] = pcHVG_data.obsm['X_pca']
adata.obsm['X_umap'] = pcHVG_data.obsm['X_umap']
adata.uns['neighbors'] = pcHVG_data.uns['neighbors']
adata.uns['pca'] = pcHVG_data.uns['pca']
# adata.varm['PCs'] = pcHVG_data.varm['PCs']
adata.obsp['distances'] = pcHVG_data.obsp['distances']
adata.var['highly_variable'] =  False
adata.uns =  pcHVG_data.uns
adata.write_h5ad('/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/adata_filtered_norm_clusters.h5ad')

# %% OBTAIN LIST OF EXCLUSIVE DEG GENES PER CLUSTER 
import pandas as pd
# Step 1: Count how many clusters each gene appears in
gene_counts = markers_def['names'].value_counts()
# Step 2: Filter for genes appearing in exactly one group
exclusive_genes = gene_counts[gene_counts == 1].index
# Step 3: Filter the original dataframe to keep only those exclusive genes
exclusive_markers = markers_def[markers_def['names'].isin(exclusive_genes)]
exclusive_markers.to_excel("/Users/Sebastian/Documentos/SLCU_lab/Projects/scRNA-seq/repositories/moreno_etal_2024/data/marker_genes_reg_final_exclusive.xlsx")






