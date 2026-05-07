#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 17:39:26 2025

@author: Sebastian Moreno-Ramirez
"""
import pandas as pd 
import matplotlib.cm
import matplotlib.pyplot as plt  
import seaborn as sns
import scanpy as sc
import numpy as np
import os 
import matplotlib.colors as mcolors
from pathlib import Path
import matplotlib as mpl
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

def add_suffix_to_duplicates(lst):
    counts = {}
    result = []
    for item in lst:
        if item in counts:
            counts[item] += 1
            result.append(f"{item}_{counts[item]}")
        else:
            counts[item] = 0
            result.append(item)
    return result
   
def darken_cmap(cmap, factor=0.9):
    """Darkens a colormap by multiplying RGB values by a factor."""
    colors = cmap(np.linspace(0, 1, 256))  # Sample original colormap
    darkened_colors = np.clip(colors[:, :3] * factor, 0, 1)  # Darken RGB
    return mcolors.LinearSegmentedColormap.from_list("darkened_cmap", darkened_colors)

original_cmap = plt.get_cmap("BuPu")
darkened_cmap = darken_cmap(original_cmap, factor=1)

# %% Load data
adata = sc.read_h5ad(OUTPUT_DATA / 'adata_filtered_norm_clusters.h5ad') ## normalized counts single-nuclei
genes = pd.read_excel(OUTPUT_DATA / 'marker_genes_reg_final.xlsx', index_col= 0)

# %% Add GENE NAMES 
adata.raw.var['gene_ids'] = names_changes_list(adata.raw.var.index)

# %% epidermis 
sc.set_figure_params(frameon=False, dpi=300)
epidermis1 = [ 'ATML1', 'MYB4', 'ABCG13','HTH','YAB1', 'PDF1','HTH',  'YAB3' ] 
sc.pl.matrixplot(adata, var_names = epidermis1 , groupby= 'leiden', 
                 swap_axes = True , cmap = 'BuPu',
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "epidermis1_cluster_matrixplot.pdf")

#%% 
epidermis2 = ['CER1', 'CER3','MLP28', 'FAR3', 'ABCG18', 'JAL23', 'LTP7', 'CHS',] 
sc.pl.matrixplot(adata, var_names = epidermis2 , groupby= 'leiden', 
                 swap_axes = True , cmap = 'BuPu',
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "epidermis2_cluster_matrixplot.pdf")

# %%  Rib zone
rib_zone = ['PME53', 'PCO2',  'REM10',  'AED3', 'SPL4', 'GAMT2', 'AT1G12805', 'S-ACP-DES6', ] 
sc.pl.matrixplot(adata, var_names = rib_zone , groupby= 'leiden', 
                 swap_axes = True , cmap = 'BuPu',
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "ribzone_related_cluster_matrixplot.pdf")

# %%  Floral identity
floral_identity = ['PI',  'AP3', 'SUP', 'AG', 'SEP1', 'SEP2', 'SEP3'] 
sc.pl.matrixplot(adata, var_names = floral_identity ,
                  groupby= 'leiden', 
                  swap_axes = True , cmap = darkened_cmap, 
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "floral_identity_related_cluster_matrixplot.pdf")


# %%  vasculature
vasculature2 = ['TMO5', 'FAF1', 'FAF3', 'CLE46', 'LAMP1', 'GH3.6', 'XTH4', 'WRKY70'] 
sc.pl.matrixplot(adata, var_names = vasculature2 , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap, 
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "vasculature2_cluster_matrixplot.pdf")

# %% 
vasculature = ['ATHB8', 'SMXL5' , 'PFA5',   'BRXL4', 'TMO6', 'SHR', 'ACL5'] 
sc.pl.matrixplot(adata, var_names = vasculature , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap,
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black')

# %% 

iSC = [ 'WUS','CLV3', 'CLV1', 'AT1G26680', 'AT3G59270','CKX3', 'AHK4',  'LOG4', 'AIL7','MCT2'  ] 
sc.pl.matrixplot(adata, var_names = iSC , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap, 
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (9,2),
                num_categories = 5, 
                colorbar_title ='Standarized  \nmean expression  \nper cluster',
                linewidth = 0.5, color = 'black',
                save =  "iSC_related_cluster_matrixplot.pdf")


sc.set_figure_params(figsize=(4,4),frameon=False, dpi=500)
sc.pl.umap(adata, color=['AT2G17950', 'AT2G27250', 'AT1G26680', 'AT3G59270'], title = ['WUS', 'CLV3', 'AT1G26680', 'CKX3'],
            add_outline=False,alpha =1, size =  200,
        frameon=True, legend_fontsize=7, legend_fontoutline=5,
        cmap="magma", 
        save = "iSC_UMAP-related_cluster_matrixplot.pdf")

# %%  Obtain AUXIN related genes as marker genes ####
auxin = pd.read_excel(INPUT_DATA / 'auxin_related_genes_curated.xlsx')
genes = pd.read_excel(OUTPUT_DATA / 'marker_genes_reg_final.xlsx', index_col= 0)

auxin = list(set(auxin['AGI']) & set(genes.names.tolist()))  ## Auxin-related genes as marker genes 
# auxin = list(set(auxin) & set(genes['names'].tolist()))
auxin = names_changes_list(auxin)
sc.set_figure_params(figsize=(4,4),frameon=False, dpi=300)

sc.pl.matrixplot(adata, var_names = auxin,  groupby='leiden', #dendrogram=True,
                gene_symbols= 'gene_ids', swap_axes = True, figsize= (5,5),
                standard_scale='var', 
                use_raw = True,
                colorbar_title ='Standarized  \nmean expression  \nper cluster',
                cmap = darkened_cmap,
                linewidth = 0.5, color = 'black',
                )
# %%  auxin reordered 
auxin_reorderd = ['ARF6',
                  'ETT',
                 'IAA30',
                 'IAA20',
                 'PIN6',
                 'GH3.3',
                 'GH3.5',
                 'AUX1',
                 'LAX1',
                 'GH3.2',
                 'PIN7',
                 'ARF11',         
                 'YUC1',
                 'YUC4',
                 'ARF5',
                 'ARF19',
                 'PIN3',
                 'IAA7',
                 'PIN4',
                 'LAX2',
                 'GRH1',
                 'ARF4',
                 'PIN1',
                 'GH3.6', 
                 ]
sc.set_figure_params(frameon=False, dpi=500)
sc.pl.matrixplot(adata, var_names = auxin_reorderd,  groupby='leiden', #dendrogram=True,
                gene_symbols= 'gene_ids', swap_axes = True, figsize= (5,6),
                standard_scale='var', 
                use_raw = True,
                colorbar_title ='Standarized  \nmean expression  \nper cluster',
                cmap = 'BuPu',
                linewidth = 0.5, color = 'black',
                )

# %% CK 
ck = pd.read_excel(INPUT_DATA / 'CK-related_genes_curated.xlsx', header= 0)
ck = list(set(ck['LOCUS']) & set(adata.var.index.tolist()))
ck = list(set(ck) & set(genes['names'].tolist()))
ck = names_changes_list(ck)
# CK_related = ['AHK4','CKX3','UGT76C2', 'ARR7','CKX6', 'CKX5',  'ARR1', 'AHP6',   'LOG5',  'ARR5', 'ARR10',  'ARR4', 'PUP14'] 

sc.pl.matrixplot(adata, var_names = ck ,
                  groupby= 'leiden', swap_axes = True , cmap = darkened_cmap,
                  figsize= (7,4),
                  use_raw = True,
                  gene_symbols= 'gene_ids', 
                standard_scale='var', 
                linewidth = 0.5, color = 'black', alpha = 1,
                # gene_symbols = 'gene_ids',
                colorbar_title ='Standardized  \ngene expression',
                save =  "CK_related_matrixplot.pdf",
                )

# %% 
ga = pd.read_excel(INPUT_DATA / 'GA-related_genes_curated.xlsx', header= 0)
ga = list(set(ga['LOCUS']) & set(adata.var.index.tolist()))
ga = list(set(ga) & set(genes['names'].tolist()))
ga = names_changes_list(ga)

sc.set_figure_params(frameon=False, dpi=300)
GA = ga 
sc.pl.matrixplot(adata, var_names = GA, 
                  groupby= 'leiden', swap_axes = True , cmap = 'BuPu',
                  figsize= (7,1.5),
                  gene_symbols= 'gene_ids', 
                standard_scale='var', 
                linewidth = 0.5, color = 'black', alpha = 1,
                colorbar_title ='Standardized  \ngene expression',
                )


# %% mesophyll
mesophyll = ['LHCB3', 'LHCB2.2','GER3', 'LHCB1.1', 'GLP1', 'LHCB5'] 
sc.pl.matrixplot(adata, var_names = mesophyll , groupby= 'leiden', 
                 swap_axes = True , cmap = 'BuPu',
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,  1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "myrosin_related_cluster_matrixplot.pdf")

# %% 
s_phase = ['HIS4',  'HTB1', 'HTB11', 'HTA5', 'HTR11', 'HTA13', 'HTA10'] 
sc.pl.matrixplot(adata, var_names = s_phase , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap,
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "s-phase_related_cluster_matrixplot.pdf")

# %% 
m_phase = ['CDC20-1',  'SCL28', 'FZR3', 'KIN7A', 'CYCB2-2', 'KIN12D', 'IMK2', 'TOP2'] 
sc.pl.matrixplot(adata, var_names = m_phase , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap, 
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "m-phase-phase_related_cluster_matrixplot.pdf")


# %% 
pollen = ['AMS',  'TA1', 'ABCG26', 'CYP704B1', 'TKPR1', 'SWEET3', 'PKSB', 'MYB99'] 
sc.pl.matrixplot(adata, var_names = pollen , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap, 
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "pollen_related_cluster_matrixplot.pdf")


# %% 
cell_wall = ['CLE41', 'PME5',  'PME48', 'AGP6', 'RALFL8', 'RALFL9', 'ATA20', 'PGA3', 'GRP19'] 
sc.pl.matrixplot(adata, var_names = cell_wall , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap, 
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black')
                # save =  "pollen_related_cluster_matrixplot.pdf")


# %%  CORTEX
cortex = ['JKD',   'AED3', 'EXT3', 'VSP2', 'MT2A', 'MGP','C/VIF2'] 
sc.pl.matrixplot(adata, var_names = cortex ,linewidth = 0.5, color = 'black',
                 groupby= 'leiden', swap_axes = True , cmap = darkened_cmap,
                 figsize= (7,1.5),
                standard_scale='var', 
                gene_symbols = 'gene_ids',
                colorbar_title ='Standardized  \ngene expression',
                save =  "cortex_matrixplot.pdf")


# %% BD
boundary = ['CUC3', 'CUC2', 'WOX12', 'NPF1.1', 'SOD7', 'NAC041', 'NPF1.2'] 
# sc.pl.umap(pcHVG_data, color=boundary, #legend_loc='on data',
#            add_outline=True, cmap = 'magma', size = 100, 
#        frameon=False, legend_fontsize=7, legend_fontoutline=2, 
#        palette="Set1", save= "boundary" + ".pdf",
#        title=[ 'clusters', 'SOD7', 'CUC2', 'CUC3'  , 'WOX12'])
sc.pl.matrixplot(adata, var_names = boundary , groupby= 'leiden', 
                 swap_axes = True , cmap = darkened_cmap,
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black')
# %%  EP 
early_primordia = ['AHP6', 'LAX1', 'TPPJ','DRNL', 'DOF5.8', 'IAA20', 'SAP', 'PLT3'] 
sc.pl.matrixplot(adata, var_names = early_primordia , groupby= 'leiden', 
                 swap_axes = True , cmap =  darkened_cmap,
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7,1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "EP_cluster_matrixplot.pdf")

# %% 
stress_response = ['ERF109', 'WRKY40', 'LOX3','LOX4', 'ERF018', 'XTH22', 'BAP1'] 
sc.pl.matrixplot(adata, var_names = stress_response , groupby= 'leiden', 
                 swap_axes = True , cmap =  darkened_cmap,
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (7, 1.5),
                colorbar_title ='Standardized  \ngene expression',
                linewidth = 0.5, color = 'black',
                save =  "stress_response_cluster_matrixplot.pdf")

# %%  hypoxia
hypoxia = ['AT5G65510',   'AT5G06300', 'AT5G07930', 'AT3G59270', 'AT5G56970', 'AT1G26680'] 
sc.pl.matrixplot(adata, var_names = hypoxia , groupby= 'leiden', 
                 swap_axes = True , cmap =  'BuPu',
                standard_scale='var', 
                gene_symbols= 'gene_ids', 
                figsize= (8,2),
                colorbar_title ='Standarized  \nmean expression  \nper cluster',
                linewidth = 0.5, color = 'black',
                save =  "hypoxia_cluster_matrixplot.pdf")


# %%  PLOTTING UMAPS 

adata = sc.read_h5ad(OUTPUT_DATA / 'adata_filtered_norm_clusters.h5ad') ## normalized counts single-nuclei

adata.var['gene_ids'] = names_changes_list(adata.var.index)

sc.set_figure_params(figsize=(4,4),frameon=False, dpi=300)
sc.pl.umap(adata, color=['SOD7', 'CUC2', 'CUC3' , 'NAC041', 'leiden'], 
           # title = ['WUS', 'ABI5', 'ABI2', 'ABI1'],
           size =  100,  gene_symbols= 'gene_ids', 
       add_outline=False , 
       frameon=True, legend_fontsize=7, legend_fontoutline=5,
       cmap="inferno")

sc.set_figure_params(figsize=(4,4),frameon=False, dpi=300)
sc.pl.umap(adata, color=['ATHB-8', 'PFA5', 'PFA3' , 'DOF2.5'], 
           # title = ['WUS', 'ABI5', 'ABI2', 'ABI1'],
           size =  100,  gene_symbols= 'gene_ids', 
       add_outline=False , 
       frameon=True, legend_fontsize=7, legend_fontoutline=5,
       cmap="inferno")

sc.set_figure_params(figsize=(4,4),frameon=False, dpi=300)
sc.pl.umap(adata, color=['AHP6', 'DRNL', 'TPPJ', 'IAA20'], 
           # title = ['WUS', 'ABI5', 'ABI2', 'ABI1'],
           size =  100,  gene_symbols= 'gene_ids', 
       add_outline=False , 
       frameon=True, legend_fontsize=7, legend_fontoutline=5,
       cmap="inferno")

sc.set_figure_params(figsize=(4,4),frameon=False, dpi=500)
sc.pl.umap(adata, color=['AHP6', 'GRXC7', 'DRNL', 'TPPJ'] , 
           # title = ['WUS', 'ABI5', 'ABI2', 'ABI1'],
           size =  100,  gene_symbols= 'gene_ids', 
       add_outline=False , 
       frameon=True, legend_fontsize=7, legend_fontoutline=5,
       cmap="inferno",
       save = 'EP_UMAP.pdf')

