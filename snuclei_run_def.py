#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:49:23 2022

@author: Sebastian
"""


import glob
import pandas as pd 
import os 
import matplotlib.pyplot as plt  
import seaborn as sns

import scanpy as sc
from sklearn.neighbors import NearestNeighbors 
# from sklearn.neighbors import DistanceMetric
from sklearn.decomposition import PCA
#from kneed import KneeLocator as kl 
from scipy import sparse
import numpy as np

import scrublet as scr

import matplotlib as mpl
import anndata as ad

# import scvelo as scv
# import scvi

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor = "white")


def scrublet_processing(adata):
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    scrub.call_doublets(threshold=0.25)
    scrub.plot_histogram()
    adata.obs["Doublet"] = predicted_doublets
    #list_doublets = list(adata.obs[adata.obs["Doublet"] == False].index)
    adatab = adata[adata.obs.Doublet == False, :] ## Removing doublets barcodes 
    adatab.write_h5ad('adata_scrublet.h5ad')
    number = len(adata.obs) - len(adatab.obs)
    print(str(number) + ' doublets detected by scrublet')
    return adatab

# def scvi_doublet_removal(adata, threshold_value ):
#     sc.pp.filter_genes(adata, min_cells = 5)
#     sc.pp.highly_variable_genes(adata, n_top_genes= 2000, subset = True, flavor = 'seurat_v3')
#     scvi.model.SCVI.setup_anndata(adata)
#     vae = scvi.model.SCVI(adata)
#     vae.train()
#     solo = scvi.external.SOLO.from_scvi_model(vae)
#     solo.train()
#     df = solo.predict()
#     df['prediction'] = solo.predict(soft = False)
#     df.index = df.index.map(lambda x: x[:-2])
#     df['dif'] = df.doublet - df.singlet
#     sns.distplot(df[df.prediction == 'doublet'], x = df[df.prediction == 'doublet']['dif'])
#     doublets = df[(df.prediction == 'doublet') & (df.dif > threshold_value)]
#     return doublets 


def qc_visualization(adata):  
    adata.var["atc"] = adata.var_names.str.startswith("ATC")
    sc.pp.calculate_qc_metrics(adata,qc_vars=["atc"], percent_top = None, log1p=False, inplace=True)
    adata.var["atm"] = adata.var_names.str.startswith("ATM")
    ath_ribo = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/ribosome_units_ath.xlsx')
    ath_ribo = ath_ribo[ath_ribo['CHROMOSOME NO.'].str.contains('AT', na= False)]
    ath_ribo['CHROMOSOME NO.']  = ath_ribo['CHROMOSOME NO.'].apply(lambda x: str(x).replace(u'\xa0', u''))
    adata.var['ribo'] = adata.var_names.isin(list(ath_ribo['CHROMOSOME NO.'])) 
    sc.pp.calculate_qc_metrics(adata,qc_vars=["atm", 'ribo'], percent_top = None, log1p=False, inplace=True)
    ax = sc.pl.scatter(adata,x = "total_counts", y = "pct_counts_atc")
    ax = sc.pl.scatter(adata,x = "total_counts", y = "pct_counts_atm")
    ax = sc.pl.scatter(adata,x = "total_counts", y = "n_genes_by_counts")
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_atm"], 
                 jitter=  4, multi_panel=True)
    # ax = sns.histplot(adata.obs["total_counts"])
    return adata
   
def qc_processing(adata1, min_genes, min_cells, max_genes, max_reads, min_reads, atm ,atc):   
    sc.pp.filter_cells(adata1, min_genes= min_genes)
    sc.pp.filter_genes(adata1, min_cells = min_cells)  ## Genes not consistant 
    sc.pp.filter_cells(adata1, max_genes=max_genes)
    adata1 = adata1[adata1.obs.total_counts < max_reads, :] ## Cells with less than 50000 reads
    adata1 = adata1[adata1.obs.total_counts > min_reads, :]  ##cells more than 1000 reads 
    adata1 = adata1[adata1.obs.pct_counts_atc < atc, :] ## Chloroplast
    adata1 = adata1[adata1.obs.pct_counts_atm < atm, :] ## Mitochondria 
    ax = sc.pl.scatter(adata1,x = "total_counts", y = "n_genes_by_counts")
    return adata1



def hvg_pca(adata):
    ##HCV: genes are binned by their mean expression, and the genes with the highest variance-to-mean- ratio are selected as HVGs in each bin. 
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) ## default values 
    HVG_data = adata.copy()
    HVG_data.raw = HVG_data
    HVG_data = HVG_data[:,HVG_data.var.highly_variable]
    sc.tl.pca(HVG_data,svd_solver="arpack" ) ### COMPUTING PCA 
    sc.pl.pca_variance_ratio(HVG_data) ## To decide how many PCA compute
    return HVG_data


def umap(adata, resolution,n_pcs, n_neighbors, method, output):
    ## reduce the dimensionaility 
    sc.pp.neighbors(adata, n_pcs = n_pcs , n_neighbors = n_neighbors) #Compute a neighborhood graph of observations
    # visualizing data with umap
    sc.tl.umap(adata) 
    sc.tl.leiden(adata, resolution = resolution)
    sc.tl.rank_genes_groups(adata, 'leiden', method = method)
    sc.pl.rank_genes_groups(adata, n_genes = 10, sharey = False, 
                                save = "rank_groups:" + str(resolution)+ ".pdf")
    output_rank = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(10)
    output_rank = names_change_dataframe(output_rank)
    cluster2annotation = {} 
    for col in list(output_rank.columns): 
        temp1 = output_rank[col]
        temp2 = [x for x in list(temp1) if str(x) != 'nan']
        temp3 = str(temp2)[1:-1] 
        cluster2annotation[str(col)] = temp3 
    name = 'rank_genes_groups_wilcoxon. res:'+ str(resolution)
    adata.obs[name] = adata.obs['leiden'].map(cluster2annotation).astype('category')
    sc.pl.umap(adata, color=name, #legend_loc='on data',
                  add_outline=True,
                  frameon=False, legend_fontsize=7, legend_fontoutline=2,
                  palette="tab20", save= str(name) + ".pdf")
    output_rank.to_excel(output + "output_rank_genes_group_10" +str(resolution)+ ".xlsx")    



def names_change_dataframe(dataframe): 
    pd_temp = dataframe
    names = pd.read_csv("data/gene_names_for_scrna.csv", sep = ",", index_col = 0)
    ab = dataframe.values.tolist()
    ab = list((set([y for x in ab for y in x])))
    cleanedList = [x for x in ab if str(x) != 'nan']
    for g in cleanedList: 
        print(g)
        temp1 = names[names["gene_ids"] == g ]
        if len(temp1) > 0:
            nn = temp1.iloc[0,1]
            print(nn)
            pd_temp = pd_temp.replace(g ,nn)
    return pd_temp

def regressing_out_cell_cycle(adata): 
    ## Alternative pathway, remove cell-cycle effect before clustering 
    cc = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/SITTC6/attempts/cell_cycle_df_phase_curated.xlsx")
    cc = cc[['AGI', 'name' ,'phase']]
    g2m = cc[cc['phase'] == 'G2/M']
    g2m = g2m.drop_duplicates()
    g2m = g2m.reset_index(drop=True)
    #g2m = names_change_dataframe(g2m)
    g2m_genes = list(g2m.AGI.unique())

    s = cc[cc['phase'] == 'G2/M']
    s = s.drop_duplicates()
    s = s.reset_index(drop=True)
    #s = names_change_dataframe(s)
    s_genes = list(s.AGI.unique())
    

    cell_cycle_genes = list(set(s_genes +g2m_genes))

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes= g2m_genes)

    return adata, s_genes, g2m_genes, cell_cycle_genes
## scale to unit variance of PCA , clip alues exceding var 10. 
#sc.pp.scale(pcHVG_data) # "Scaling has the effect that all genes are weighted equally for downstream analysis. There is currently no consensus on whether or not to perform normalization over genes. While the popular Seurat tutorials (Butler et al, 2018) generally apply gene scaling, the authors of the Slingshot method opt against scaling over genes in their tutorial 

def population_per_domain(adata_normalized, output): 
    cell_population_perc = pd.DataFrame()
    cell_types = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/cell_types_SAM.xlsx")
    for ct in cell_types.cell_type.unique().tolist():
        temp1 = cell_types[cell_types['cell_type']== ct]
        barcodes_ct = []
        for g in temp1.AGI.unique().tolist():
            temp2 = adata_normalized[adata_normalized[: , g].X > 0, :].obs.index.tolist()
            barcodes_ct = barcodes_ct + temp2
        barcodes_ct   = list(set(barcodes_ct))
        d = {'cell_type': ct, 'percentage': len(barcodes_ct) / len(adata_normalized) * 100 ,
             'genes' :str('|'.join(temp1.gene_name.tolist())), 
             'number_cells' :len(barcodes_ct) }
        df = pd.DataFrame(d, index = [0])
        cell_population_perc = cell_population_perc.append(df)
    cell_population_perc = cell_population_perc.reset_index()
    plt.figure(figsize=(3,3))
    sns.set(style="ticks")
    ax = sns.barplot(x="cell_type", y="percentage", 
                            data=cell_population_perc,
                           palette = 'tab20')
    for index, row in cell_population_perc.iterrows():
        ax.text(row.name, row.percentage +1 , row.number_cells, color='black', ha='center')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    ax.tick_params(axis='x', rotation=90)
    ax.set(xlabel='Cell domains', ylabel='Percentage (%)')
    plt.setp(ax.lines, color='k') 
    sns.despine(offset=3, trim=False)
    plt.savefig(output + "cell_domains_percentage"+  ".pdf",  bbox_inches='tight')
    plt.show()    
    return cell_population_perc
    
        
def separate_floral_meristem(adata_normalized): 
    cell_types = pd.read_excel("/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/nuclei/sequencing/cell_types_SAM.xlsx")
    fm = cell_types[cell_types['cell_type']=='floral_meristem']
    barcodes_ct = []
    for g in fm.AGI.unique().tolist():
        temp2 = adata_normalized[adata_normalized[: , g].X > 0, :].obs.index.tolist()
        barcodes_ct = barcodes_ct + temp2
    barcodes_ct   = list(set(barcodes_ct))
    adata_fm = adata_normalized[adata_normalized.obs.index.isin(barcodes_ct)]
    adata_im = adata_normalized[~adata_normalized.obs.index.isin(barcodes_ct)]
    return adata_fm, adata_im
 

def removing_outlier_clusters(adata, cluster_number): 
    barcode_cluster = adata[adata.obs['leiden'] ==cluster_number].obs.index.tolist()
    return barcode_cluster
         
         
def change_names_adata(adata1): 
    # gene_desc = pd.read_csv('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/gene_desc_location.csv', sep = ',')
    gene_desc = pd.read_csv("data/gene_names_for_scrna.csv", sep = ",", index_col = 0)

    # gene_desc['Gene name'] = gene_desc['Gene name'].fillna('unknown')
    new_gene_names_in_correct_order = []
    for g in list(adata1.var_names):
        temp1 = gene_desc[gene_desc['gene_ids'] == g]
        if len(temp1) == 0:
            new_gene_names_in_correct_order.append(g)  
            print(g, 'no match')
        if len(temp1) != 0:
            if list(temp1['new_names'])[0]== g:
                new_gene_names_in_correct_order.append(g) 
                print(g, 'unknown')
            if list(temp1['new_names'])[0] != g:
                new_gene_names_in_correct_order.append(list(temp1['new_names'])[0])
                print(g, list(temp1['new_names'])[0])
    seen = set()
    result = []
    suffix = 1
    for name in new_gene_names_in_correct_order:
        while name in seen:
            name = name + str(suffix)
            suffix += 1
        seen.add(name)
        result.append(name)
        suffix = 1
    return  result    

def scar_ambient_denoising(adata, adata_raw): 
    from scar import setup_anndata, model
    setup_anndata(
        adata = adata,
        raw_adata = adata_raw,
        prob = 0.999,
        kneeplot = True
        )
    ambient_profile0 = adata.uns["ambient_profile_Gene Expression"]  ##ambient profile per gene 
    scarObj = model(raw_count = adata.to_df(),
                      ambient_profile = ambient_profile0, # In the case of default None, the ambient_profile will be calculated by averaging pooled cells
                      feature_type='mRNA',
                      sparsity=1 )  # initialize scar model
    scarObj.train(epochs=400,
                    batch_size=64,
                    verbose=True)  # start training
    scarObj.inference(batch_size = None)  # inference
    denoised_count = pd.DataFrame(scarObj.native_counts, index=adata.to_df().index, columns=adata.to_df().columns)
    return denoised_count

    
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


    
    

def coexpression_cluster(adata_normalized, adata_with_clusters, cluster, top_expressed, bon_pvalue):
    from scipy import stats
    adata_normalized.obs['cluster'] = adata_with_clusters.obs['leiden'].tolist()  ## Add cluster name to normalized dataset with barcodes generated with HVG. 
    # adata3 = change_names_adata(adata3)
    # data = adata_normalized.X.toarray()
    out = []
    adata_tmp1 = adata_normalized[adata_normalized.obs['cluster'] == cluster]  
    data = pd.DataFrame(adata_tmp1.X.toarray(), columns = adata_tmp1.var_names)
    mean = data.mean()
    genes = mean.sort_values(ascending = False)[0:top_expressed].index.tolist() 
    data = data[[c for c in data.columns if c in genes]]
    data.index = adata_tmp1.obs.index
    for g1 in list(range(len(data.columns))): 
        for g2 in  list(range(len(data.columns))):
            print(str(g1) + ' from '  + str(len(data.columns)))
            res = stats.pearsonr(data.iloc[:,g1], data.iloc[:,g2])
            out.append([data.iloc[:,g1].name,data.iloc[:,g2].name,res[0], res[1], cluster])
    df = pd.DataFrame(out, columns = ['gene1', 'gene2','r','p','cluster' ])       
    df['bon'] = df.p * len(df)
    df['-log10_p'] = -np.log10(df.p)
    df_filter = df[df.bon < bon_pvalue].sort_values('bon').reset_index(drop = True)
    return df_filter


def coexpression_arboreto(adata_normalized, algorithm ):
    from arboreto.algo import grnboost2, genie3
    # from arboreto.utils import load_tf_names    
    data = pd.DataFrame(adata_normalized.X, columns = adata_normalized.var_names, index = adata_normalized.obs.index )
    tf = pd.read_csv("data/Ath_TF_list.txt", sep = "\t")
    tf = tf['Gene_ID']
    if algorithm == 'grnboost2':
        network = grnboost2(expression_data=data,
                            tf_names=tf.tolist())
    if algorithm == 'genie3':
        network = genie3(expression_data=data,
                            tf_names=tf.tolist())
    return network


def coexpression_arboreto_by_cluster(adata_normalized_with_cluster, algorithm,  cluster):
    from arboreto.algo import grnboost2, genie3
    # g_cluster = dataframe  # from arboreto.utils import load_tf_names 
    # g_cluster = g_cluster[cluster].tolist()
    # # adata_with_cluster.obs['cluster'] = adata_with_cluster.obs['cluster']
    # adata_normalized.obs['cluster']=  adata_with_cluster.obs['cluster']
    # g_cluster = 
    atemp1 = adata_normalized_with_cluster[adata_normalized_with_cluster.obs['cluster'] == cluster]
    # non_mito_genes_list = [name for name in adata.var_names if not name.startswith('MT-')]
    # atemp2 = atemp1[:, g_cluster]
    data = pd.DataFrame(atemp1.X.toarray(), columns = atemp1.var_names, index = atemp1.obs.index )
    tf = pd.read_csv("data/Ath_TF_list.txt", sep = "\t")
    tf = tf['Gene_ID']
    if algorithm == 'grnboost2':
          network = grnboost2(expression_data=data,
                           tf_names=tf.tolist())     
    if algorithm == 'genie3':
        network = genie3(expression_data=data,
                         tf_names=tf.tolist())
    network['cluster'] = cluster
    return network


def umap_3d(anndata, path): 
    sc.tl.umap(anndata, n_components = 3)
    color = dict(zip(range(0,20), plt.cm.tab20(range(0,20))))
    umap = anndata.obsm['X_umap']

    for i in range(0,360,4):
        fig = plt.figure(figsize = (10,10))
        ax = fig.add_subplot(projection = '3d')
        ax.scatter(umap[:,0], umap[:,1], umap[:,2], c = anndata.obs.leiden.astype('int').map(color))
        x_center = (umap[:,0].max() + umap[:,0].min())/2	
        y_center = (umap[:,1].max() + umap[:,1].min())/2
        z_center = (umap[:,2].max() + umap[:,2].min())/2
        
        
        ax.plot([x_center,x_center],[y_center,y_center], [umap[:,2].min() -0.5, umap[:2].max() +0.5], c = 'k' , lw = 2 )
        ax.plot([x_center,x_center], [umap[:,1].min() -0.5, umap[:,1].max() +0.5], [z_center, z_center], c = 'k' , lw = 2 )
        ax.plot([umap[:,0].min() -0.5, umap[:,0].max() +0.5], [y_center,y_center],[z_center, z_center], c = 'k' , lw = 2 )
        ax.view_init(20,i)
        
        ax.axis('off')
        plt.savefig(path +  str(i) , dpi = 100, facecolor = 'white')
        plt.show()
        from natsort import natsorted
        import imageio
        images = []
        for file_name in natsorted(os.listdir(path)):
            if file_name.endswith('.png'):
                file_path = os.path.join(path, file_name)
                images.append(imageio.imread(file_path))
        imageio.mimsave('animation_3d_umap.gif', images, duration=0.5)



def cell_per_cluster(adata, output): 
    temp1 = pd.DataFrame(adata.obs['leiden'].value_counts())
    temp1['cluster'] = temp1.index
    temp1['percentage'] = temp1['count']/len(adata)*100
    temp1['cluster'] = temp1.index
    plt.figure(figsize=(3,3))
    sns.set(style="ticks")
    ax = sns.barplot(x='cluster', y="percentage", 
                            data=temp1,
                           palette = 'tab20')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    ax.tick_params(axis='x', rotation=90)
    ax.set(xlabel='Cell type', ylabel='Percentage (%)')
    plt.setp(ax.lines, color='k') 
    sns.despine(offset=3, trim=False)
    plt.savefig(output + "cell_domains_percentage"+  ".pdf",  bbox_inches='tight')
    plt.show()    
    return temp1

    
    
    
    
    






 